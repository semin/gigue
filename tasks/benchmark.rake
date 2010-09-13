require 'fork_manager'

include Gigue

$logger.level = Logger::INFO

max_fork  = 6
cur_dir   = Pathname.new(__FILE__).dirname.realpath
tmp_dir   = cur_dir + "../tmp"
ali_dir   = tmp_dir + "ali"
data_dir  = cur_dir + "../data"
mat_dir   = data_dir + "mats"
list_dir  = mat_dir + "lists"
dna_tems  = Pathname::glob(data_dir + "./bipa/scop/rep/dna/*/dnamodsalign.tem")
rna_tems  = Pathname::glob(data_dir + "./bipa/scop/rep/rna/*/rnamodsalign.tem")
envs      = [64, 128]
pids      = (40..90).step(10)
sigmas    = [0.0001, 0.001, 0.01, 0.1, 3, 5]

# load SCOP description file
sunid_to_sccs = {}
sid_to_sccs   = {}
sid_to_sunid  = {}
scop_des_file = data_dir + "scop/dir.des.scop.txt_1.75"

IO.foreach(scop_des_file) do |line|
  line.chomp!
  if line =~ /^#/
    next
  elsif (cols = line.split("\t")).size == 5
    sunid_to_sccs[cols[0]]  = cols[2]
    sid_to_sccs[cols[3]]    = cols[2] if cols[1] == "px"
    sid_to_sunid[cols[3]]   = cols[0] if cols[1] == "px"
  else
    warn "?: #{line}"
  end
end

def match_positions(aas1, aas2)
  len1 = aas1.length
  len2 = aas2.length
  gap1 = 0
  gap2 = 0
  idx1 = 0
  idx2 = 0
  aafnd1 = false
  aafns2 = false
  cols = []

  for i in (0...len1)
    if (aas1[i] == '-') && (aas2[i] == '-')
      gap1 += 1
      gap2 += 1
    elsif (aas1[i] == '-') && (aas2[i] != '-')
      gap1 += 1
    elsif (aas1[i] != '-') && (aas2[i] == '-')
      gap2 += 1
    elsif (aas1[i] != '-') && (aas2[i] != '-')
      cols << [i-gap1, i-gap2]
    else
      warn "Unknown combination!"
      exit 1
    end
  end
  cols
end


namespace :bench do
  namespace :ali do
    desc "Test alignment accuracy of Needleman-Wunsch algorithm using BLOSUM62 for DNA/RNA-bindnig protein families"
    task :nw do
      fm = ForkManager.new(max_fork)
      fm.manage do
        %w[dna rna].each do |na|
          tems = Pathname::glob(data_dir + "./bipa/scop/rep/#{na}/*/#{na}modsalign.tem")
          tems.each do |tem|
            fm.fork do
              fam_sunid     = tem.to_s.match(/(\d+)/)[1]
              fam_sccs      = sunid_to_sccs[fam_sunid]
              bio_tem       = Bio::FlatFile.auto(tem)
              entries       = bio_tem.entries
              entries_hash  = {}
              csv_file      = ali_dir + "nw/#{na}/#{fam_sccs}.csv"
              csv_file.open("w") do |file|
                entries.each do |entry|
                  if !entries_hash.has_key?(entry.entry_id)
                    entries_hash[entry.entry_id] = {}
                  end
                  entries_hash[entry.entry_id][entry.definition] = entry
                end

                entries_hash.keys.each do |ref_entry_id|
                  ref     = entries_hash[ref_entry_id]
                  others  = entries_hash.reject { |k, v| k == ref_entry_id }
                  others.each do |qry_entry_id, qry_entry_hash|
                    ref_seq_ent       = ref['sequence']
                    ref_qry_ent       = qry_entry_hash['sequence']
                    ref_str_seq       = Sequence.new(ref_seq_ent.aaseq, ref_seq_ent.entry_id, ref_seq_ent.definition)
                    ref_qry_seq       = Sequence.new(ref_qry_ent.aaseq, ref_qry_ent.entry_id, ref_qry_ent.definition)
                    ref_seq_seq_cols  = match_positions(ref_seq_ent.aaseq.split(''), ref_qry_seq.amino_acids)
                    ref_seq_seq_pid   = Sequence::calculate_pid_cpp(ref_seq_ent.aaseq.split(''), ref_qry_seq.amino_acids)

                    raw_ref_seq       = Sequence.new(ref_seq_ent.aaseq.gsub('-',''), ref_seq_ent.entry_id, ref_seq_ent.definition)
                    raw_qry_seq       = Sequence.new(ref_qry_ent.aaseq.gsub('-',''), ref_qry_ent.entry_id, ref_qry_ent.definition)
                    ssa               = SequenceSequenceAligner.new(raw_ref_seq, raw_qry_seq)
                    gaf               = ssa.global_alignment_affine_gap_cpp
                    gig_seq_seq_cols  = match_positions(gaf.aligned_amino_acids1, gaf.aligned_amino_acids2)
                    gig_seq_seq_pid   = Sequence::calculate_pid_cpp(gaf.aligned_amino_acids1, gaf.aligned_amino_acids2)

                    mat_cols  = (ref_seq_seq_cols & gig_seq_seq_cols).size
                    ali_score = mat_cols / Float(ref_seq_seq_cols.size)

                    elems = [
                      "#{ref_entry_id} (#{sunid_to_sccs[ref_entry_id]})",
                      "#{qry_entry_id} (#{sunid_to_sccs[qry_entry_id]})",
                      ref_seq_seq_pid,
                      gig_seq_seq_pid,
                      ali_score
                    ]

                    file.puts elems.join(", ")
                    $logger.info "%-20s %-20s %8.6f %8.6f %10.8f" % elems
                  end
                end
              end
            end
          end
        end
      end
    end
  end
end


namespace :bench do
  namespace :ali do
    desc "Test alignment accuracy of GIGUE using DNA/RNA-ESSTs, multipe PID cutoffs and sigma values"
    task :gig do

      fm = ForkManager.new(max_fork)
      fm.manage do
        %w[dna rna].each do |na|
          tems = Pathname::glob(data_dir + "./bipa/scop/rep/#{na}/*/#{na}modsalign.tem")
          tems.each do |tem|
            fm.fork do
              fam_sunid     = tem.to_s.match(/(\d+)/)[1]
              fam_sccs      = sunid_to_sccs[fam_sunid]
              bio_tem       = Bio::FlatFile.auto(tem)
              entries       = bio_tem.entries
              entries_hash  = {}
              csv_file      = ali_dir + "gig/#{na}/#{fam_sccs}.csv"
              csv_file.open("w") do |file|
                entries.each do |entry|
                  if !entries_hash.has_key?(entry.entry_id)
                    entries_hash[entry.entry_id] = {}
                  end
                  entries_hash[entry.entry_id][entry.definition] = entry
                end

                entries_hash.keys.each do |ref_entry_id|
                  ref = entries_hash[ref_entry_id]
                  joy = ""
                  ref.values.each do |e|
                    joy += %Q{>P1;#{e.entry_id}\n#{e.definition}\n#{e.data.gsub('-','').gsub("\n",'')}*\n}
                  end

                  others = entries_hash.reject { |k, v| k == ref_entry_id }
                  others.each do |qry_entry_id, qry_entry_hash|
                    envs.each do |env|
                      pids.each do |pid|
                        sigmas.each do |sigma|
                          esst              = mat_dir + "ulla-#{na}#{env}-pid#{pid}-sigma#{sigma}-out2.mat"
                          str_prf           = StructuralProfile.create_from_joy_tem_and_essts(joy, esst)
                          ref_seq_ent       = ref['sequence']
                          ref_qry_ent       = qry_entry_hash['sequence']
                          ref_str_seq       = Sequence.new(ref_seq_ent.aaseq, ref_seq_ent.entry_id, ref_seq_ent.definition)
                          ref_qry_seq       = Sequence.new(ref_qry_ent.aaseq, ref_qry_ent.entry_id, ref_qry_ent.definition)
                          ref_str_seq_cols  = match_positions(ref_seq_ent.aaseq.split(''), ref_qry_seq.amino_acids)
                          ref_str_seq_pid   = Sequence::calculate_pid_cpp(ref_seq_ent.aaseq.split(''), ref_qry_seq.amino_acids)

                          raw_qry_seq       = Sequence.new(ref_qry_ent.aaseq.gsub('-',''), ref_qry_ent.entry_id, ref_qry_ent.definition)
                          psa               = ProfileSequenceAligner.new(str_prf, raw_qry_seq)
                          gal               = psa.global_alignment_affine_gap_cpp
                          gig_str_seq_cols  = match_positions(gal.aligned_structural_profile_positions.map(&:probe), gal.aligned_amino_acids)
                          gig_str_seq_pid   = Sequence::calculate_pid_cpp(gal.aligned_structural_profile_positions.map(&:probe), gal.aligned_amino_acids)
                          gig_prf_seq       = Sequence.new(gal.aligned_structural_profile_positions.map(&:probe).join(''), ref_entry_id, ref_seq_ent.definition)
                          gig_qry_seq       = Sequence.new(gal.aligned_amino_acids.join(''), qry_entry_hash, ref_qry_ent.definition)

                          mat_cols  = (ref_str_seq_cols & gig_str_seq_cols).size
                          ali_score = mat_cols / Float(ref_str_seq_cols.size)

                          elems = [
                            "#{ref_entry_id} (#{sunid_to_sccs[ref_entry_id]})",
                            "#{qry_entry_id} (#{sunid_to_sccs[qry_entry_id]})",
                              File.basename(esst, ".mat").gsub("ulla-ent-dna-","").gsub("-logo",""),
                                ref_str_seq_pid,
                                gig_str_seq_pid,
                                ali_score
                          ]

                          file.puts elems.join(", ")
                          $logger.info "%-20s %-20s %-40s %8.6f %8.6f %10.8f" % elems
                        end
                      end
                    end
                  end
                end
              end
            end
          end
        end
      end
    end
  end
end


namespace :bench do
  namespace :rec do
    desc "Parse PSI-BLAST results and create CSV files"
    task :psires do

      fm = ForkManager.new(8)
      fm.manage do
        %w[dna rna].each do |na|
          psi_dir = cur_dir + "../tmp/rec/psi/#{na}"
          xmls    = Pathname::glob(psi_dir + "*.xml")
          xmls.each do |xml|
            fm.fork do
              hits  = []
              rs    = Bio::Blast.reports(File.open(xml))
              rs.each do |r|
                r.each do |hit|
                  if (hit.target_id !~ /^UniRef/)
                    qry_sunid = hit.query_def
                    qry_sccs  = sunid_to_sccs[qry_sunid].split('.')
                    qry_class = qry_sccs[0]
                    qry_fold  = qry_sccs[0..1].join('.')
                    qry_sfam  = qry_sccs[0..2].join('.')
                    qry_fam   = qry_sccs[0..3].join('.')

                    tgt_sid   = hit.target_id.gsub(/^[g|e]/, 'd')
                    tgt_sunid = sid_to_sunid[tgt_sid]
                    tgt_sccs  = sid_to_sccs[tgt_sid].split('.')
                    tgt_class = tgt_sccs[0]
                    tgt_fold  = tgt_sccs[0..1].join('.')
                    tgt_sfam  = tgt_sccs[0..2].join('.')
                    tgt_fam   = tgt_sccs[0..3].join('.')

                    next if qry_class != tgt_class

                    cols = [
                      qry_sunid, qry_class, qry_fold, qry_sfam, qry_fam,
                      tgt_sunid, tgt_class, tgt_fold, tgt_sfam, tgt_fam,
                      hit.bit_score, hit.evalue
                    ]

                    log_fmt = "%8d %2s %5s %9s %12s %8d %2s %5s %9s %12s %10.3f %g"
                    $logger.info log_fmt % cols
                    hits << cols
                  end
                end
              end

              sorted_hits = hits.sort_by { |h| h[-1] }
              res_file    = psi_dir + "#{xml.basename('.xml')}.csv"
              res_file.open("w") do |file|
                sorted_hits.each do |sorted_hit|
                  file.puts sorted_hit.join(", ")
                end
              end
            end
          end
        end
      end
    end
  end
end


namespace :bench do
  namespace :rec do
    desc "Parse PSI-BLAST results and create CSV files"
    task :blares do

      fm = ForkManager.new(8)
      fm.manage do
        %w[dna rna].each do |na|
          bla_dir = cur_dir + "../tmp/rec/bla/#{na}"
          xmls    = Pathname::glob(bla_dir + "*.xml")
          xmls.each do |xml|
            fm.fork do
              hits  = []
              rs    = Bio::Blast.reports(File.open(xml))
              rs.each do |r|
                r.each do |hit|
                  if (hit.target_id !~ /^UniRef/)
                    qry_sunid = hit.query_def
                    qry_sccs  = sunid_to_sccs[qry_sunid].split('.')
                    qry_class = qry_sccs[0]
                    qry_fold  = qry_sccs[0..1].join('.')
                    qry_sfam  = qry_sccs[0..2].join('.')
                    qry_fam   = qry_sccs[0..3].join('.')

                    tgt_sid   = hit.target_id.gsub(/^[g|e]/, 'd')
                    tgt_sunid = sid_to_sunid[tgt_sid]
                    tgt_sccs  = sid_to_sccs[tgt_sid].split('.')
                    tgt_class = tgt_sccs[0]
                    tgt_fold  = tgt_sccs[0..1].join('.')
                    tgt_sfam  = tgt_sccs[0..2].join('.')
                    tgt_fam   = tgt_sccs[0..3].join('.')

                    next if qry_class != tgt_class

                    cols = [
                      qry_sunid, qry_class, qry_fold, qry_sfam, qry_fam,
                      tgt_sunid, tgt_class, tgt_fold, tgt_sfam, tgt_fam,
                      hit.bit_score, hit.evalue
                    ]

                    log_fmt = "%8d %2s %5s %9s %12s %8d %2s %5s %9s %12s %10.3f %g"
                    $logger.info log_fmt % cols
                    hits << cols
                  end
                end
              end

              sorted_hits = hits.sort_by { |h| h[-1] }
              res_file    = bla_dir + "#{xml.basename('.xml')}.csv"
              res_file.open("w") do |file|
                sorted_hits.each do |sorted_hit|
                  file.puts sorted_hit.join(", ")
                end
              end
            end
          end
        end
      end
    end
  end
end


namespace :bench do
  namespace :rec do
    desc "Test recognition performance of BLAST for DNA/RNA-binding protein families"
    task :bla do

      blastdb = cur_dir + "../data/scop/astral-scopdom-seqres-gd-sel-gs-bib-40-1.75.fa"

      %w[dna rna].each do |na|
        bla_dir = cur_dir + "../tmp/rec/bla/#{na}"

        fm = ForkManager.new(8)
        fm.manage do
          tems = Pathname::glob(data_dir + "./bipa/scop/rep/#{na}/*/#{na}modsalign.tem")
          tems.each_with_index do |tem, i|
            sunid = tem.to_s.match(/(\d+)/)[1]
            sccs  = sunid_to_sccs[sunid]
            unless sccs.match(/^[a|b|c|d|e|f|g]/)
              $logger.warn "Skip #{sccs}, it's not a true SCOP class"
              next
            end

            fm.fork do
              bio_tem = Bio::FlatFile.auto(tem)
              entries = bio_tem.entries
              entries.each do |e|
                if e.definition =~ /^sequence/
                  # 1. create an input file for sequences in a FASTA format
                  fasta       = ">#{e.entry_id}\n#{e.aaseq.gsub('-','')}"
                  fasta_file  = bla_dir + "#{e.entry_id}.fa"
                  fasta_file.open('w') { |f| f.puts fasta }

                  # 2. Run BLAST
                  blastout = bla_dir + "#{sccs}.xml"
                  cmd = "blastall -p blastp -i #{fasta_file} -e 1e10000 -d #{blastdb} -m 7 -o #{blastout}"
                  system cmd

                  $logger.info "BLAST searching against #{na.upcase}-bindnig SCOP family, #{sccs} (#{i+1}/#{tems.size}): done"
                end
              end
            end
          end
        end
      end
    end
  end
end


namespace :bench do
  namespace :rec do
    desc "Test recognition performance of PSI-BLAST for DNA/RNA-binding protein families"
    task :psi do

      #blastdb = cur_dir + "../data/scop/astral40-uniref90.fa"
      blastdb = cur_dir + "../data/scop/astral-scopdom-seqres-gd-sel-gs-bib-40-1.75.fa"
      #blastdb = cur_dir + "../data/scop/astral40.fa"

      %w[dna rna].each do |na|
        psi_dir = cur_dir + "../tmp/rec/psi/#{na}"

        fm = ForkManager.new(8)
        fm.manage do
          tems = Pathname::glob(data_dir + "./bipa/scop/rep/#{na}/*/#{na}modsalign.tem")
          tems.each_with_index do |tem, i|
            sunid = tem.to_s.match(/(\d+)/)[1]
            sccs  = sunid_to_sccs[sunid]
            unless sccs.match(/^[a|b|c|d|e|f|g]/)
              $logger.warn "Skip #{sccs}, it's not a true SCOP class"
              next
            end

            fm.fork do
              bio_tem       = Bio::FlatFile.auto(tem)
              entries       = bio_tem.entries
              entries_hash  = {}

              entries.each do |entry|
                if !entries_hash.has_key?(entry.entry_id)
                  entries_hash[entry.entry_id] = {}
                end
                entries_hash[entry.entry_id][entry.definition] = entry
              end

              # 1. create an input file for the first sequence in a FASTA format
              first_entry_id  = entries[0].entry_id
              first_entry_aa  = entries_hash[first_entry_id]['sequence']
              fasta           = ">#{first_entry_id}\n#{first_entry_aa.aaseq.gsub('-','')}"
              fasta_file      = psi_dir + "#{first_entry_id}.fa"
              fasta_file.open('w') { |f| f.puts fasta }

              # 2. create an input msa file in a ClustalW format
              clustalw = []
              entries_hash.keys.each do |entry_id|
                seq_entry = entries_hash[entry_id]['sequence']
                clustalw  << "%-10s     %s" % [entry_id, seq_entry.aaseq]
              end
              clustalw_file = psi_dir + "#{sccs}.msa"
              clustalw_file.open('w') { |f| f.puts clustalw.join("\n") }

              # 3. Run PSI-BLAST
              blastout = psi_dir + "#{sccs}.xml"
              cmd = "/BiO/Install/Blast/bin/blastpgp -i #{fasta_file} -B #{clustalw_file} -e 1e100 -t 1 -j 5 -d #{blastdb} -m 7 -o #{blastout}"
              system cmd

              $logger.info "PSI-BLAST searching against #{na.upcase}-bindnig SCOP family, #{sccs} (#{i+1}/#{tems.size}): done"
            end
          end
        end
      end
    end
  end
end


namespace :bench do
  namespace :rec do
    namespace :gig do
      desc "Test recognition performance of GIGUE with linear gap penalties for DNA-binding protein families"
      task :lindnasgl do

        seqs    = []
        scop40  = cur_dir + "../data/scop/astral-scopdom-seqres-gd-sel-gs-bib-40-1.75.fa"

        Bio::FlatFile.open(Bio::FastaFormat, scop40) do |ff|
          ff.each do |entry|
            if sid_to_sccs[entry.entry_id.gsub(/^[g|e]/, 'd')].match(/^[a|b|c|d|e|f|g]/)
              seqs << Sequence.new(entry.aaseq, entry.entry_id.gsub(/^[g|e]/, 'd'), "sequence;#{entry.definition}")
            else
              $logger.warn "Skip #{entry.definition}, not a true SCOP class"
            end
          end
        end

        fm = ForkManager.new(8)
        fm.manage do
          tems = Pathname::glob(data_dir + "./bipa/scop/rep/dna/*/dnamodsalign.tem")
          tems.each_with_index do |tem, i|
            sunid = tem.to_s.match(/(\d+)/)[1]
            sccs  = sunid_to_sccs[sunid]
            unless sccs.match(/^[a|b|c|d|e|f|g]/)
              $logger.warn "Skip #{sccs}, it's not a true SCOP class"
              next
            end
            $logger.info  "Searching against DNA-bindnig SCOP family, #{sccs} (#{i+1}/#{tems.size})"

            bio_tem = Bio::FlatFile.auto(tem)
            entries = bio_tem.entries
            entries_hash  = {}
            entries.each do |entry|
              if !entries_hash.has_key?(entry.entry_id)
                entries_hash[entry.entry_id] = {}
              end
              entries_hash[entry.entry_id][entry.definition] = entry
            end

            [ mat_dir + "ulla-dna64-pid60-sigma0.001-out2.mat",
              mat_dir + "ulla-dna128-pid60-sigma0.001-out2.mat"].each do |esst|
              fm.fork do

                hits  = []
                env   = esst.to_s.match(/(dna\d+)/)[1]

                entries_hash.keys.each do |entry_id|
                  ents = entries_hash[entry_id]
                  joy = ""
                  ents.values.each do |e|
                    joy += %Q{>P1;#{e.entry_id}\n#{e.definition}\n#{e.data.gsub('-','').gsub("\n",'')}*\n}
                  end

                  stp = StructuralProfile::create_from_joy_tem_and_essts(joy, esst)
                  seqs.each do |seq|
                    stp_sccs  = sccs.split('.')
                    stp_class = stp_sccs[0]
                    stp_fold  = stp_sccs[0..1].join('.')
                    stp_sfam  = stp_sccs[0..2].join('.')
                    stp_fam   = stp_sccs[0..3].join('.')
                    seq_sccs  = sid_to_sccs[seq.code].split('.')
                    seq_class = seq_sccs[0]
                    seq_fold  = seq_sccs[0..1].join('.')
                    seq_sfam  = seq_sccs[0..2].join('.')
                    seq_fam   = seq_sccs[0..3].join('.')

                    next if stp_class != seq_class

                    psa = ProfileSequenceAligner.new(stp, seq)
                    begin
                      gal = psa.global_alignment_linear_gap_cpp
                    rescue
                      $logger.debug "Problematic prf:"
                      $logger.debug ">#{entry_id}"
                      $logger.debug "Problematic seq:"
                      $logger.debug ">#{seq.code}"
                      $logger.debug "#{seq.data}"
                    end

                    if gal.raw_score < gal.reverse_score
                      $logger.debug "Skip #{sccs} <=> #{seq.code}, reverse raw score is greater than raw score."
                      next
                    end

                    result = [ entry_id, stp.length, stp_class, stp_fold, stp_sfam, stp_fam,
                      sid_to_sunid[seq.code], seq.code, seq.length, seq_class, seq_fold, seq_sfam, seq_fam,
                      gal.raw_score, gal.reverse_score, gal.z_score, psa.algorithm.to_s
                    ]

                    log_fmt = "%8d %5d %2s %5s %9s %12s %8d %-10s %5d %2s %5s %9s %12s %7d %7d %10.7f %10s"
                    $logger.info log_fmt % result
                    hits << result
                  end
                end

                sorted_hits = hits.sort_by { |h| h[-2] }.reverse
                res_file    = cur_dir + "../tmp/rec/gig/dna/single_str-single_seq-linear-#{sccs}-#{env}.csv"
                res_file.open("w") do |file|
                  sorted_hits.each do |sorted_hit|
                    file.puts sorted_hit.join(", ")
                  end
                end
              end
            end
          end
        end
      end
    end
  end
end


namespace :bench do
  namespace :rec do
    namespace :gig do
      desc "Test recognition performance of GIGUE with linear gap penalties for RNA-binding protein families"
      task :linrnasgl do

        seqs    = []
        scop40  = cur_dir + "../data/scop/astral-scopdom-seqres-gd-sel-gs-bib-40-1.75.fa"

        Bio::FlatFile.open(Bio::FastaFormat, scop40) do |ff|
          ff.each do |entry|
            if sid_to_sccs[entry.entry_id.gsub(/^[g|e]/, 'd')].match(/^[a|b|c|d|e|f|g]/)
              seqs << Sequence.new(entry.aaseq, entry.entry_id.gsub(/^[g|e]/, 'd'), "sequence;#{entry.definition}")
            else
              $logger.warn "Skip #{entry.definition}, not a true SCOP class"
            end
          end
        end

        fm = ForkManager.new(8)
        fm.manage do
          tems = Pathname::glob(data_dir + "./bipa/scop/rep/rna/*/rnamodsalign.tem")
          tems.each_with_index do |tem, i|
            sunid = tem.to_s.match(/(\d+)/)[1]
            sccs  = sunid_to_sccs[sunid]
            unless sccs.match(/^[a|b|c|d|e|f|g]/)
              $logger.warn "Skip #{sccs}, it's not a true SCOP class"
              next
            end
            $logger.info  "Searching against RNA-bindnig SCOP family, #{sccs} (#{i+1}/#{tems.size})"

            bio_tem = Bio::FlatFile.auto(tem)
            entries = bio_tem.entries
            entries_hash  = {}
            entries.each do |entry|
              if !entries_hash.has_key?(entry.entry_id)
                entries_hash[entry.entry_id] = {}
              end
              entries_hash[entry.entry_id][entry.definition] = entry
            end

            [ mat_dir + "ulla-rna64-pid60-sigma0.001-out2.mat",
              mat_dir + "ulla-rna128-pid60-sigma0.001-out2.mat"].each do |esst|
              fm.fork do

                hits  = []
                env   = esst.to_s.match(/(rna\d+)/)[1]

                entries_hash.keys.each do |entry_id|
                  ents = entries_hash[entry_id]
                  joy = ""
                  ents.values.each do |e|
                    joy += %Q{>P1;#{e.entry_id}\n#{e.definition}\n#{e.data.gsub('-','').gsub("\n",'')}*\n}
                  end

                  stp = StructuralProfile::create_from_joy_tem_and_essts(joy, esst)
                  seqs.each do |seq|
                    stp_sccs  = sccs.split('.')
                    stp_class = stp_sccs[0]
                    stp_fold  = stp_sccs[0..1].join('.')
                    stp_sfam  = stp_sccs[0..2].join('.')
                    stp_fam   = stp_sccs[0..3].join('.')
                    seq_sccs  = sid_to_sccs[seq.code].split('.')
                    seq_class = seq_sccs[0]
                    seq_fold  = seq_sccs[0..1].join('.')
                    seq_sfam  = seq_sccs[0..2].join('.')
                    seq_fam   = seq_sccs[0..3].join('.')

                    next if stp_class != seq_class

                    psa = ProfileSequenceAligner.new(stp, seq)
                    begin
                      gal = psa.global_alignment_linear_gap_cpp
                    rescue
                      $logger.debug "Problematic prf:"
                      $logger.debug ">#{entry_id}"
                      $logger.debug "Problematic seq:"
                      $logger.debug ">#{seq.code}"
                      $logger.debug "#{seq.data}"
                    end

                    if gal.raw_score < gal.reverse_score
                      $logger.debug "Skip #{sccs} <=> #{seq.code}, reverse raw score is greater than raw score."
                      next
                    end

                    result = [ entry_id, stp.length, stp_class, stp_fold, stp_sfam, stp_fam,
                      sid_to_sunid[seq.code], seq.code, seq.length, seq_class, seq_fold, seq_sfam, seq_fam,
                      gal.raw_score, gal.reverse_score, gal.z_score, psa.algorithm.to_s
                    ]

                    log_fmt = "%8d %5d %2s %5s %9s %12s %8d %-10s %5d %2s %5s %9s %12s %7d %7d %10.7f %10s"
                    $logger.info log_fmt % result
                    hits << result
                  end
                end

                sorted_hits = hits.sort_by { |h| h[-2] }.reverse
                res_file    = cur_dir + "../tmp/rec/gig/rna/single_str-single_seq-linear-#{sccs}-#{env}.csv"
                res_file.open("w") do |file|
                  sorted_hits.each do |sorted_hit|
                    file.puts sorted_hit.join(", ")
                  end
                end
              end
            end
          end
        end
      end
    end
  end
end


namespace :bench do
  namespace :rec do
    namespace :gig do
      desc "Test recognition performance of GIGUE with linear gap penalties for DNA-binding protein families"
      task :lindna do

        seqs    = []
        scop40  = cur_dir + "../data/scop/astral-scopdom-seqres-gd-sel-gs-bib-40-1.75.fa"

        Bio::FlatFile.open(Bio::FastaFormat, scop40) do |ff|
          ff.each do |entry|
            if sid_to_sccs[entry.entry_id.gsub(/^[g|e]/, 'd')].match(/^[a|b|c|d|e|f|g]/)
              seqs << Sequence.new(entry.aaseq, entry.entry_id.gsub(/^[g|e]/, 'd'), "sequence;#{entry.definition}")
            else
              $logger.warn "Skip #{entry.definition}, not a true SCOP class"
            end
          end
        end

        fm = ForkManager.new(8)
        fm.manage do
          dna_tems.each_with_index do |dna_tem, i|
            sunid = dna_tem.to_s.match(/(\d+)/)[1]
            sccs  = sunid_to_sccs[sunid]
            unless sccs.match(/^[a|b|c|d|e|f|g]/)
              $logger.warn "Skip #{sccs}, it's not a true SCOP class"
              next
            end
            $logger.info  "Searching against DNA-bindnig SCOP family, #{sccs} (#{i+1}/#{dna_tems.size})"

            fm.fork do
              [ mat_dir + "ulla-dna64-pid60-sigma0.001-out2.mat",
                mat_dir + "ulla-dna128-pid60-sigma0.001-out2.mat"].each do |esst|
                env   = esst.to_s.match(/(dna\d+)/)[1]
                stp   = StructuralProfile::create_from_joy_tem_and_essts(dna_tem, esst)
                hits  = []
                seqs.each do |seq|
                  stp_sccs  = sccs.split('.')
                  stp_class = stp_sccs[0]
                  stp_fold  = stp_sccs[0..1].join('.')
                  stp_sfam  = stp_sccs[0..2].join('.')
                  stp_fam   = stp_sccs[0..3].join('.')
                  seq_sccs  = sid_to_sccs[seq.code].split('.')
                  seq_class = seq_sccs[0]
                  seq_fold  = seq_sccs[0..1].join('.')
                  seq_sfam  = seq_sccs[0..2].join('.')
                  seq_fam   = seq_sccs[0..3].join('.')

                  next if stp_class != seq_class

                  psa = ProfileSequenceAligner.new(stp, seq)
                  begin
                    gal = psa.global_alignment_linear_gap_cpp
                  rescue
                    $logger.debug "Problematic prf:"
                    $logger.debug ">#{dna_tem}"
                    $logger.debug "Problematic seq:"
                    $logger.debug ">#{seq.code}"
                    $logger.debug "#{seq.data}"
                  end

                  if gal.raw_score < gal.reverse_score
                    $logger.debug "Skip #{sccs} <=> #{seq.code}, reverse raw score is greater than raw score."
                    next
                  end

                  result = [ sunid, stp.length, stp_class, stp_fold, stp_sfam, stp_fam,
                    sid_to_sunid[seq.code], seq.code, seq.length, seq_class, seq_fold, seq_sfam, seq_fam,
                    gal.raw_score, gal.reverse_score, gal.z_score, psa.algorithm.to_s
                  ]

                  log_fmt = "%8d %5d %2s %5s %9s %12s %8d %-10s %5d %2s %5s %9s %12s %7d %7d %10.7f %10s"
                  $logger.info log_fmt % result
                  hits << result
                end

                sorted_hits = hits.sort_by { |h| h[-2] }.reverse
                res_file    = cur_dir + "../tmp/rec/gig/dna/multi_str-single_seq-linear-#{sccs}-#{env}.csv"
                res_file.open("w") do |file|
                  sorted_hits.each do |sorted_hit|
                    file.puts sorted_hit.join(", ")
                  end
                end
              end
            end
          end
        end
      end
    end
  end
end


namespace :bench do
  namespace :rec do
    namespace :gig do
      desc "Test recognition performance of GIGUE with linear gap penalties for RNA-binding protein families"
      task :linrna do

        seqs    = []
        scop40  = cur_dir + "../data/scop/astral-scopdom-seqres-gd-sel-gs-bib-40-1.75.fa"

        Bio::FlatFile.open(Bio::FastaFormat, scop40) do |ff|
          ff.each do |entry|
            if sid_to_sccs[entry.entry_id.gsub(/^[g|e]/, 'd')].match(/^[a|b|c|d|e|f|g]/)
              seqs << Sequence.new(entry.aaseq, entry.entry_id.gsub(/^[g|e]/, 'd'), "sequence;#{entry.definition}")
            else
              $logger.warn "Skip #{entry.definition}, not a true SCOP class"
            end
          end
        end

        fm = ForkManager.new(8)
        fm.manage do
          rna_tems.each_with_index do |rna_tem, i|
            sunid = rna_tem.to_s.match(/(\d+)/)[1]
            sccs  = sunid_to_sccs[sunid]
            unless sccs.match(/^[a|b|c|d|e|f|g]/)
              $logger.warn "Skip #{sccs}, it's not a true SCOP class"
              next
            end
            $logger.info  "Searching against RNA-bindnig SCOP family, #{sccs} (#{i+1}/#{rna_tems.size})"

            fm.fork do
              [ mat_dir + "ulla-rna64-pid60-sigma0.001-out2.mat",
                mat_dir + "ulla-rna128-pid60-sigma0.001-out2.mat"].each do |esst|
                env   = esst.to_s.match(/(rna\d+)/)[1]
                stp   = StructuralProfile::create_from_joy_tem_and_essts(rna_tem, esst)
                hits  = []
                seqs.each do |seq|
                  stp_sccs  = sccs.split('.')
                  stp_class = stp_sccs[0]
                  stp_fold  = stp_sccs[0..1].join('.')
                  stp_sfam  = stp_sccs[0..2].join('.')
                  stp_fam   = stp_sccs[0..3].join('.')
                  seq_sccs  = sid_to_sccs[seq.code].split('.')
                  seq_class = seq_sccs[0]
                  seq_fold  = seq_sccs[0..1].join('.')
                  seq_sfam  = seq_sccs[0..2].join('.')
                  seq_fam   = seq_sccs[0..3].join('.')

                  next if stp_class != seq_class

                  psa = ProfileSequenceAligner.new(stp, seq)
                  begin
                    gal = psa.global_alignment_linear_gap_cpp
                  rescue
                    $logger.debug "Problematic prf:"
                    $logger.debug ">#{rna_tem}"
                    $logger.debug "Problematic seq:"
                    $logger.debug ">#{seq.code}"
                    $logger.debug "#{seq.data}"
                  end

                  if gal.raw_score < gal.reverse_score
                    $logger.debug "Skip #{sccs} <=> #{seq.code}, reverse raw score is greater than raw score."
                    next
                  end

                  result = [ sunid, stp.length, stp_class, stp_fold, stp_sfam, stp_fam,
                    sid_to_sunid[seq.code], seq.code, seq.length, seq_class, seq_fold, seq_sfam, seq_fam,
                    gal.raw_score, gal.reverse_score, gal.z_score, psa.algorithm.to_s
                  ]

                  log_fmt = "%8d %5d %2s %5s %9s %12s %8d %-10s %5d %2s %5s %9s %12s %7d %7d %10.7f %10s"
                  $logger.info log_fmt % result
                  hits << result
                end
                sorted_hits = hits.sort_by { |h| h[-2] }.reverse
                res_file    = cur_dir + "../tmp/rec/gig/rna/multi_str-single_seq-linear-#{sccs}-#{env}.csv"
                res_file.open("w") do |file|
                  sorted_hits.each do |sorted_hit|
                    file.puts sorted_hit.join(", ")
                  end
                end
              end
            end
          end
        end
      end
    end
  end
end


namespace :bench do
  namespace :rec do
    namespace :gig do
      desc "Test recognition performance of GIGUE with affine gap penalties for DNA-binding protein families"
      task :affdna do

        seqs    = []
        scop40  = cur_dir + "../data/scop/astral-scopdom-seqres-gd-sel-gs-bib-40-1.75.fa"

        Bio::FlatFile.open(Bio::FastaFormat, scop40) do |ff|
          ff.each do |entry|
            if sid_to_sccs[entry.entry_id.gsub(/^[g|e]/,'d')].match(/^[a|b|c|d|e|f|g]/)
              seqs << Sequence.new(entry.aaseq, entry.entry_id.gsub(/^[g|e]/,'d'), "sequence;#{entry.definition}")
            else
              $logger.warn "Skip #{entry.definition}, not a true SCOP class"
            end
          end
        end

        fm = ForkManager.new(8)
        fm.manage do
          dna_tems.each_with_index do |dna_tem, i|
            sunid = dna_tem.to_s.match(/(\d+)/)[1]
            sccs  = sunid_to_sccs[sunid]
            unless sccs.match(/^[a|b|c|d|e|f|g]/)
              $logger.warn "Skip #{sccs}, it's not a true SCOP class"
              next
            end
            $logger.info  "Searching against DNA-bindnig SCOP family, #{sccs} (#{i+1}/#{dna_tems.size})"

            fm.fork do
              [ mat_dir + "ulla-dna64-pid60-sigma0.001-out2.mat",
                mat_dir + "ulla-dna128-pid60-sigma0.001-out2.mat"].each do |esst|
                env   = esst.to_s.match(/(dna\d+)/)[1]
                stp   = StructuralProfile::create_from_joy_tem_and_essts(dna_tem, esst)
                hits  = []
                seqs.each do |seq|
                  stp_sccs  = sccs.split('.')
                  stp_class = stp_sccs[0]
                  stp_fold  = stp_sccs[0..1].join('.')
                  stp_sfam  = stp_sccs[0..2].join('.')
                  stp_fam   = stp_sccs[0..3].join('.')
                  seq_sccs  = sid_to_sccs[seq.code].split('.')
                  seq_class = seq_sccs[0]
                  seq_fold  = seq_sccs[0..1].join('.')
                  seq_sfam  = seq_sccs[0..2].join('.')
                  seq_fam   = seq_sccs[0..3].join('.')

                  next if stp_class != seq_class

                  psa = ProfileSequenceAligner.new(stp, seq)
                  begin
                    gal = psa.global_alignment_affine_gap_cpp
                  rescue
                    $logger.debug "Problematic prf:"
                    $logger.debug ">#{dna_tem}"
                    $logger.debug "Problematic seq:"
                    $logger.debug ">#{seq.code}"
                    $logger.debug "#{seq.data}"
                  end

                  if gal.raw_score < gal.reverse_score
                    $logger.debug "Skip #{sccs} <=> #{seq.code}, reverse raw score is greater than raw score."
                    next
                  end

                  result = [ sunid, stp.length, stp_class, stp_fold, stp_sfam, stp_fam,
                    sid_to_sunid[seq.code], seq.code, seq.length, seq_class, seq_fold, seq_sfam, seq_fam,
                    gal.raw_score, gal.reverse_score, gal.z_score, psa.algorithm.to_s
                  ]

                  log_fmt = "%8d %5d %2s %5s %9s %12s %8d %-10s %5d %2s %5s %9s %12s %7d %7d %10.7f %10s"
                  $logger.info log_fmt % result
                  hits << result
                end

                sorted_hits = hits.sort_by { |h| h[-2] }.reverse
                res_file    = cur_dir + "../tmp/rec/gig/dna/multi_str-single_seq-affine-#{sccs}-#{env}.csv"
                res_file.open("w") do |file|
                  sorted_hits.each do |sorted_hit|
                    file.puts sorted_hit.join(", ")
                  end
                end
              end
            end
          end
        end
      end
    end
  end
end


namespace :bench do
  namespace :rec do
    namespace :gig do
      desc "Test recognition performance of GIGUE with affine gap penalties for RNA-binding protein families"
      task :affrna do

        seqs    = []
        scop40  = cur_dir + "../data/scop/astral-scopdom-seqres-gd-sel-gs-bib-40-1.75.fa"

        Bio::FlatFile.open(Bio::FastaFormat, scop40) do |ff|
          ff.each do |entry|
            if sid_to_sccs[entry.entry_id.gsub(/^[g|e]/,'d')].match(/^[a|b|c|d|e|f|g]/)
              seqs << Sequence.new(entry.aaseq, entry.entry_id.gsub(/^[g|e]/,'d'), "sequence;#{entry.definition}")
            else
              $logger.warn "Skip #{entry.definition}, not a true SCOP class"
            end
          end
        end

        fm = ForkManager.new(8)
        fm.manage do
          rna_tems.each_with_index do |rna_tem, i|
            sunid = rna_tem.to_s.match(/(\d+)/)[1]
            sccs  = sunid_to_sccs[sunid]
            unless sccs.match(/^[a|b|c|d|e|f|g]/)
              $logger.warn "Skip #{sccs}, it's not a true SCOP class"
              next
            end
            $logger.info  "Searching against RNA-bindnig SCOP family, #{sccs} (#{i+1}/#{rna_tems.size})"

            fm.fork do
              [ mat_dir + "ulla-rna64-pid60-sigma0.001-out2.mat",
                mat_dir + "ulla-rna128-pid60-sigma0.001-out2.mat"].each do |esst|
                env   = esst.to_s.match(/(rna\d+)/)[1]
                stp   = StructuralProfile::create_from_joy_tem_and_essts(rna_tem, esst)
                hits  = []
                seqs.each do |seq|
                  stp_sccs  = sccs.split('.')
                  stp_class = stp_sccs[0]
                  stp_fold  = stp_sccs[0..1].join('.')
                  stp_sfam  = stp_sccs[0..2].join('.')
                  stp_fam   = stp_sccs[0..3].join('.')
                  seq_sccs  = sid_to_sccs[seq.code].split('.')
                  seq_class = seq_sccs[0]
                  seq_fold  = seq_sccs[0..1].join('.')
                  seq_sfam  = seq_sccs[0..2].join('.')
                  seq_fam   = seq_sccs[0..3].join('.')

                  next if stp_class != seq_class

                  psa = ProfileSequenceAligner.new(stp, seq)
                  begin
                    gal = psa.global_alignment_affine_gap_cpp
                  rescue
                    $logger.debug "Problematic prf:"
                    $logger.debug ">#{rna_tem}"
                    $logger.debug "Problematic seq:"
                    $logger.debug ">#{seq.code}"
                    $logger.debug "#{seq.data}"
                  end

                  if gal.raw_score < gal.reverse_score
                    $logger.debug "Skip #{sccs} <=> #{seq.code}, reverse raw score is greater than raw score."
                    next
                  end

                  result = [
                    sunid, stp.length, stp_class, stp_fold, stp_sfam, stp_fam,
                    sid_to_sunid[seq.code], seq.code, seq.length, seq_class, seq_fold, seq_sfam, seq_fam,
                    gal.raw_score, gal.reverse_score, gal.z_score, psa.algorithm.to_s
                  ]

                  log_fmt = "%8d %5d %2s %5s %9s %12s %8d %-10s %5d %2s %5s %9s %12s %7d %7d %10.7f %10s"
                  $logger.info log_fmt % result
                  hits << result
                end
                sorted_hits = hits.sort_by { |h| h[-2] }.reverse
                res_file    = cur_dir + "../tmp/rec/gig/rna/multi_str-single_seq-affine-#{sccs}-#{env}.csv"
                res_file.open("w") do |file|
                  sorted_hits.each do |sorted_hit|
                    file.puts sorted_hit.join(", ")
                  end
                end
              end
            end
          end
        end
      end
    end
  end
end
