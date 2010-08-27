include Gigue
require 'fork_manager'

cur_dir = Pathname.new(__FILE__).dirname

# logger settings
$logger.level   = Logger::INFO

# load SCOP description file
sunid_to_sccs = {}
sid_to_sccs   = {}
sid_to_sunid  = {}
scop_des_file = cur_dir + "../data/dir.des.scop.txt_1.75"

IO.foreach(scop_des_file) do |line|
  line.chomp!
  if line =~ /^#/
    next
  elsif (cols = line.split("\t")).size == 5
    sunid_to_sccs[cols[0]] = cols[2]
    sid_to_sccs[cols[3]] = cols[2] if cols[1] == "px"
    sid_to_sunid[cols[3]] = cols[0] if cols[1] == "px"
  else
    warn "? #{line}"
  end
end

dna_tems = Pathname::glob(cur_dir + "../data/bipa/scop/rep/dna/*/dnamodsalign.tem")
rna_tems = Pathname::glob(cur_dir + "../data/bipa/scop/rep/rna/*/rnamodsalign.tem")

dna_esst64_sigma5       = cur_dir + "../data/ulla-ent-dna-esst64-pid60-sigma5-logo.mat"
dna_esst128_sigma5      = cur_dir + "../data/ulla-ent-dna-esst128-pid60-sigma5-logo.mat"
dna_esst64_sigma1       = cur_dir + "../data/ulla-ent-dna-esst64-pid60-sigma1-logo.mat"
dna_esst128_sigma1      = cur_dir + "../data/ulla-ent-dna-esst128-pid60-sigma1-logo.mat"
dna_esst64_sigma0_1     = cur_dir + "../data/ulla-ent-dna-esst64-pid60-sigma0.1-logo.mat",
dna_esst128_sigma0_1    = cur_dir + "../data/ulla-ent-dna-esst128-pid60-sigma0.1-logo.mat",
dna_esst64_sigma0_01    = cur_dir + "../data/ulla-ent-dna-esst64-pid60-sigma0.01-logo.mat",
dna_esst128_sigma0_01   = cur_dir + "../data/ulla-ent-dna-esst128-pid60-sigma0.01-logo.mat",
dna_esst64_sigma0_001   = cur_dir + "../data/ulla-ent-dna-esst64-pid60-sigma0.001-logo.mat",
dna_esst128_sigma0_001  = cur_dir + "../data/ulla-ent-dna-esst128-pid60-sigma0.001-logo.mat",
dna_esst64_sigma0_002   = cur_dir + "../data/ulla-ent-dna-esst64-pid60-sigma0.002-logo.mat",
dna_esst128_sigma0_002  = cur_dir + "../data/ulla-ent-dna-esst128-pid60-sigma0.002-logo.mat",

rna_esst64_sigma5       = cur_dir + "../data/ulla-ent-rna-esst64-pid60-sigma5-logo.mat"
rna_esst128_sigma5      = cur_dir + "../data/ulla-ent-rna-esst128-pid60-sigma5-logo.mat"
rna_esst64_sigma1       = cur_dir + "../data/ulla-ent-rna-esst64-pid60-sigma1-logo.mat"
rna_esst128_sigma1      = cur_dir + "../data/ulla-ent-rna-esst128-pid60-sigma1-logo.mat"
rna_esst64_sigma0_1     = cur_dir + "../data/ulla-ent-rna-esst64-pid60-sigma0.1-logo.mat",
rna_esst128_sigma0_1    = cur_dir + "../data/ulla-ent-rna-esst128-pid60-sigma0.1-logo.mat",
rna_esst64_sigma0_01    = cur_dir + "../data/ulla-ent-rna-esst64-pid60-sigma0.01-logo.mat",
rna_esst128_sigma0_01   = cur_dir + "../data/ulla-ent-rna-esst128-pid60-sigma0.01-logo.mat",
rna_esst64_sigma0_001   = cur_dir + "../data/ulla-ent-rna-esst64-pid60-sigma0.001-logo.mat",
rna_esst128_sigma0_001  = cur_dir + "../data/ulla-ent-rna-esst128-pid60-sigma0.001-logo.mat",
rna_esst64_sigma0_002   = cur_dir + "../data/ulla-ent-rna-esst64-pid60-sigma0.002-logo.mat",
rna_esst128_sigma0_002  = cur_dir + "../data/ulla-ent-rna-esst128-pid60-sigma0.002-logo.mat",

def calculate_pid(aas1, aas2)
  cols  = aas1.zip(aas2)
  gap   = '-'
  align = 0.0 # no. of aligned columns
  ident = 0.0 # no. of identical columns
  intgp = 0.0 # no. of internal gaps
  aafnd = false

  cols.each do |col|
    if (col[0] != gap) && (col[1] != gap)
      align += 1
      aafnd = true
      ident += 1 if col[0] == col[1]
    elsif (((col[0] == gap) && (col[1] != gap)) ||
           ((col[0] != gap) && (col[1] == gap)))
      intgp += 1 if aafnd
    elsif (col[0] == gap) && (col[1] == gap)
      next
    else
      warn "Unknown combination!"
      exit 1
    end
  end
  Float(ident) / (align + intgp)
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

  desc "Accuracy of alignments"
  task :align2 do
    rna_tems.each do |tem|
      fam_sunid     = tem.to_s.match(/(\d+)/)[1]
      fam_sccs      = sunid_to_sccs[fam_sunid]
      bio_tem       = Bio::FlatFile.auto(tem)
      entries       = bio_tem.entries
      entries_hash  = {}
      csv_file      = cur_dir + "../tmp/align/rna/#{fam_sccs}.csv"

      csv_file.open("w") do |file|
        entries.each do |entry|
          if !entries_hash.has_key?(entry.entry_id)
            entries_hash[entry.entry_id] = {}
          end
          entries_hash[entry.entry_id][entry.definition] = entry
        end

        log_title = "%-20s %-20s %-40s %8s %8s %10s" % %w[PRF SEQ ESST REF_PID GIG_PID S0]
        puts log_title
        puts "-" * 111

        entries_hash.keys.each do |ref_entry_id|
          ref = entries_hash[ref_entry_id]
          joy = ""
          ref.values.each do |e|
            joy += %Q{>P1;#{e.entry_id}\n#{e.definition}\n#{e.data.gsub('-','').gsub("\n",'')}*\n}
          end

          others = entries_hash.reject { |k, v| k == ref_entry_id }
          others.each do |qry_entry_id, qry_entry_hash|
            [
              cur_dir + "../data/ulla-ent-rna-esst64-pid60-sigma5-logo.mat",
              cur_dir + "../data/ulla-ent-rna-esst128-pid60-sigma5-logo.mat",
              cur_dir + "../data/ulla-ent-rna-esst64-pid60-sigma1-logo.mat",
              cur_dir + "../data/ulla-ent-rna-esst128-pid60-sigma1-logo.mat",
              cur_dir + "../data/ulla-ent-rna-esst64-pid60-sigma0.1-logo.mat",
              cur_dir + "../data/ulla-ent-rna-esst128-pid60-sigma0.1-logo.mat",
              cur_dir + "../data/ulla-ent-rna-esst64-pid60-sigma0.01-logo.mat",
              cur_dir + "../data/ulla-ent-rna-esst128-pid60-sigma0.01-logo.mat",
              #cur_dir + "../data/ulla-ent-rna-esst64-pid60-sigma0.001-logo.mat",
              #cur_dir + "../data/ulla-ent-rna-esst128-pid60-sigma0.001-logo.mat",
              cur_dir + "../data/ulla-ent-rna-esst64-pid60-sigma0.002-logo.mat",
              cur_dir + "../data/ulla-ent-rna-esst128-pid60-sigma0.002-logo.mat",
              #cur_dir + "../data/ulla-logo-smooth-toccata-maskA.mat",
              #cur_dir + "../data/ulla-mlr-logo-toccata-maskA.mat",
            ].each do |esst|
              str_prf     = StructuralProfile.new(joy, esst)
              ref_seq_ent = ref['sequence']
              ref_qry_ent = qry_entry_hash['sequence']
              ref_str_seq = Sequence.new(ref_seq_ent.aaseq, ref_seq_ent.entry_id, ref_seq_ent.definition)
              ref_qry_seq = Sequence.new(ref_qry_ent.aaseq, ref_qry_ent.entry_id, ref_qry_ent.definition)
              ref_str_seq_cols  = match_positions(ref_seq_ent.aaseq.split(''), ref_qry_seq.amino_acids)
              ref_str_seq_pid   = calculate_pid(ref_seq_ent.aaseq.split(''), ref_qry_seq.amino_acids)
              #ref_str_seq.to_flatfile
              #ref_qry_seq.to_flatfile

              raw_qry_seq = Sequence.new(ref_qry_ent.aaseq.gsub('-',''), ref_qry_ent.entry_id, ref_qry_ent.definition)
              psa         = ProfileSequenceAligner.new(str_prf, raw_qry_seq)
              gal         = psa.global_alignment_affine_gap_cpp
              gig_str_seq_cols  = match_positions(gal.aligned_structural_profile_positions.map(&:probe), gal.aligned_amino_acids)
              gig_str_seq_pid   = calculate_pid(gal.aligned_structural_profile_positions.map(&:probe), gal.aligned_amino_acids)
              gig_prf_seq = Sequence.new(gal.aligned_structural_profile_positions.map(&:probe).join(''), ref_entry_id, ref_seq_ent.definition)
              gig_qry_seq = Sequence.new(gal.aligned_amino_acids.join(''), qry_entry_hash, ref_qry_ent.definition)
              #gig_prf_seq.to_flatfile
              #gig_qry_seq.to_flatfile

              mat_cols  = (ref_str_seq_cols & gig_str_seq_cols).size
              ali_score = mat_cols / Float(ref_str_seq_cols.size)

              elems = [
                "#{ref_entry_id} (#{sunid_to_sccs[ref_entry_id]})",
                "#{qry_entry_id} (#{sunid_to_sccs[qry_entry_id]})",
                  File.basename(esst, ".mat").gsub("ulla-ent-rna-","").gsub("-logo",""),
                    ref_str_seq_pid,
                    gig_str_seq_pid,
                    ali_score
              ]

              file.puts elems.join(", ")

              log_res   = "%-20s %-20s %-40s %8.6f %8.6f %10.8f" % elems
              puts log_res
            end
            puts "-" * 111
          end
        end
      end # csv.open
    end
  end

  desc "Accuracy of alignments"
  task :align do
    dna_tems[2..-1].each do |tem|
      fam_sunid     = tem.to_s.match(/(\d+)/)[1]
      fam_sccs      = sunid_to_sccs[fam_sunid]
      bio_tem       = Bio::FlatFile.auto(tem)
      entries       = bio_tem.entries
      entries_hash  = {}
      csv_file      = cur_dir + "../tmp/align/dna/#{fam_sccs}.csv"

      csv_file.open("w") do |file|
        entries.each do |entry|
          if !entries_hash.has_key?(entry.entry_id)
            entries_hash[entry.entry_id] = {}
          end
          entries_hash[entry.entry_id][entry.definition] = entry
        end

        log_title = "%-20s %-20s %-40s %8s %8s %10s" % %w[PRF SEQ ESST REF_PID GIG_PID S0]
        puts log_title
        puts "-" * 115

        entries_hash.keys.each do |ref_entry_id|
          ref = entries_hash[ref_entry_id]
          joy = ""
          ref.values.each do |e|
            joy += %Q{>P1;#{e.entry_id}\n#{e.definition}\n#{e.data.gsub('-','').gsub("\n",'')}*\n}
          end

          others = entries_hash.reject { |k, v| k == ref_entry_id }
          others.each do |qry_entry_id, qry_entry_hash|
            [
              #cur_dir + "../data/subst-ent-dna-esst128-pid60-sigma5-logo.mat",
              #cur_dir + "../data/subst-ent-dna-esst64-pid60-sigma5-logo.mat",
              cur_dir + "../data/ulla-ent-dna-esst64-pid60-sigma1-logo.mat",
              cur_dir + "../data/ulla-ent-dna-esst128-pid60-sigma1-logo.mat",
              cur_dir + "../data/ulla-ent-dna-esst64-pid60-sigma0.1-logo.mat",
              cur_dir + "../data/ulla-ent-dna-esst128-pid60-sigma0.1-logo.mat",
              cur_dir + "../data/ulla-ent-dna-esst64-pid60-sigma0.01-logo.mat",
              cur_dir + "../data/ulla-ent-dna-esst128-pid60-sigma0.01-logo.mat",
              #cur_dir + "../data/ulla-ent-dna-esst64-pid60-sigma0.001-logo.mat",
              #cur_dir + "../data/ulla-ent-dna-esst128-pid60-sigma0.001-logo.mat",
              cur_dir + "../data/ulla-ent-dna-esst64-pid60-sigma0.002-logo.mat",
              cur_dir + "../data/ulla-ent-dna-esst128-pid60-sigma0.002-logo.mat",
              #cur_dir + "../data/ulla-logo-smooth-toccata-maskA.mat",
              #cur_dir + "../data/ulla-mlr-logo-toccata-maskA.mat",
            ].each do |esst|
              str_prf     = StructuralProfile.new(joy, esst)
              ref_seq_ent = ref['sequence']
              ref_qry_ent = qry_entry_hash['sequence']
              ref_str_seq = Sequence.new(ref_seq_ent.aaseq, ref_seq_ent.entry_id, ref_seq_ent.definition)
              ref_qry_seq = Sequence.new(ref_qry_ent.aaseq, ref_qry_ent.entry_id, ref_qry_ent.definition)
              ref_str_seq_cols  = match_positions(ref_seq_ent.aaseq.split(''), ref_qry_seq.amino_acids)
              ref_str_seq_pid   = calculate_pid(ref_seq_ent.aaseq.split(''), ref_qry_seq.amino_acids)
              #ref_str_seq.to_flatfile
              #ref_qry_seq.to_flatfile

              raw_qry_seq = Sequence.new(ref_qry_ent.aaseq.gsub('-',''), ref_qry_ent.entry_id, ref_qry_ent.definition)
              psa         = ProfileSequenceAligner.new(str_prf, raw_qry_seq)
              gal         = psa.global_alignment_affine_gap_cpp
              gig_str_seq_cols  = match_positions(gal.aligned_structural_profile_positions.map(&:probe), gal.aligned_amino_acids)
              gig_str_seq_pid   = calculate_pid(gal.aligned_structural_profile_positions.map(&:probe), gal.aligned_amino_acids)
              gig_prf_seq = Sequence.new(gal.aligned_structural_profile_positions.map(&:probe).join(''), ref_entry_id, ref_seq_ent.definition)
              gig_qry_seq = Sequence.new(gal.aligned_amino_acids.join(''), qry_entry_hash, ref_qry_ent.definition)
              #gig_prf_seq.to_flatfile
              #gig_qry_seq.to_flatfile

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

              log_res   = "%-20s %-20s %-40s %8.6f %8.6f %10.8f" % elems
              puts log_res
            end
            puts "-" * 115
          end
        end
      end # csv.open
    end
  end

  desc "Remote homology recognition performance"
  task :recog_dna do

    seqs = []

    scop20 = cur_dir + "../data/astral-scopdom-seqres-gd-sel-gs-bib-20-1.75.fa"
    Bio::FlatFile.open(Bio::FastaFormat, scop20) do |ff|
      ff.each do |entry|
        if sid_to_sccs[entry.entry_id.gsub(/^g/, 'd').gsub(/^e/, 'd')].match(/^[a|b|c|d|e|f|g]/)
          seqs << Sequence.new(entry.aaseq.gsub('X',''),
                               entry.entry_id.gsub(/^g/, 'd').gsub(/^e/, 'd'),
                               "sequence;#{entry.definition}")
        else
          $logger.warn "Skip #{entry.definition}, not a true SCOP class"
        end
      end
    end

    #dna_tems.each do |dna_tem|
      #Bio::FlatFile.auto(dna_tem) do |pir_file|
        #pir_file.each do |entry|
          #if (entry.definition =~ /^sequence/)
            #args = [entry.aaseq.gsub('-',''), entry.entry_id, "sequence;#{entry.definition}"]
            #seqs << Sequence.new(*args)
          #end
        #end
      #end
    #end

    log_fmt = "%8d %5d %2s %6s %10s %14s %8d %8s %5d %2s %6s %10s %14s %7d %7d %6.4 %10s"

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
          [ dna_esst64_sigma0_002, dna_esst128_sigma0_002 ].each do |esst|
            env   = esst.to_s.match(/(esst\d+)/)[1]
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

              if stp_class != seq_class
                #$logger.warn "Skip comparison between #{stp_sccs} and #{seq_sccs}"
                next
              end

              psa = ProfileSequenceAligner.new(stp, seq)
              begin
                gal = psa.global_alignment_linear_gap_cpp
              rescue
                puts "Problematic prf:"
                puts ">#{dna_tem}"
                puts "Problematic seq:"
                puts ">#{seq.code}"
                puts "#{seq.data}"
              end
              result    = [
                sunid,
                stp.length,
                stp_class,
                stp_fold,
                stp_sfam,
                stp_fam,
                sid_to_sunid[seq.code],
                seq.code,
                seq.length,
                seq_class,
                seq_fold,
                seq_sfam,
                seq_fam,
                gal.raw_score,
                gal.reverse_score,
                gal.z_score,
                psa.algorithm.to_s
              ]
              $logger.info log_fmt % result
              hits << result
            end

            sorted_hits = hits.sort_by { |h| h[-2] }.reverse
            res_file    = cur_dir + "../tmp/recog/dna/#{sccs}-multi_str-single_seq-#{env}.csv"
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


  desc "Remote homology recognition performance"
  task :recog_rna do

    seqs    = []
    scop20  = cur_dir + "../data/astral-scopdom-seqres-gd-sel-gs-bib-20-1.75.fa"

    Bio::FlatFile.open(Bio::FastaFormat, scop20) do |ff|
      ff.each do |entry|
        if sid_to_sccs[entry.entry_id.gsub(/^g/, 'd').gsub(/^e/, 'd')].match(/^[a|b|c|d|e|f|g]/)
          seqs << Sequence.new(entry.aaseq.gsub('X',''),
                               entry.entry_id.gsub(/^g/, 'd').gsub(/^e/, 'd'),
                               "sequence;#{entry.definition}")
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
          [ rna_esst64_sigma0_002, rna_esst128_sigma0_002 ].each do |esst|
            env   = esst.to_s.match(/(esst\d+)/)[1]
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

              if stp_class != seq_class
                #$logger.warn "Skip comparison between #{stp_sccs} and #{seq_sccs}"
                next
              end

              psa = ProfileSequenceAligner.new(stp, seq)
              begin
                gal = psa.global_alignment_linear_gap_cpp
              rescue
                puts "Problematic prf:"
                puts ">#{rna_tem}"
                puts "Problematic seq:"
                puts ">#{seq.code}"
                puts "#{seq.data}"
              end
              result    = [
                sunid,
                stp.length,
                stp_class,
                stp_fold,
                stp_sfam,
                stp_fam,
                sid_to_sunid[seq.code],
                seq.code,
                seq.length,
                seq_class,
                seq_fold,
                seq_sfam,
                seq_fam,
                gal.raw_score,
                gal.reverse_score,
                gal.z_score,
                psa.algorithm.to_s
              ]
              #    46562    73  a    a.2      a.2.2        a.2.2.1    15737  d1coja1    89  a    a.2     a.2.11       a.2.11.1   -2588   -2854 1.9448     global

              log_fmt = "%8d %5d %2s %5s %9s %12s %8d %-10s %5d %2s %5s %9s %12s %7d %7d %10.7f %10s"
              $logger.info log_fmt % result
              hits << result
            end

            sorted_hits = hits.sort_by { |h| h[-2] }.reverse
            res_file    = cur_dir + "../tmp/recog/rna/#{sccs}-multi_str-single_seq-#{env}.csv"
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

  desc "Remote homology recognition performance"
  task :recog_dna_affine do

    seqs    = []
    scop40  = cur_dir + "../data/astral-scopdom-seqres-gd-sel-gs-bib-40-1.75.fa"

    Bio::FlatFile.open(Bio::FastaFormat, scop40) do |ff|
      ff.each do |entry|
        if sid_to_sccs[entry.entry_id.gsub(/^g/, 'd').gsub(/^e/, 'd')].match(/^[a|b|c|d|e|f|g]/)
          seqs << Sequence.new(entry.aaseq.gsub('X','').gsub('Z','Q'),
                               entry.entry_id.gsub(/^g/, 'd').gsub(/^e/, 'd'),
                               "sequence;#{entry.definition}")
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
          [ dna_esst64_sigma0_002, dna_esst128_sigma0_002 ].each do |esst|
            env   = esst.to_s.match(/(esst\d+)/)[1]
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

              if stp_class != seq_class
                $logger.warn "Skip comparison between #{stp_sccs} and #{seq_sccs}"
                next
              end

              psa = ProfileSequenceAligner.new(stp, seq)
              begin
                gal = psa.global_alignment_affine_gap_cpp
              rescue
                puts "Problematic prf:"
                puts ">#{dna_tem}"
                puts "Problematic seq:"
                puts ">#{seq.code}"
                puts "#{seq.data}"
              end
              result    = [
                sunid,
                stp.length,
                stp_class,
                stp_fold,
                stp_sfam,
                stp_fam,
                sid_to_sunid[seq.code],
                seq.code,
                seq.length,
                seq_class,
                seq_fold,
                seq_sfam,
                seq_fam,
                gal.raw_score,
                gal.reverse_score,
                gal.z_score,
                psa.algorithm.to_s
              ]
              #    46562    73  a    a.2      a.2.2        a.2.2.1    15737  d1coja1    89  a    a.2     a.2.11       a.2.11.1   -2588   -2854 1.9448     global

              log_fmt = "%8d %5d %2s %5s %9s %12s %8d %-10s %5d %2s %5s %9s %12s %7d %7d %10.7f %10s"
              $logger.info log_fmt % result
              hits << result
            end

            sorted_hits = hits.sort_by { |h| h[-2] }.reverse
            res_file    = cur_dir + "../tmp/recog/dna/#{sccs}-multi_str-single_seq-#{env}-affine.csv"
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
