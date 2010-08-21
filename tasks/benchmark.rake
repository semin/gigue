include Gigue
require 'forkoff'

cur_dir = Pathname.new(__FILE__).dirname

namespace :bench do

  desc "Remote homology recognition performance"
  task :recog do

    $logger.level   = Logger::WARN
    log_fmt_title   = "%-27s %-30s %7s %7s %7s %7s %7s"
    log_fmt_result  = "%-27s %-30s %7d %7d %7d %7d %7.2f"
    log_titles      = %w[PROFILE SEQUENCE PLEN SLEN RAWS RVN ZSCORE]

    # load SCOP description file
    sunid_to_sccs = {}
    sid_to_sccs   = {}
    scop_des = cur_dir + "../data/dir.des.scop.txt_1.75"

    IO.foreach(scop_des) do |line|
      line.chomp!
      if line =~ /^#/
        next
      elsif (cols = line.split("\t")).size == 5
        sunid_to_sccs[cols[0]] = cols[2]
        sid_to_sccs[cols[3]] = cols[2] if cols[1] == "px"
      else
        warn "? #{line}"
      end
    end

    # load structural profiles
    dna_esst64  = cur_dir + "../data/ulla-ent-dna-esst64-pid60-sigma5-logo.mat"
    dna_esst128 = cur_dir + "../data/ulla-ent-dna-esst128-pid60-sigma5-logo.mat"
    dna_tems    = Pathname::glob(cur_dir + "../data/bipa/scop/rep/dna/*/dnamodsalign.tem")

    ## load sequences
    seqs = []

    #scop40 = cur_dir + "../data/astral-scopdom-seqres-gd-sel-gs-bib-40-1.75.fa"
    #Bio::FlatFile.open(Bio::FastaFormat, scop40) do |fasta_file|
      #fasta_file.each do |entry|
        #seqs << Sequence.new(entry.aaseq, entry.entry_id, "sequence;#{entry.definition}")
      #end
    #end

    dna_tems.each do |dna_tem|
      Bio::FlatFile.auto(dna_tem) do |pir_file|
        pir_file.each do |entry|
          if (entry.definition =~ /^sequence/)
            args = [entry.aaseq.gsub('-',''), entry.entry_id, entry.definition]
            seqs << Sequence.new(*args)
          end
        end
      end
    end

    dna_tems.each_with_index do |dna_tem, i|
      hits  = []
      sunid = dna_tem.to_s.match(/(\d+)/)[1]
      stp   = StructuralProfile.new(dna_tem, dna_esst64)
      $logger.warn  "Searching against DNA-bindnig SCOP family, #{sunid_to_sccs[sunid]} (#{i+1}/#{dna_tems.size})"
      $logger.warn log_fmt_title % log_titles

      hits = seqs.forkoff!(:processes => 2) do |seq|
        psa = ProfileSequenceAligner.new(stp, seq)
        gal = psa.global_alignment_linear_gap_cpp
        result = [
          "#{sunid} (#{sunid_to_sccs[sunid]})",
          "#{seq.code} (#{sunid_to_sccs[seq.code]})",
          stp.length,
          seq.length,
          gal.raw_score,
          gal.reverse_score,
          gal.z_score
        ]
        $logger.warn log_fmt_result % result
        result
      end

      puts log_fmt_title % log_titles
      hits = hits.sort_by { |h| h[-1] }
      hits.each do |hit|
        puts log_fmt_result % hit
      end
    end

  end
end
