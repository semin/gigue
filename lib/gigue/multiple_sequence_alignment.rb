module Gigue
  class MultipleSequenceAlignment

    def self.create_from_psiblast_output_style6(file, weighting=:va)
      # find maximum interation no.
      max_iter = 0

      IO.foreach(file) do |line|
        if line =~ /^Results\s+from\s+round\s+(\d+)/
          iter      = Integer($1)
          max_iter  = iter if iter > max_iter
        end
      end

      # parse the final multiple sequence alignment
      parse_tag     = false
      code_to_data  = {}
      code_to_start = {}
      code_to_end   = {}

      IO.foreach(file) do |line|
        if    line =~ /^Results\s+from\s+round\s+#{max_iter}/
          parse_tag = true
        elsif line =~ /^(\S+)\s+(\d+|)\s+(\S+)\s*(\d+|)\s*$/ && parse_tag
          seq_code  = $1
          seq_start = $2.to_i
          seq_data  = $3
          seq_end   = $4.to_i
          if code_to_data.has_key? seq_code
            code_to_data[seq_code] += seq_data
          else
            code_to_data[seq_code] = seq_data
          end
          if code_to_start.has_key?(seq_code) && code_to_start[seq_code] > 0
            if seq_start < code_to_start[seq_code]
              code_to_start[seq_code] = seq_start
            end
          else
            code_to_start[seq_code] = seq_start
          end
          if code_to_end.has_key?(seq_code) && code_to_end[seq_code] > 0
            if seq_end > code_to_end[seq_code]
              code_to_end[seq_code] = seq_end
            end
          else
            code_to_end[seq_code] = seq_end
          end
        elsif line =~ /^\s+Database:/ && parse_tag
          break
        end
      end

      # create MultipleSequenceAlignment instance using Sequence instances
      seqs = []
      code_to_data.each do |code, data|
        seqs << Sequence.new(data, "#{code}|#{code_to_start[code]}-#{code_to_end[code]}", "sequence")
      end
      self.new(seqs, weighting)
    end


    attr_reader :sequences, :length, :depth, :columns

    def initialize(seqs, weighting=:va)
      @sequences  = seqs
      @depth      = seqs.size
      @length     = seqs[0].length
      @columns    = (0...@length).map do |mi|
        MultipleSequenceAlignmentColumn.new(@sequences.inject('') { |p, s| p + s.data[mi] })
      end
      @weighting  = weighting
      if    weighting == :va
        Sequence::calculate_va_weights_cpp(@sequences)
      elsif weighting == :blosum
        Sequence::calculate_blosum_weights(@sequences)
      else
        $logger.error "Unknown weighting scheme!"
        exit 1
      end
    end

    def shuffle
      seqs = @sequences.dup
      for i in 0...@length
        r = Kernel.rand(@length-i)+i
        seqs.each do |seq|
          seq[r], seq[i] = seq[i], seq[r]
        end
      end
      self.class.new(seqs, weighting)
    end

    def to_sequence_profile
      prf_pss = []

      @columns.each do |col|
        probe       = col.probe
        aa_raw_frqs = Hash.new(0)
        aa_frqs     = Hash.new(0.0)
        aa_prbs     = Hash.new(0.0)
        aa_rel_prbs = Hash.new(0.0)

        probe.split('').each_with_index do |aa, si|
          if AMINO_ACIDS.include?(aa) || aa == '-'
            aa_frqs[aa] += 1 * @sequences[si].weight
            aa_raw_frqs[aa] += 1
          else
            $logger.warn "#{aa} is a unknown type of amino acid and ignored."
          end
        end

        sum     = aa_frqs.values.sum
        aa_sum  = sum - aa_frqs['-']
        aa_frqs.each do |aa, f|
          aa_prbs[aa]     = f / sum
          aa_rel_prbs[aa] = f / aa_sum if aa != '-'
        end

        prf_pss << SequenceProfilePosition.new(probe, aa_raw_frqs, aa_frqs, aa_prbs, aa_rel_prbs)
      end

      SequenceProfile.new(self, prf_pss)
    end

  end
end
