module Gigue
  class MultipleSequenceAlignment

    def self.create_from_psiblast_output_style6(file)
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
      self.new(seqs)
    end


    attr_reader :sequences, :length, :depth, :columns

    def initialize(seqs)
      @sequences  = seqs
      @depth      = seqs.size
      @length     = seqs[0].length
      @columns    = (0...@length).map do |mi|
        MultipleSequenceAlignmentColumn.new(@sequences.inject('') { |p, s| p + s.data[mi] })
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
      SequenceProfile::create_from_multiple_sequence_alignment(self)
    end

  end
end
