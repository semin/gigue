module Gigue
  class FugueProfile

    attr_reader :command, :length, :weighting,
                :sequences, :multiple_factor, :format,
                :rowsymbols, :colsymbols, :envsymbols,
                :gap_open_ins_term, :gap_open_del_term,
                :gap_ext_ins_term, :gap_ext_del_term,
                :aa_colnames, :gap_colnames, :env_colnames,
                :positions

    def initialize(file)
      @sequences  = []
      @positions  = []
      parse_tag   = nil

      IO.foreach(file) do |line|
        if    line =~ /^Command:\s+(.*)/
          @command = $1.strip
        elsif line =~ /^Profile_length:\s+(\d+)/
          @length  = Integer($1)
        elsif line =~ /^Sequence_in_profile:\s+(\d+)/
          @seq_cnt = Integer($1)
        elsif line =~ /^Weighting:\s+(\d+)/
          @weighting  = 1
          parse_tag   = :weighting
        elsif line =~ /^\s+(\S+)\s+(\S+)\s*$/ && parse_tag == :weighting
          @sequences << Sequence.new('', $1)
          @sequences[-1].weight = Float($2)
        elsif line =~ /^Multiple_factor:\s+(\S+)/
          @multiple_factor = Float($1)
        elsif line =~ /^Profile_format:\s+(\d+)\s+(\S+)/
          @format = "#{$1}-#{$2}"
        elsif line =~ /^Symbol_in_row\(sequence\):\s+(\S+)/ # for sequence
          @rowsymbols = $1.split('')
        elsif line =~ /^Symbol_in_column\(structure\):\s+(\S+)/ # for structure
          @colsymbols = $1.split('')
        elsif line =~ /^Symbol_structural_feature:\s+(\S+)/
          @envsymbols = $1.split('')
        elsif line =~ /^GapInsOpenTerminal\s+(\S+)/
          @gap_open_ins_term = Float($1)
        elsif line =~ /^GapDelOpenTerminal\s+(\S+)/
          @gap_open_del_term = Float($1)
        elsif line =~ /^GapInsExtTerminal\s+(\S+)/
          @gap_ext_ins_term = Float($1)
        elsif line =~ /^GapDelExtTerminal\s+(\S+)/
          @gap_ext_del_term = Float($1)
        elsif line =~ /^\s+Seq\s+\/(.*)/
          colnames      = $1.split('/')
          @aa_colnames  = colnames[0].strip.split(/\s+/)
          @gap_colnames = colnames[1].strip.split(/\s+/)
          @env_colnames = colnames[2].strip.split(/\s+/)
        elsif line =~ /^START/
          parse_tag = :profile
        elsif line =~ /^THEEND/
          break
        elsif parse_tag == :profile
          cols      = line.chomp.split(/\s+/)
          probe     = cols[0]
          mat_score = cols[1,@aa_colnames.size].map { |s| Integer(s) }
          gap_score = cols[1+@aa_colnames.size,@gap_colnames.size].map { |s| Integer(s) }
          env_score = cols[1+@aa_colnames.size+@gap_colnames.size,@env_colnames.size].map { |s| Integer(s) }
          @positions << FugueProfilePosition.new(
            probe,
            Hash[*@aa_colnames.zip(mat_score).flatten],
            Hash[*@gap_colnames.zip(gap_score).flatten],
            Hash[*@env_colnames.zip(env_score).flatten]
          )
          @sequences.each_with_index { |s, i| s.data += probe[i] }
        end
      end
    end

  end
end
