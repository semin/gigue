class FugueProfile

  attr_reader :command, :length, :seq_cnt, :weighting,
              :entry_names, :entry_weights, :multiple_factor,
              :format, :row_symbols, :col_symbols, :env_symbols,
              :gap_open_ins_term, :gap_open_del_term,
              :gap_ext_ins_term, :gap_ext_del_term,
              :aa_colnames, :gap_colnames, :env_colnames,
              :probes, :mat_scores, :gap_scores, :env_scores


  def initialize(file)
    @entry_names    = []
    @entry_weights  = []
    @probes         = []
    @mat_scores     = []
    @gap_scores     = []
    @env_scores     = []
    parse_tag       = nil
    IO.foreach(file) do |line|
      if    line =~ /^Command:\s+(.*)/
        @command = $1.strip
      elsif line =~ /^Profile_length:\s+(\d+)/
        @length = Integer($1)
      elsif line =~ /^Sequence_in_profile:\s+(\d+)/
        @seq_cnt = Integer($1)
      elsif line =~ /^Weighting:\s+(\d+)/
        @weighting  = 1
        parse_tag = :weighting
      elsif line =~ /^\s+(\S+)\s+(\S+)\s*$/ && parse_tag == :weighting
        @entry_names    << $1
        @entry_weights  << Float($2)
      elsif line =~ /^Multiple_factor:\s+(\S+)/
        @multiple_factor = Float($1)
      elsif line =~ /^Profile_format:\s+(\d+)\s+(\S+)/
        @format = "#{$1}-#{$2}"
      elsif line =~ /^Symbol_in_row\(sequence\):\s+(\S+)/ # for sequence
        @row_symbols = $1.split('')
      elsif line =~ /^Symbol_in_column\(structure\):\s+(\S+)/ # for structure
        @col_symbols = $1.split('')
      elsif line =~ /^Symbol_structural_feature:\s+(\S+)/
        @env_symbols = $1.split('')
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
        cols        = line.chomp.split(/\s+/)
        @probes     << cols[0]
        @mat_scores << cols[1,@aa_colnames.size]
        @gap_scores << cols[1+@aa_colnames.size,@gap_colnames.size]
        @env_scores << cols[1+@aa_colnames.size+@gap_colnames.size,@env_colnames.size]
      end
    end
  end
end
