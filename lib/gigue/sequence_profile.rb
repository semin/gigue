module Gigue
  class SequenceProfile

    def self.create_from_multiple_sequence_alignment(msa, weighting=:va)
      if    weighting == :va
        Sequence::calculate_va_weights_cpp(msa.sequences)
      elsif weighting == :blosum
        Sequence::calculate_blosum_weights(msa.sequences)
      elsif weighting == :equal
        Sequence::calculate_equal_weights(msa.sequences)
      else
        $logger.error "Unknown weighting scheme!"
        exit 1
      end

      prf_pss = []

      msa.columns.each do |col|
        probe       = col.probe
        aa_raw_frqs = Hash.new(0)
        aa_frqs     = Hash.new(0.0)
        aa_prbs     = Hash.new(0.0)
        aa_rel_prbs = Hash.new(0.0)

        probe.split('').each_with_index do |aa, si|
          aa_frqs[aa] += 1 * msa.sequences[si].weight
          aa_raw_frqs[aa] += 1
        end

        sum     = aa_frqs.values.sum
        aa_sum  = sum - aa_frqs['-']
        aa_frqs.each do |aa, f|
          aa_prbs[aa]     = f / sum
          aa_rel_prbs[aa] = f / aa_sum if aa != '-'
        end

        prf_pss << SequenceProfilePosition.new(probe, aa_raw_frqs, aa_frqs, aa_prbs, aa_rel_prbs)
      end

      inp_opts = { :msa => msa, :weighting => :va }
      self.new(prf_pss, inp_opts)
    end

    attr_reader :positions, :depth, :length, :sequences

    def initialize(pss, options={})
      @opts       = options
      @positions  = pss
      @depth      = pss[0].probe.length
      @length     = pss.size
      @msa        = @opts[:msa]
      @sequences  = @msa.sequences
    end

    #def sequences
      ##seqs = Array.new(depth, '')
      ##(0...length).each do |pi|
        ##(0...depth).each do |si|
          ##seqs[si] += @positions[pi].probe[si]
        ##end
      ##end
      ##seqs
      #@msa.sequences
    #end

    #def depth
      #@positions[0].probe.length
    #end

    def shuffle
      self.class.new(@positions.shuffle, @opts)
    end

    def reverse
      self.class.new(@positions.reverse, @opts)
    end

    def purge
      pss = []
      @positions.each do |ps|
        w = 0
        ps.probe.split('').each_with_index do |aa, si|
          w += sequences[si].weight if aa == '-'
        end
        pss << ps if w < 0.5
      end
      self.class.new(pss, @opts)
    end

    def to_fug(os = STDOUT)
      os = File.open(os, 'w') if os.is_a? String
      os.puts "%-30s %-s" % ['Command:', "gigue"]
      os.puts
      os.puts "%-30s %-5d %s" % ['Profile_length:', "#{@length}", 'alignment positions']
      os.puts "%-30s %-5d %s" % ['Sequence_in_profile:', "#{@depth}", 'sequences']
      os.puts "%-30s %-5d %s" % ['Real_Structure:', 0, 'structures']
      os.puts "%-30s %-5d %s" % ['Enhance_Num:', "#{@depth}", 'sequences'] # ???
      os.puts "%-30s %-d"     % ['Enhance_Div:', 0] # ???
      os.puts
      if    @opts[:weighting] == :va
        os.puts "%-30s %-5d %s" % ['Weighting:', 2, 'VAWeight -- Vingon and Argo weight']
      elsif @opts[:weighting] == :blosum
        os.puts "%-30s %-5d %s" % ['Weighting:', 1, 'BlosumWeight -- weighting scheme based on single linkage clustering']
      elsif @opts[:weighting] == :equal
        os.puts "%-30s %-5d %s" % ['Weighting:', 0, 'EqualWeight -- weight each sequence equally']
      else
        $logger.error "Unknown weighting scheme!"
      end
      os.puts "%-30s %-d"     % ['Weighting_threshold:', 0]
      os.puts "%-30s %-d"     % ['Weighting_seed:', 0]
      os.puts
      for s in @sequences
        os.puts "%-30s %-f"   % [' ' * 10 + s.code, s.weight]
      end
      os.puts
      os.puts "%-30s %-5d %s" % ['Profile_format:', 0, 'FUGUE']
      os.puts "%-30s %-5s"    % ['Similarity_matrix:', 'OFF'] # ???
      os.puts "%-30s %-5s"    % ['Similarity_matrix_offset:', 'OFF'] # ???
      os.puts
      os.puts "%-30s %-5s"    % ['Symbol_in_row(sequence):', AMINO_ACIDS]
      os.puts "%-30s %-5s"    % ['Symbol_in_column(sequence):', AMINO_ACIDS]
      os.puts "\n\n\n"
      format = "%-#{depth}s/" + ("%5s " * AMINO_ACIDS.size)
      os.puts format % [' Seq', *AMINO_ACIDS.split('')]
      os.puts
      os.puts "START"
      format = "%-#{depth}s " + ("%5.3f " * AMINO_ACIDS.size)
      for ps in @positions
        os.puts format % [ps.probe, *AMINO_ACIDS.split('').map { |a| ps.relative_probability_of(a) }]
      end
      os.puts "THEEND"
      os.close if [File, String].include? os.class
    end

  end
end
