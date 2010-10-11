module Gigue
  class StructuralProfile

    # Default gap penalties from Shi et al. JMB (2001)
    VVLi, VLi, Li, Hi = 8, 20, 22, 28
    VVLe, VLe, Le, He = 2, 2, 4, 8

    attr_reader :joytem, :essts, :sequences,
                :length, :depth, :positions,
                :amino_acids

    def self.create_from_joy_tem_and_essts(joy, essts, options={})
      joytem    = JoyTem.new(joy)
      essts     = Essts.new(essts)
      sequences = joytem.sequences
      length    = joytem.alignment_length
      depth     = sequences.size
      positions = []
      opts      = {
        :weighting          => :va,
        :multi              => 10,
        :ignore_gap_weight  => true,
        :gap_ins_open_term  => VVLi,
        :gap_ins_ext_term   => VVLe,
        :gap_del_open_term  => VVLi,
        :gap_del_ext_term   => VVLe
      }.merge!(options)

      # 1. Weighting structures
      if    opts[:weighting] == :va
        Sequence::calculate_va_weights(sequences)
      elsif opts[:weighting] == :blosum
        Sequence::calculate_blosum_weights(sequences)
      elsif opts[:weighting] == :equal
        Sequence::calculate_equal_weights(sequences)
      end

      # 2. Annotate structures with ESSTs and JOY
      sequences.each do |seq|
        seq.environments = Array.new(length, '')
        essts.environments.each do |env|
          next if env.silent
          values = joytem.entries[seq.code][env.name].split('')
          values.each_with_index do |value, vi|
            if value == '-'
              seq.environments[vi] += '-'
            else
              seq.environments[vi] += env.labels[env.values.index(value)]
            end
          end
        end
        seq.ungapped_data         = seq.data.gsub('-', '')
        seq.ungapped_length       = seq.ungapped_data.length
        seq.ungapped_environments = seq.environments.select { |e| e !~ /-/ }
      end

      # 3. Generating PSSMs
      (0...length).each do |pi|
        probe       = ''
        score       = NMatrix.int(1, essts.rownames.size).fill!(0)
        gap_score   = Hash.new(0)
        sequences.each do |seq|
          seq.gap_cnt = 0 if pi == 0
          aa          = seq.data[pi]
          env         = seq.environments[pi]
          probe       += aa
          if aa == '-'
            seq.gap_cnt += 1
          else
            # 3.1. Scoring matrix
            score += opts[:multi] * seq.weight * essts[env].scores_from(aa)

            # 3.2. Gap penalties
            ui = pi - seq.gap_cnt # ungapped index

            # 3.2.1 Gap penalties for the first and the last column
            if (ui == 0) || (ui == seq.ungapped_length-1)
              gap_score['DelO'] += opts[:multi] * seq.weight * opts[:gap_del_open_term]
              gap_score['DelE'] += opts[:multi] * seq.weight * opts[:gap_del_ext_term]
              gap_score['InsO'] += opts[:multi] * seq.weight * opts[:gap_ins_open_term]
              gap_score['InsE'] += opts[:multi] * seq.weight * opts[:gap_ins_ext_term]
              next
            end

            # 3.2.2. Gap penalties for deletion
            if    ((seq.ungapped_environments[ui]   =~ /^C/) &&       # i:    COIL
                   (seq.ungapped_environments[ui+1] =~ /^[H|E|P]/))   # i+1:  SSE
              gap_score['DelO'] += opts[:multi] * seq.weight * VLi
              gap_score['DelE'] += opts[:multi] * seq.weight * VLe
            elsif ((seq.ungapped_environments[ui-1] =~ /^C/) &&       # i-1:  COIL
                   (seq.ungapped_environments[ui]   =~ /^[H|E|P]/))   # i:    SSE
              gap_score['DelO'] += opts[:multi] * seq.weight * Li
              gap_score['DelE'] += opts[:multi] * seq.weight * Le
            elsif ((seq.ungapped_environments[ui-1] =~ /^[H|E|P]/) && # i-1:  SSE
                   (seq.ungapped_environments[ui]   =~ /^[H|E|P]/) && # i:    SSE
                   (seq.ungapped_environments[ui+1] =~ /^[H|E|P]/))   # i+1:  SSE
              gap_score['DelO'] += opts[:multi] * seq.weight * Hi
              gap_score['DelE'] += opts[:multi] * seq.weight * He
            elsif ((seq.ungapped_environments[ui]   =~ /^[H|E|P]/) && # i:    SSE
                   (seq.ungapped_environments[ui+1] =~ /^C/))         # i+1:  COIL
              gap_score['DelO'] += opts[:multi] * seq.weight * Li
              gap_score['DelE'] += opts[:multi] * seq.weight * Le
            elsif ((seq.ungapped_environments[ui-1] =~ /^[H|E|P]/) && # i-1:  SSE
                   (seq.ungapped_environments[ui]   =~ /^C/))         # i+1:  COIL
              gap_score['DelO'] += opts[:multi] * seq.weight * VLi
              gap_score['DelE'] += opts[:multi] * seq.weight * VLe
            elsif ((seq.ungapped_environments[ui-1] =~ /^C/) &&       # i-1:  COIL
                   (seq.ungapped_environments[ui]   =~ /^C/) &&       # i:    COIL
                   (seq.ungapped_environments[ui+1] =~ /^C/))         # i+1:  COIL
              gap_score['DelO'] += opts[:multi] * seq.weight * VLi
              gap_score['DelE'] += opts[:multi] * seq.weight * VLe
            else
              $logger.error [
                "Unknown gap deletion category:",
                  "#{ui-1}: #{seq.ungapped_environments[ui-1]},",
                  "#{ui  }: #{seq.ungapped_environments[ui  ]},",
                  "#{ui+1}: #{seq.ungapped_environments[ui+1]}",
                  "in #{seq.code}"
              ].join(' ')
              exit 1
            end

            # 3.2.3. Gap penalties for insertion
            if    ((seq.ungapped_environments[ui]   =~ /^C/) &&       # i:    COIL
                   (seq.ungapped_environments[ui+1] =~ /^[H|E|P]/))   # i+1:  SSE
              gap_score['InsO'] += opts[:multi] * seq.weight * VLi
              gap_score['InsE'] += opts[:multi] * seq.weight * VLe
            elsif ((seq.ungapped_environments[ui]   =~ /^[H|E|P]/) && # i:    SSE
                   (seq.ungapped_environments[ui+1] =~ /^[H|E|P]/))   # i+1:  SSE
              gap_score['InsO'] += opts[:multi] * seq.weight * Hi
              gap_score['InsE'] += opts[:multi] * seq.weight * He
            elsif ((seq.ungapped_environments[ui]   =~ /^[H|E|P]/) && # i:    SSE
                   (seq.ungapped_environments[ui+1] =~ /^C/))         # i+1:  COIL
              gap_score['InsO'] += opts[:multi] * seq.weight * VLi
              gap_score['InsE'] += opts[:multi] * seq.weight * VLe
            elsif ((seq.ungapped_environments[ui]   =~ /^C/) &&       # i:    COIL
                   (seq.ungapped_environments[ui+1] =~ /^C/))         # i+1:  COIL
              gap_score['InsO'] += opts[:multi] * seq.weight * VLi
              gap_score['InsE'] += opts[:multi] * seq.weight * VLe
            else
              $logger.error [
                "Unknown gap insertion category:",
                  "#{ui-1}: #{seq.ungapped_environments[ui-1]},",
                  "#{ui  }: #{seq.ungapped_environments[ui  ]},",
                  "#{ui+1}: #{seq.ungapped_environments[ui+1]}",
                  "in #{seq.code}"
              ].join(' ')
              exit 1
            end
          end
        end
        mat_score = essts.rownames.to_hash(score.round.flatten.to_a)
        gap_score.each { |k, v| gap_score[k] = v.round }
        positions << StructuralProfilePosition.new(probe, mat_score, gap_score)
      end
      add_opts = {
        :joytem    => joytem,
        :essts     => essts,
        :sequences => sequences,
      }
      self.new(positions, opts.merge!(add_opts))
    end

    def initialize(pss, options={})
      @opts       = {}.merge!(options)
      @positions  = pss
      @joytem     = @opts[:joytem]
      @essts      = @opts[:essts]
      @sequences  = @opts[:sequences]
      @length     = @positions.length
      @depth      = @sequences.size
    end

    def purge(cutoff=0.5)
      pss = []
      @positions.each do |ps|
        w = 0
        ps.probe.split('').each_with_index do |aa, si|
          w += @sequences[si].weight if aa == '-'
        end
        pss << ps if w < cutoff
      end
      self.class.new(pss, @opts)
    end

    def to_gig(os = STDOUT)
      os = File.open(os, 'w') if os.is_a? String
      os.puts "%-30s %-s" % ['Command:', "gigue build -t #{@joytem.file} -s #{@essts.file}"]
      os.puts
      os.puts "%-30s %-5d %s" % ['Profile_length:', "#{@length}", 'alignment positions']
      os.puts "%-30s %-5d %s" % ['Sequence_in_profile:', "#{@sequences.size}", 'sequences']
      os.puts "%-30s %-5d %s" % ['Real_Structure:', "#{@sequences.size}", 'structures']
      os.puts "%-30s %-5d %s" % ['Enhance_Num:', "#{@sequences.size}", 'sequences'] # ???
      os.puts "%-30s %-d"     % ['Enhance_Div:', 0] # ???
      os.puts
      if    @opts[:weighting] == :va
        os.puts "%-30s %-5d %s" % ['Weighting:', 2, 'VAWeight -- Vingron and Argo weight']
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
      os.puts "%-30s %-d"     % ['Multiple_factor:', @opts[:multi]]
      os.puts "%-30s %-5d %s" % ['Profile_format:', 0, 'FUGUE']
      os.puts "%-30s %-5s"    % ['Similarity_matrix:', 'OFF'] # ???
      os.puts "%-30s %-5s"    % ['Similarity_matrix_offset:', 'OFF'] # ???
      os.puts "%-30s %-5s"    % ['Ignore_gap_weight:', @opts[:ignore_gap_weight] == true ? 'ON' : 'OFF']
      os.puts
      os.puts "%-30s %-5s"    % ['Symbol_in_row(sequence):', @essts.rownames.join]
      os.puts "%-30s %-5s"    % ['Symbol_in_column(structure):', @essts.colnames.join]
      os.puts "%-30s %-5s"    % ['Symbol_structural_feature:', @essts.environments.map { |e| e.labels.join if !e.silent }.join]
      os.puts
      os.puts "%-30s %-5d"    % ['GapInsOpenTerminal',  @opts[:multi] * @opts[:gap_ins_open_term]]
      os.puts "%-30s %-5d"    % ['GapDelOpenTerminal',  @opts[:multi] * @opts[:gap_del_open_term]]
      os.puts "%-30s %-5d"    % ['GapInsExtTerminal',   @opts[:multi] * @opts[:gap_ins_ext_term]]
      os.puts "%-30s %-5d"    % ['GapDelExtTerminal',   @opts[:multi] * @opts[:gap_del_ext_term]]
      os.puts
      os.puts "%-30s %-5d"    % ['EVD', 0] # ???
      os.puts "\n\n\n"
      format = "%-#{@sequences.size}s/" + "%5s" * @essts.rownames.size + " / " + "%-5s" * 4
      os.puts format % [' Seq', *@essts.rownames, 'InsO', 'InsE', 'DelO', 'DelE']
      os.puts
      os.puts "START"
      format = "%-#{@sequences.size}s " + "%5s" * @essts.rownames.size + "   " + "%-5s" * 4
      for ps in @positions
        os.puts format % [ps.probe, *@essts.rownames.map { |a| ps.mat_score(a) }.flatten, ps.gap_ins_open, ps.gap_ins_ext, ps.gap_del_open, ps.gap_del_ext]
      end
      os.puts "THEEND"
      os.close if [File, String].include? os.class
    end

  end
end
