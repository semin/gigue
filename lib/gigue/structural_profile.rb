module Gigue
  class StructuralProfile

    # Default gap penalties from Shi et al. JMB (2001)
    VVLi, VLi, Li, Hi = 8, 20, 22, 28
    VVLe, VLe, Le, He = 2, 2, 4, 8

    attr_reader :joytem, :essts, :weights, :seq_cnt,
                :sequences, :entry_names, :length, :positions

    def initialize(joy, essts, opts = {})
      @joytem       = JoyTem.new(joy)
      @essts        = Essts.new(essts)
      @sequences    = @joytem.sequences
      @length       = @joytem.alignment_length
      @seq_cnt      = @sequences.size
      @entry_names  = @sequences.map(&:code)
      @options      = {
        :weighting          => :blosum,
        :multi              => 10,
        :ignore_gap_weight  => true,
        :gap_ins_open_term  => VVLi,
        :gap_ins_ext_term   => VVLe,
        :gap_del_open_term  => VVLi,
        :gap_del_ext_term   => VVLe
      }.merge!(opts)

      #
      # 1. Weighting structures
      #
      if    @options[:weighting] == :blosum
        @weights = calculate_blosum_weights
      elsif @options[:weighting] == :va
        @weights = calculate_va_weights
      end

      #
      # 2. Annotate structures with ESSTs and JOY
      #
      @sequences.each do |seq|
        seq.environments = Array.new(@length, '')
        @essts.environments.each do |env|
          next if env.silent
          values = @joytem.entries[seq.code][env.name].split('')
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

      #
      # 3. Generating PSSMs
      #
      @positions = []

      (0...@length).each do |pi|
        probe       = ''
        score       = NMatrix.int(1, @essts.rownames.size).fill!(0)
        gap_score   = Hash.new(0)
        @sequences.each do |seq|
          seq.gap_cnt = 0 if pi == 0
          aa          = seq.data[pi]
          env         = seq.environments[pi]
          probe       += aa
          if aa == '-'
            seq.gap_cnt += 1
          else
            # 3.1. Scoring matrix
            #
            score += @options[:multi] * @weights[seq.code] * @essts[env].scores_from(aa)
            # 3.2. Gap penalties
            #
            ui = pi - seq.gap_cnt # ungapped index
            # 3.2.1. Gap penalties for terminal region
            #
            if (ui == 0) || (ui == seq.ungapped_length-1)
              gap_score['DelO'] = @options[:multi] * @weights[seq.code] * @options[:gap_del_open_term]
              gap_score['DelE'] = @options[:multi] * @weights[seq.code] * @options[:gap_del_ext_term]
              gap_score['InsO'] = @options[:multi] * @weights[seq.code] * @options[:gap_ins_open_term]
              gap_score['InsE'] = @options[:multi] * @weights[seq.code] * @options[:gap_ins_ext_term]
            else
              # 3.2.2. Gap penalties for deletion
              #
              if    ((seq.ungapped_environments[ui]   =~ /^C/) &&       # i:    COIL
                     (seq.ungapped_environments[ui+1] =~ /^[H|E|P]/))   # i+1:  SSE
                gap_score['DelO'] = @options[:multi] * @weights[seq.code] * VLi
                gap_score['DelE'] = @options[:multi] * @weights[seq.code] * VLe
              elsif ((seq.ungapped_environments[ui-1] =~ /^C/) &&       # i-1:  COIL
                     (seq.ungapped_environments[ui]   =~ /^[H|E|P]/))   # i:    SSE
                gap_score['DelO'] = @options[:multi] * @weights[seq.code] * Li
                gap_score['DelE'] = @options[:multi] * @weights[seq.code] * Le
              elsif ((seq.ungapped_environments[ui-1] =~ /^[H|E|P]/) && # i-1:  SSE
                     (seq.ungapped_environments[ui]   =~ /^[H|E|P]/) && # i:    SSE
                     (seq.ungapped_environments[ui+1] =~ /^[H|E|P]/))   # i+1:  SSE
                gap_score['DelO'] = @options[:multi] * @weights[seq.code] * Hi
                gap_score['DelE'] = @options[:multi] * @weights[seq.code] * He
              elsif ((seq.ungapped_environments[ui]   =~ /^[H|E|P]/) && # i:    SSE
                     (seq.ungapped_environments[ui+1] =~ /^C/))         # i+1:  COIL
                gap_score['DelO'] = @options[:multi] * @weights[seq.code] * Li
                gap_score['DelE'] = @options[:multi] * @weights[seq.code] * Le
              elsif ((seq.ungapped_environments[ui-1] =~ /^[H|E|P]/) && # i-1:  SSE
                     (seq.ungapped_environments[ui]   =~ /^C/))         # i+1:  COIL
                gap_score['DelO'] = @options[:multi] * @weights[seq.code] * VLi
                gap_score['DelE'] = @options[:multi] * @weights[seq.code] * VLe
              elsif ((seq.ungapped_environments[ui-1] =~ /^C/) &&       # i-1:  COIL
                     (seq.ungapped_environments[ui]   =~ /^C/) &&       # i:    COIL
                     (seq.ungapped_environments[ui+1] =~ /^C/))         # i+1:  COIL
                gap_score['DelO'] = @options[:multi] * @weights[seq.code] * VLi
                gap_score['DelE'] = @options[:multi] * @weights[seq.code] * VLe
              else
                $logger.error [
                  "Unknown gap deletion category:",
                   "#{ui-1}: #{ungapped_environments[ui-1]},",
                   "#{ui  }: #{ungapped_environments[ui  ]},",
                   "#{ui+1}: #{ungapped_environments[ui+1]}",
                   "in #{seq.code}"
                ].join(' ')
                exit 1
              end
              # 3.2.3. Gap penalties for insertion
              #
              if    ((seq.ungapped_environments[ui]   =~ /^C/) &&       # i:    COIL
                     (seq.ungapped_environments[ui+1] =~ /^[H|E|P]/))   # i+1:  SSE
                gap_score['InsO'] = @options[:multi] * @weights[seq.code] * VLi
                gap_score['InsE'] = @options[:multi] * @weights[seq.code] * VLe
              elsif ((seq.ungapped_environments[ui]   =~ /^[H|E|P]/) && # i:    SSE
                     (seq.ungapped_environments[ui+1] =~ /^[H|E|P]/))   # i+1:  SSE
                gap_score['InsO'] = @options[:multi] * @weights[seq.code] * Hi
                gap_score['InsE'] = @options[:multi] * @weights[seq.code] * He
              elsif ((seq.ungapped_environments[ui]   =~ /^[H|E|P]/) && # i:    SSE
                     (seq.ungapped_environments[ui+1] =~ /^C/))         # i+1:  COIL
                gap_score['InsO'] = @options[:multi] * @weights[seq.code] * VLi
                gap_score['InsE'] = @options[:multi] * @weights[seq.code] * VLe
              elsif ((seq.ungapped_environments[ui]   =~ /^C/) &&       # i:    COIL
                     (seq.ungapped_environments[ui+1] =~ /^C/))         # i+1:  COIL
                gap_score['InsO'] = @options[:multi] * @weights[seq.code] * VLi
                gap_score['InsE'] = @options[:multi] * @weights[seq.code] * VLe
              else
                $logger.error [
                  "Unknown gap insertion category:",
                   "#{ui-1}: #{ungapped_environments[ui-1]},",
                   "#{ui  }: #{ungapped_environments[ui  ]},",
                   "#{ui+1}: #{ungapped_environments[ui+1]}",
                   "in #{seq.code}"
                ].join(' ')
                exit 1
              end
            end
          end
        end
        mat_score = @essts.rownames.to_hash(score.round.flatten.to_a)
        gap_score.each { |k, v| gap_score[k] = v.round }
        @positions << StructuralProfilePosition.new(probe, mat_score, gap_score)
      end
    end

    def calculate_pid(seq1, seq2)
      aas1  = seq1.split('')
      aas2  = seq2.split('')
      cols  = aas1.zip(aas2)
      gap   = '-'
      align = 0 # no. of aligned columns
      ident = 0 # no. of identical columns
      intgp = 0 # no. of internal gaps
      cols.each do |col|
        if (col[0] != gap) && (col[1] != gap)
          align += 1
          if col[0] == col[1]
            ident += 1
          end
        elsif (((col[0] == gap) && (col[1] != gap)) ||
                ((col[0] != gap) && (col[1] == gap)))
          intgp += 1
        end
      end
      pid = 100.0 * Float(ident) / (align + intgp)
    end

    def calculate_blosum_weights(weight = 60)
      weights   = {}
      clusters  = @sequences.map { |s| [s] }
      begin
        continue = false
        0.upto(clusters.size-2) do |i|
          indexes = []
          (i+1).upto(clusters.size-1) do |j|
            found = false
            clusters[i].each do |s1|
              clusters[j].each do |s2|
                if calculate_pid(s1.data, s2.data) >= weight
                  indexes << j
                  found = true
                  break
                end
              end
              break if found
            end
          end
          unless indexes.empty?
            continue  = true
            group     = clusters[i]
            indexes.each do |k|
              group       = group.concat(clusters[k])
              clusters[k] = nil
            end
            clusters[i] = group
            clusters.compact!
          end
        end
      end while(continue)

      clusters.each do |cluster|
        cluster.each do |seq|
          weight = cluster.size / Float(@sequences.size)
          weights[seq.code] = weight
        end
      end
      weights
    end

    def calculate_va_weights
      $logger.warn "To be implemented"
      exit 1
    end

    def to_fugue_profile(os = STDOUT, type = :fugue)
      os = File.open(os, 'w') if os.is_a? String
      os.puts "%-30s %-s" % ['Command:', "gigue build -t #{@joytem.file} -s #{@essts.file}"]
      os.puts
      os.puts "%-30s %-5d %s" % ['Profile_length:', "#{@length}", 'alignment positions']
      os.puts "%-30s %-5d %s" % ['Sequence_in_profile:', "#{@sequences.size}", 'sequences']
      os.puts "%-30s %-5d %s" % ['Real_Structure:', "#{@sequences.size}", 'structures']
      os.puts "%-30s %-5d %s" % ['Enhance_Num:', "#{@sequences.size}", 'sequences'] # ???
      os.puts "%-30s %-d"     % ['Enhance_Div:', 0] # ???
      os.puts
      os.puts "%-30s %-5d %s" % ['Weighting:', 1, 'BlosumWeight -- weighting scheme based on single linkage clustering']
      os.puts "%-30s %-d"     % ['Weighting_threshold:', 0]
      os.puts "%-30s %-d"     % ['Weighting_seed:', 0]
      os.puts
      for s in @sequences
        os.puts "%-30s %-f"   % [' ' * 10 + s.code, @weights[s.code]]
      end
      os.puts
      os.puts "%-30s %-d"     % ['Multiple_factor:', @options[:multi]]
      os.puts "%-30s %-5d %s" % ['Profile_format:', 0, 'FUGUE']
      os.puts "%-30s %-5s"    % ['Similarity_matrix:', 'OFF'] # ???
      os.puts "%-30s %-5s"    % ['Similarity_matrix_offset:', 'OFF'] # ???
      os.puts "%-30s %-5s"    % ['Ignore_gap_weight:', @options[:ignore_gap_weight] == true ? 'ON' : 'OFF']
      os.puts
      os.puts "%-30s %-5s"    % ['Symbol_in_row(sequence):', @essts.rownames.join]
      os.puts "%-30s %-5s"    % ['Symbol_in_column(structure):', @essts.colnames.join]
      os.puts "%-30s %-5s"    % ['Symbol_structural_feature:', @essts.environments.map { |e| e.labels.join if !e.silent }.join]
      os.puts
      os.puts "%-30s %-5d"    % ['GapInsOpenTerminal',  @options[:gap_ins_open_term]]
      os.puts "%-30s %-5d"    % ['GapDelOpenTerminal',  @options[:gap_del_open_term]]
      os.puts "%-30s %-5d"    % ['GapInsExtTerminal',   @options[:gap_ins_ext_term]]
      os.puts "%-30s %-5d"    % ['GapDelExtTerminal',   @options[:gap_del_ext_term]]
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
