module Gigue
  class ProfileProfileAlignment

    attr_reader :structural_profile, :sequence_profile,
                :raw_score, :reverse_raw_score,
                :aligned_structural_profile_positions,
                :aligned_sequence_profile_positions

    def initialize(str_prf, seq_prf)
      @structural_profile = str_prf
      @sequence_profile   = seq_prf
    end

    def traceback_linear_gap(max_m=nil, max_n=nil)
      str_pss     = @structural_profile.positions
      str_cnt     = @structural_profile.sequences.size
      seq_pss     = @sequence_profile.positions
      seq_cnt     = @sequence_profile.sequences.size
      ali_str_pss = []
      ali_seq_pss = []
      max_m       = @sequence_profile.length if max_m.nil?
      max_n       = @structural_profile.length if max_n.nil?
      m, n        = max_m, max_n
      log_fmt     = "%-12s : %-5s : %12s"

      loop do
        case @point[m][n]
        when NONE
          break
        when DIAG
          ali_seq_pss << seq_pss[m-1]
          ali_str_pss << str_pss[n-1]
          m -= 1
          n -= 1
          log_mat = "M[#{m}][#{n}]"
          log_prb = "#{ali_str_pss[-1].probe} <=> #{ali_seq_pss[-1].probe}"
          $logger.debug log_fmt % [log_mat, "DIAG", log_prb]
        when LEFT
          ali_seq_pss << SequenceProfilePosition.new('-'*seq_cnt)
          ali_str_pss << str_pss[n-1]
          n -= 1
          log_mat = "M[#{m}][#{n}]"
          log_prb = "#{ali_str_pss[-1].probe} <=> #{ali_seq_pss[-1].probe}"
          $logger.debug log_fmt % [log_mat, "LEFT", log_prb]
        when UP
          ali_seq_pss << seq_pss[m-1]
          ali_str_pss << StructuralProfilePosition.new('-'*str_cnt)
          m -= 1
          log_mat = "M[#{m}][#{n}]"
          log_prb = "#{ali_str_pss[-1].probe} <=> #{ali_seq_pss[-1].probe}"
          $logger.debug log_fmt % [log_mat, " UP ", log_prb]
        else
          $logger.error "Something wrong with pointer of #{log_mat}"
          exit 1
        end
      end

      @aligned_sequence_profile_positions = ali_seq_pss.reverse!
      @aligned_structural_profile_positions = ali_str_pss.reverse!
    end

    def traceback_affine_gap(max_m=nil, max_n=nil)
      str_pss     = @structural_profile.positions
      str_cnt     = @structural_profile.sequences.size
      seq_pss     = @sequence_profile.positions
      seq_cnt     = @sequence_profile.sequences.size
      ali_str_pss = []
      ali_seq_pss = []
      pre_mat     = nil
      max_m       = @sequence_profile.length if max_m.nil?
      max_n       = @structural_profile.length if max_n.nil?
      m, n        = max_m, max_n
      log_fmt     = "%-12s : %-5s : %12s"
      cur_score, cur_point, cur_jump =
        case @raw_score
        when @mat_score[m][n]
          pre_mat = :mat
          [@mat_score, @mat_point, @mat_jump]
        when @del_score[m][n]
          pre_mat = :del
          [@del_score, @del_point, @del_jump]
        when @ins_score[m][n]
          pre_mat = :ins
          [@ins_score, @ins_point, @ins_jump]
        else
          $logger.error "Maximum matrix tag doesn't match with any matrix"
          exit 1
        end

      $logger.debug "START TRACEBACK"

      loop do
        # jumping between matrices
        until cur_jump[m][n].nil?
          case pre_mat
          when :mat
            ali_str_pss << str_pss[n-1]
            ali_seq_pss << seq_pss[m-1]
            log_mat = "M[#{m}][#{n}]"
            log_prb = "#{ali_str_pss[-1].probe} <=> #{ali_seq_pss[-1].probe}"
            log_str = log_fmt % [log_mat, "JUMP", log_prb]
            $logger.debug log_str
          when :del
            ali_str_pss << str_pss[n-1]
            ali_seq_pss << SequenceProfilePosition.new('-'*seq_cnt)
            log_mat = "D[#{m}][#{n}]"
            log_prb = "#{ali_str_pss[-1].probe} <=> #{ali_seq_pss[-1].probe}"
            log_str = log_fmt % [log_mat, "JUMP", log_prb]
            $logger.debug log_str
          when :ins
            ali_str_pss << StructuralProfilePosition.new('-'*str_cnt)
            ali_seq_pss << seq_pss[m-1]
            log_mat = "I[#{m}][#{n}]"
            log_prb = "#{ali_str_pss[-1].probe} <=> #{ali_seq_pss[-1].probe}"
            log_str = log_fmt % [log_mat, "JUMP", log_prb]
            $logger.debug log_str
          end

          mat, m, n = cur_jump[m][n]
          cur_score, cur_point, cur_jump =
            case mat
            when 'M'
              pre_mat = :mat
              [@mat_score, @mat_point, @mat_jump]
            when 'D'
              pre_mat = :del
              [@del_score, @del_point, @del_jump]
            when 'I'
              pre_mat = :ins
              [@ins_score, @ins_point, @ins_jump]
            else
              $logger.warn "Something wrong in jumnng step"
              exit 1
            end
        end

        # following pointers in a matrix
        point = cur_point[m][n]

        case point
        when NONE
          $logger.debug "FINISH TRACEBACK"
          break
        when DIAG
          ali_str_pss << str_pss[n-1]
          ali_seq_pss << seq_pss[m-1]
          m -= 1
          n -= 1
          log_mat = "M[#{m}][#{n}]"
          log_prb = "#{ali_str_pss[-1].probe} <=> #{ali_seq_pss[-1].probe}"
          log_str = log_fmt % [log_mat, "DIAG", log_prb]
          $logger.debug log_str
        when LEFT
          ali_str_pss << str_pss[n-1]
          ali_seq_pss << SequenceProfilePosition.new('-'*seq_cnt)
          n -= 1
          log_mat = "D[#{m}][#{n}]"
          log_prb = "#{ali_str_pss[-1].probe} <=> #{ali_seq_pss[-1].probe}"
          log_str = log_fmt % [log_mat, "LEFT", log_prb]
          $logger.debug log_str
        when UP
          ali_str_pss << StructuralProfilePosition.new('-'*str_cnt)
          ali_seq_pss << seq_pss[m-1]
          m -= 1
          log_mat = "I[#{m}][#{n}]"
          log_prb = "#{ali_str_pss[-1].probe} <=> #{ali_seq_pss[-1].probe}"
          log_str = log_fmt % [log_mat, " UP ", log_prb]
          $logger.debug log_str
        else
          $logger.error "Something wrong with pointing stage at m: #{m}, n: #{n}"
          exit 1
        end
      end

      @aligned_sequence_profile_positions = ali_seq_pss.reverse!
      @aligned_structural_profile_positions = ali_str_pss.reverse!
    end

    # memoize Z-score
    def z_score
      @z_score ||= calculate_z_score
    end

    def reverse_raw_score
      @reverse_raw_score ||= calculate_reverse_raw_score
    end

    def to_flatfile(options={})
      opts = {
        :os     => STDOUT,
        :type   => :pir,
        :width  => 70
      }.merge!(options)

      out = opts[:os].is_a?(String) ? File.open(opts[:os], 'w') : opts[:os]

      # print aligned structural profile sequences
      @structural_profile.sequences.each_with_index do |pseq, psi|
        out.puts ">#{pseq.code}"
        out.puts "structure" if opts[:type] == :pir
        out.puts @aligned_structural_profile_positions.map_with_index { |p, pi|
          ((p == '-' ? '-' : (p.probe[psi] == 'J' ? 'C' : p.probe[psi])) +
           (pi > 0 && (pi+1) % opts[:width] == 0 ? "\n" : ''))
        }.join('')
      end

      # print aligned sequence profile sequences
      @sequence_profile.sequences.each_with_index do |pseq, psi|
        out.puts ">#{pseq.code}"
        out.puts "sequence" if opts[:type] == :pir
        out.puts @aligned_sequence_profile_positions.map_with_index { |p, pi|
          ((p == '-' ? '-' : (p.probe[psi] == 'J' ? 'C' : p.probe[psi])) +
           (pi > 0 && (pi+1) % opts[:width] == 0 ? "\n" : ''))
        }.join('')
      end

      out.close if [File, String].include?(out.class)
    end

  end

  class ProfileProfileLocalAlignmentLinearGap < ProfileProfileAlignment

    attr_reader :score, :point, :max_m, :max_n

    def initialize(str_prf, seq_prf, score, point, max_m, max_n)
      super(str_prf, seq_prf)
      @score      = score
      @point      = point
      @max_m      = max_m
      @max_n      = max_n
      @raw_score  = @score[max_m][max_n]
      traceback_linear_gap(max_m, max_n)
    end

    def calculate_reverse_raw_score
      aligner = ProfileSequenceAligner.new(@structural_profile, @sequence_profile.reverse)
      begin
        aligner.local_alignment_linear_gap_cpp.raw_score
      rescue
        aligner.local_alignment_linear_gap_rb.raw_score
      end
    end

    def calculate_z_score(iter=100)
      scores = NArray.int(iter)
      (0...iter).each do |i|
        aligner   = ProfileProfileAligner.new(@structural_profile, @sequence_profile.shuffle)
        scores[i] = aligner.local_alignment_linear_gap.raw_score
      end
      (@raw_score - scores.mean) / scores.stddev
    end

  end

  class ProfileProfileLocalAlignmentAffineGap < ProfileProfileAlignment

    attr_reader :mat_score, :mat_point, :mat_jump,
                :del_score, :del_point, :del_jump,
                :ins_score, :ins_point, :ins_jump

    def initialize(str_prf, seq_prf,
                   mat_score, mat_point, mat_jump,
                   del_score, del_point, del_jump,
                   ins_score, ins_point, ins_jump,
                   max_m, max_n)
      super(str_prf, seq_prf)
      @mat_score  = mat_score
      @mat_point  = mat_point
      @mat_jump   = mat_jump
      @del_score  = del_score
      @del_point  = del_point
      @del_jump   = del_jump
      @ins_score  = ins_score
      @ins_point  = ins_point
      @ins_jump   = ins_jump
      @raw_score  = [
        @mat_score[max_m][max_n],
        @del_score[max_m][max_n],
        @ins_score[max_m][max_n]
      ].max
      traceback_affine_gap(max_m, max_n)
    end

    def calculate_reverse_raw_score
      aligner = ProfileProfileAligner.new(@structural_profile, @sequence_profile.reverse)
      begin
        aligner.local_alignment_affine_gap_cpp.raw_score
      rescue
        aligner.local_alignment_affine_gap_rb.raw_score
      end
    end

    def calculate_z_score(iter=100)
      scores = NArray.int(iter)
      (0...iter).each do |i|
        aligner   = ProfileProfileAligner.new(@structural_profile, @sequence_profile.shuffle)
        scores[i] = aligner.local_alignment_affine_gap.raw_score
      end
      (@raw_score - scores.mean) / scores.stddev
    end

  end

  class ProfileProfileGlobalAlignmentLinearGap < ProfileProfileAlignment

    attr_reader :score, :point

    def initialize(str_prf, seq_prf, score, point)
      super(str_prf, seq_prf)
      @score      = score
      @point      = point
      @raw_score  = @score[-1][-1]
      traceback_linear_gap
    end

    def calculate_reverse_raw_score
      aligner = ProfileProfileAligner.new(@structural_profile, @sequence_profile.reverse)
      begin
        aligner.global_alignment_linear_gap_cpp.raw_score
      rescue
        aligner.global_alignment_linear_gap_rb.raw_score
      end
    end

    def calculate_z_score(iter=100)
      scores = NArray.int(iter)
      (0...iter).each do |i|
        aligner   = ProfileProfileAligner.new(@structural_profile, @sequence_profile.shuffle)
        scores[i] = aligner.global_alignment_linear_gap.raw_score
      end
      (@raw_score - scores.mean) / scores.stddev
    end

  end

  class ProfileProfileGlobalAlignmentAffineGap < ProfileProfileAlignment

    attr_reader :mat_score, :mat_point, :mat_jump,
                :del_score, :del_point, :del_jump,
                :ins_score, :ins_point, :ins_jump

    def initialize(str_prf, seq_prf,
                   mat_score, mat_point, mat_jump,
                   del_score, del_point, del_jump,
                   ins_score, ins_point, ins_jump)
      super(str_prf, seq_prf)
      @mat_score  = mat_score
      @mat_point  = mat_point
      @mat_jump   = mat_jump
      @del_score  = del_score
      @del_point  = del_point
      @del_jump   = del_jump
      @ins_score  = ins_score
      @ins_point  = ins_point
      @ins_jump   = ins_jump
      @raw_score  = [
        @mat_score[-1][-1],
        @del_score[-1][-1],
        @ins_score[-1][-1]
      ].max
      traceback_affine_gap
    end

    def calculate_reverse_raw_score
      aligner = ProfileProfileAligner.new(@structural_profile, @sequence_profile.reverse)
      begin
        aligner.global_alignment_affine_gap_cpp.raw_score
      rescue
        aligner.global_alignment_affine_gap_rb.raw_score
      end
    end

    def calculate_z_score(iter=100)
      scores = NArray.int(iter)
      (0...iter).each do |i|
        aligner   = ProfileProfileAligner.new(@structural_profile, @sequence_profile.shuffle)
        scores[i] = aligner.global_alignment_affine_gap.raw_score
      end
      (@raw_score - scores.mean) / scores.stddev
    end

  end
end
