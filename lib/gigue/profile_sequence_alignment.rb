module Gigue
  class ProfileSequenceAlignment

    attr_reader :structural_profile, :sequence,
                :raw_score, :reverse_score,
                :aligned_structural_profile_positions,
                :aligned_amino_acids

    def initialize(prf, seq)
      @structural_profile = prf
      @sequence           = seq
    end

    # memoize Z-score
    def z_score
      @z_score ||= calculate_z_score
    end

    def reverse_score
      @reverse_score ||= calculate_reverse_score
    end

    def to_flatfile(options={})
      opts = {
        :os     => STDOUT,
        :type   => :pir,
        :width  => 70
      }.merge!(options)

      out = opts[:os].is_a?(String) ? File.open(opts[:os], 'w') : opts[:os]

      # 1. print aligned structural profile sequences
      @structural_profile.sequences.each_with_index do |pseq, psi|
        out.puts ">#{pseq.code}"
        out.puts "structure" if opts[:type] == :pir
        out.puts @aligned_structural_profile_positions.map_with_index { |p, pi|
          ((p == '-' ? '-' : (p.probe[psi] == 'J' ? 'C' : p.probe[psi])) +
           (pi > 0 && (pi+1) % opts[:width] == 0 ? "\n" : ''))
        }.join('')
      end

      # 2. print aligned query sequence
      out.puts ">#{@sequence.code}"
      out.puts "sequence" if opts[:type] == :pir
      out.puts @aligned_amino_acids.map_with_index { |a, ai|
        a + (ai > 0 && (ai+1) % opts[:width] == 0 ? "\n" : '')
      }.join('')

      out.close if [File, String].include?(out.class)
    end

  end

  class ProfileSequenceLocalAlignmentLinearGap < ProfileSequenceAlignment

    attr_reader :score, :point, :max_m, :max_n

    def initialize(prf, seq, score, point, max_m, max_n)
      super(prf, seq)
      @score      = score
      @point      = point
      @max_m      = max_m
      @max_n      = max_n
      @raw_score  = @score[max_m][max_n]
      traceback(max_m, max_n)
    end

    def traceback(max_m, max_n)
      pss     = @structural_profile.positions
      aas     = @sequence.amino_acids
      str_cnt = @structural_profile.sequences.size
      ali_pss = []
      ali_aas = []
      m, n    = max_m, max_n
      log_fmt = "%-12s : %-5s : %12s"

      loop do
        case @point[m][n]
        when NONE
          break
        when DIAG
          ali_aas << aas[m-1]
          ali_pss << pss[n-1]
          m -= 1
          n -= 1
          log_mat = "M[#{m}][#{n}]"
          log_prb = "#{ali_pss[-1].probe} <=> #{ali_aas[-1]}"
          $logger.debug log_fmt % [log_mat, "DIAG", log_prb]
        when LEFT
          ali_pss << pss[n-1]
          ali_aas << '-'
          n -= 1
          log_mat = "M[#{m}][#{n}]"
          log_prb = "#{ali_pss[-1].probe} <=> #{ali_aas[-1]}"
          $logger.debug log_fmt % [log_mat, "LEFT", log_prb]
        when UP
          ali_pss << StructuralProfilePosition.new('-'*str_cnt)
          ali_aas << aas[m-1]
          m -= 1
          log_mat = "M[#{m}][#{n}]"
          log_prb = "#{ali_pss[-1].probe} <=> #{ali_aas[-1]}"
          $logger.debug log_fmt % [log_mat, " UP ", log_prb]
        else
          $logger.error "Something wrong with pointer of #{log_mat}"
          exit 1
        end
      end

      @aligned_amino_acids = ali_aas.reverse!
      @aligned_structural_profile_positions = ali_pss.reverse!
    end

    def calculate_reverse_score
      aligner = ProfileSequenceAligner.new(@structural_profile, @sequence.reverse)
      begin
        aligner.local_alignment_linear_gap_cpp.raw_score
      rescue
        aligner.local_alignment_linear_gap_rb.raw_score
      end
    end

    def calculate_z_score(iter=100)
      scores = NArray.int(iter)
      (0...iter).each do |i|
        aligner   = ProfileSequenceAligner.new(@structural_profile, @sequence.shuffle)
        begin
          scores[i] = aligner.local_alignment_linear_gap_cpp.raw_score
        rescue
          scores[i] = aligner.local_alignment_linear_gap_rb.raw_score
        end
      end
      (@raw_score - scores.mean) / scores.stddev
    end

  end

  class ProfileSequenceLocalAlignmentAffineGap < ProfileSequenceAlignment

    attr_reader :mat_score, :mat_point, :mat_jump,
                :del_score, :del_point, :del_jump,
                :ins_score, :ins_point, :ins_jump

    def initialize(prf, seq,
                   mat_score, mat_point, mat_jump,
                   del_score, del_point, del_jump,
                   ins_score, ins_point, ins_jump,
                   max_m, max_n)
      super(prf, seq)
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
      traceback(max_m, max_n)
    end

    def traceback(max_m, max_n)
      pss     = @structural_profile.positions
      aas     = @sequence.amino_acids
      str_cnt = @structural_profile.sequences.size
      ali_pss = []
      ali_aas = []
      pre_mat = nil
      m, n    = max_m, max_n
      log_fmt = "%-12s : %-5s : %12s"
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
            ali_pss << pss[n-1]
            ali_aas << aas[m-1]
            log_mat = "M[#{m}][#{n}]"
            log_prb = "#{ali_pss[-1].probe} <=> #{ali_aas[-1]}"
            log_str = log_fmt % [log_mat, "JUMP", log_prb]
            $logger.debug log_str
          when :del
            ali_pss << pss[n-1]
            ali_aas << '-'
            log_mat = "D[#{m}][#{n}]"
            log_prb = "#{ali_pss[-1].probe} <=> #{ali_aas[-1]}"
            log_str = log_fmt % [log_mat, "JUMP", log_prb]
            $logger.debug log_str
          when :ins
            ali_pss << StructuralProfilePosition.new('-'*str_cnt)
            ali_aas << aas[m-1]
            log_mat = "I[#{m}][#{n}]"
            log_prb = "#{ali_pss[-1].probe} <=> #{ali_aas[-1]}"
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
          ali_aas << aas[m-1]
          ali_pss << pss[n-1]
          m -= 1
          n -= 1
          log_mat = "M[#{m}][#{n}]"
          log_prb = "#{ali_pss[-1].probe} <=> #{ali_aas[-1]}"
          log_str = log_fmt % [log_mat, "DIAG", log_prb]
          $logger.debug log_str
        when LEFT
          ali_aas << '-'
          ali_pss << pss[n-1]
          n -= 1
          log_mat = "D[#{m}][#{n}]"
          log_prb = "#{ali_pss[-1].probe} <=> #{ali_aas[-1]}"
          log_str = log_fmt % [log_mat, "LEFT", log_prb]
          $logger.debug log_str
        when UP
          ali_aas << aas[m-1]
          ali_pss << StructuralProfilePosition.new('-'*str_cnt)
          m -= 1
          log_mat = "I[#{m}][#{n}]"
          log_prb = "#{ali_pss[-1].probe} <=> #{ali_aas[-1]}"
          log_str = log_fmt % [log_mat, " UP ", log_prb]
          $logger.debug log_str
        else
          $logger.error "Something wrong with pointing stage at m: #{m}, n: #{n}"
          exit 1
        end
      end

      @aligned_amino_acids = ali_aas.reverse!
      @aligned_structural_profile_positions = ali_pss.reverse!
    end

    def calculate_reverse_score
      aligner = ProfileSequenceAligner.new(@structural_profile, @sequence.reverse)
      begin
        aligner.local_alignment_affine_gap_cpp.raw_score
      rescue
        aligner.local_alignment_affine_gap_rb.raw_score
      end
    end

    def calculate_z_score(iter=100)
      scores = NArray.int(iter)
      (0...iter).each do |i|
        aligner   = ProfileSequenceAligner.new(@structural_profile, @sequence.shuffle)
        begin
          scores[i] = aligner.local_alignment_affine_gap_cpp.raw_score
        rescue
          scores[i] = aligner.local_alignment_affine_gap_rb.raw_score
        end
      end
      (@raw_score - scores.mean) / scores.stddev
    end

  end

  class ProfileSequenceGlobalAlignmentLinearGap < ProfileSequenceAlignment

    attr_reader :score, :point

    def initialize(prf, seq, score, point)
      super(prf, seq)
      @score      = score
      @point      = point
      @raw_score  = @score[-1][-1]
      begin
        traceback_cpp
      rescue
        traceback
      end
    end

    def traceback
      aas     = @sequence.amino_acids
      pss     = @structural_profile.positions
      str_cnt = @structural_profile.sequences.size
      ali_aas = []
      ali_pss = []
      m, n    = @sequence.length, @structural_profile.length
      log_fmt = "%-12s : %-5s : %12s"

      loop do
        log_mat = "M[#{m}][#{n}]"
        log_prb = "#{ali_pss[-1].probe} <=> #{ali_aas[-1]}"

        case @point[m][n]
        when NONE
          break
        when DIAG
          ali_aas << aas[m-1]
          ali_pss << pss[n-1]
          m -= 1
          n -= 1
          $logger.debug log_fmt % [log_mat, "DIAG", log_prb]
        when LEFT
          ali_aas << '-'
          ali_pss << pss[n-1]
          n -= 1
          $logger.debug log_fmt % [log_mat, "LEFT", log_prb]
        when UP
          ali_aas << aas[m-1]
          ali_pss << StructuralProfilePosition.new('-'*str_cnt)
          m -= 1
          $logger.debug log_fmt % [log_mat, " UP ", log_prb]
        else
          $logger.error "Something wrong with pointer of #{log_mat}"
          exit 1
        end
      end

      @aligned_amino_acids = ali_aas.reverse!
      @aligned_structural_profile_positions = ali_pss.reverse!
    end

    inline(:C) do |builder|
      builder.add_compile_flags '-x c++', '-lstdc++'
      builder.c_raw <<-EOCPP
        static VALUE traceback_cpp(int argc, VALUE *argv, VALUE self) {
          VALUE pnt = rb_iv_get(self, "@point");
          VALUE stp = rb_iv_get(self, "@structural_profile");
          VALUE seq = rb_iv_get(self, "@sequence");
          VALUE pss = rb_funcall(stp, rb_intern("positions"), 0);
          VALUE aas = rb_funcall(seq, rb_intern("amino_acids"), 0);
          VALUE str_cnt = rb_funcall(rb_funcall(stp, rb_intern("sequences"), 0), rb_intern("size"), 0);
          VALUE ali_pss = rb_ary_new();
          VALUE ali_aas = rb_ary_new();

          long pi = NUM2LONG(rb_funcall(stp, rb_intern("length"), 0));
          long si = NUM2LONG(rb_funcall(seq, rb_intern("length"), 0));

          while(1) {
            int p = FIX2INT(rb_ary_entry(rb_ary_entry(pnt, si), pi));
            if (p == 1) {
              // UP
              VALUE gap = rb_str_new2("-");
              for (int i=0; i<NUM2INT(str_cnt); i++) {
                rb_str_concat(gap, rb_str_new2("-"));
              }
              VALUE args[1];
              args[0] = gap;
              rb_ary_push(ali_pss, rb_class_new_instance(1, args, rb_path2class("StructuralProfilePosition")));
              rb_ary_push(ali_aas, rb_ary_entry(aas, si-1));
              si -= 1;
            } else if (p == 2) {
              // LEFT
              rb_ary_push(ali_pss, rb_ary_entry(pss, pi-1));
              rb_ary_push(ali_aas, rb_str_new2("-"));
              pi -= 1;
            } else if (p == 3) {
              // DIAG
              rb_ary_push(ali_pss, rb_ary_entry(pss, pi-1));
              rb_ary_push(ali_aas, rb_ary_entry(aas, si-1));
              pi -= 1;
              si -= 1;
            } else if (p == 0) {
              // NONE
              break;
            } else {
              rb_fatal("Something wrong with point matrix");
            }
          }

          rb_iv_set(self, "@aligned_structural_profile_positions", rb_ary_reverse(ali_pss));
          rb_iv_set(self, "@aligned_amino_acids", rb_ary_reverse(ali_aas));

          return Qnil;
        }
      EOCPP
    end

    def calculate_reverse_score
      aligner = ProfileSequenceAligner.new(@structural_profile, @sequence.reverse)
      begin
        aligner.global_alignment_linear_gap_cpp.raw_score
      rescue
        aligner.global_alignment_linear_gap_rb.raw_score
      end
    end

    def calculate_z_score(iter=100)
      scores = NArray.int(iter)
      (0...iter).each do |i|
        aligner = ProfileSequenceAligner.new(@structural_profile, @sequence.shuffle)
        begin
          scores[i] = aligner.global_alignment_linear_gap_cpp.raw_score
        rescue
          scores[i] = aligner.global_alignment_linear_gap_rb.raw_score
        end
      end
      (@raw_score - scores.mean) / scores.stddev
    end

  end

  class ProfileSequenceGlobalAlignmentAffineGap < ProfileSequenceAlignment

    attr_reader :mat_score, :mat_point, :mat_jump,
                :del_score, :del_point, :del_jump,
                :ins_score, :ins_point, :ins_jump

    def initialize(prf, seq,
                   mat_score, mat_point, mat_jump,
                   del_score, del_point, del_jump,
                   ins_score, ins_point, ins_jump)
      super(prf, seq)
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
      traceback
    end

    def traceback
      pss     = @structural_profile.positions
      aas     = @sequence.amino_acids
      str_cnt = @structural_profile.sequences.size
      ali_pss = []
      ali_aas = []
      pre_mat = nil
      m, n    = aas.length, pss.length
      log_fmt = "%-12s : %-5s : %12s"
      cur_score, cur_point, cur_jump =
        case @raw_score
        when @mat_score[-1][-1]
          pre_mat = :mat
          [@mat_score, @mat_point, @mat_jump]
        when @del_score[-1][-1]
          pre_mat = :del
          [@del_score, @del_point, @del_jump]
        when @ins_score[-1][-1]
          pre_mat = :ins
          [@ins_score, @ins_point, @ins_jump]
        else
          $logger.error "Raw score doesn't match with any matrix"
          exit 1
        end

      $logger.debug "START TRACEBACK"

      loop do
        # jumping between matrices
        until cur_jump[m][n].nil?
          case pre_mat
          when :mat
            ali_pss << pss[n-1]
            ali_aas << aas[m-1]
            log_mat = "M[#{m}][#{n}]"
            log_prb = "#{ali_pss[-1].probe} <=> #{ali_aas[-1]}"
            log_str = log_fmt % [log_mat, "JUMP", log_prb]
            $logger.debug log_str
          when :del
            ali_pss << pss[n-1]
            ali_aas << '-'
            log_mat = "D[#{m}][#{n}]"
            log_prb = "#{ali_pss[-1].probe} <=> #{ali_aas[-1]}"
            log_str = log_fmt % [log_mat, "JUMP", log_prb]
            $logger.debug log_str
          when :ins
            ali_pss << StructuralProfilePosition.new('-'*str_cnt)
            ali_aas << aas[m-1]
            log_mat = "I[#{m}][#{n}]"
            log_prb = "#{ali_pss[-1].probe} <=> #{ali_aas[-1]}"
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
          ali_aas << aas[m-1]
          ali_pss << pss[n-1]
          m -= 1
          n -= 1
          log_mat = "M[#{m}][#{n}]"
          log_prb = "#{ali_pss[-1].probe} <=> #{ali_aas[-1]}"
          log_str = log_fmt % [log_mat, "DIAG", log_prb]
          $logger.debug log_str
        when LEFT
          ali_aas << '-'
          ali_pss << pss[n-1]
          n -= 1
          log_mat = "D[#{m}][#{n}]"
          log_prb = "#{ali_pss[-1].probe} <=> #{ali_aas[-1]}"
          log_str = log_fmt % [log_mat, "LEFT", log_prb]
          $logger.debug log_str
        when UP
          ali_aas << aas[m-1]
          ali_pss << StructuralProfilePosition.new('-'*str_cnt)
          m -= 1
          log_mat = "I[#{m}][#{n}]"
          log_prb = "#{ali_pss[-1].probe} <=> #{ali_aas[-1]}"
          log_str = log_fmt % [log_mat, " UP ", log_prb]
          $logger.debug log_str
        else
          $logger.error "Something wrong with pointing stage at m: #{m}, n: #{n}"
          exit 1
        end
      end

      @aligned_amino_acids = ali_aas.reverse!
      @aligned_structural_profile_positions = ali_pss.reverse!
    end

    def calculate_reverse_score
      aligner = ProfileSequenceAligner.new(@structural_profile, @sequence.reverse)
      begin
        aligner.global_alignment_affine_gap_cpp.raw_score
      rescue
        aligner.global_alignment_affine_gap_rb.raw_score
      end
    end

    def calculate_z_score(iter=100)
      scores = NArray.int(iter)
      (0...iter).each do |i|
        aligner   = ProfileSequenceAligner.new(@structural_profile, @sequence.shuffle)
        begin
          scores[i] = aligner.global_alignment_affine_gap_cpp.raw_score
        rescue
          scores[i] = aligner.global_alignment_affine_gap_rb.raw_score
        end
      end
      (@raw_score - scores.mean) / scores.stddev
    end

  end
end
