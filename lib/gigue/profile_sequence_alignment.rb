module Gigue
  class ProfileSequenceAlignment

    attr_reader :structural_profile, :sequence, :raw_score,
                :aligned_structural_profile_positions,
                :aligned_amino_acids

    def initialize(prf, seq)
      @structural_profile = prf
      @sequence           = seq
    end

    # memoize function for Z-score
    def z_score
      @z_score ||= calculate_z_score
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

      # print aligned query sequence
      out.puts ">#{@sequence.code}"
      out.puts "sequence" if opts[:type] == :pir
      out.puts @aligned_amino_acids.map_with_index { |a, ai|
        a + (ai > 0 && (ai+1) % opts[:width] == 0 ? "\n" : '')
      }.join('')

      out.close if [File, String].include?(out.class)
    end

  end

  class ProfileSequenceLocalAlignmentLinearGap < ProfileSequenceAlignment

    attr_reader :stmatrix

    def initialize(prf, seq, stmatrix, i, j)
      super(prf, seq)
      @stmatrix   = stmatrix
      @raw_score  = @stmatrix[i, j][:score]
      @aligned_structural_profile_positions, @aligned_amino_acids = traceback(i, j)
    end

    def traceback(i, j)
      pss     = @structural_profile.positions
      aas     = @sequence.amino_acids
      str_cnt = @structural_profile.sequences.size
      ali_prf = []
      ali_seq = []

      loop do
        case @stmatrix[i, j][:point]
        when NONE
          break
        when DIAG
          ali_prf << pss[i-1]
          ali_seq << aas[j-1]
          i -= 1
          j -= 1
          $logger.debug "%-15s : %-4s : %s" % ["MAT[#{i}, #{j}]", "DIAG", "#{ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
        when LEFT
          ali_prf << pss[i-1]
          ali_seq << '-'
          i -= 1
          $logger.debug "%-15s : %-4s : %s" % ["MAT[#{i}, #{j}]", "INS", "#{ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
        when UP
          ali_prf << StructuralProfilePosition.new('-'*str_cnt)
          ali_seq << aas[j-1]
          j -= 1
          $logger.debug "%-15s : %-4s : %s" % ["MAT[#{i}, #{j}]", "DEL", "#{ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
        else
          $logger.error "Something wrong with pointing stage at i: #{i}, j: #{j}"
          exit 1
        end
      end
      [ali_prf.reverse!, ali_seq.reverse!]
    end

    def calculate_z_score(iter=100)
      scores = NArray.int(100)
      (0...iter).each do |i|
        aligner   = ProfileSequenceAligner.new(@structural_profile, @sequence.shuffle)
        scores[i] = aligner.local_alignment_with_linear_gap_penalty.raw_score
      end
      (@raw_score - scores.mean) / scores.stddev
    end

  end

  class ProfileSequenceLocalAlignmentAffineGap < ProfileSequenceAlignment

    attr_reader :match_stmatrix, :insertion_stmatrix, :deletion_stmatrix

    def initialize(prf, seq,
                   match_stmatrix, deletion_stmatrix, insertion_stmatrix,
                   i, j, max_mat)
      super(prf, seq)
      @match_stmatrix     = match_stmatrix
      @deletion_stmatrix  = deletion_stmatrix
      @insertion_stmatrix = insertion_stmatrix
      @raw_score          = case max_mat
                            when :mat then @match_stmatrix[i,j][:score]
                            when :del then @deletion_stmatrix[i,j][:score]
                            when :ins then @insertion_stmatrix[i,j][:score]
                            end
      @aligned_structural_profile_positions, @aligned_amino_acids = traceback(i, j, max_mat)
    end

    def traceback(i, j, max_mat)
      pss     = @structural_profile.positions
      aas     = @sequence.amino_acids
      str_cnt = @structural_profile.sequences.size
      ali_prf = []
      ali_seq = []
      pre_mat = nil
      cur_mat = case max_mat
                when :mat
                  $logger.debug "START TRACEBACK FROM MAT MATRIX"
                  pre_mat = :mat
                  @match_stmatrix
                when :del
                  $logger.debug "START TRACEBACK FROM DEL MATRIX"
                  pre_mat = :del
                  @deletion_stmatrix
                when :ins
                  $logger.debug "START TRACEBACK FROM INS MATRIX"
                  pre_mat = :ins
                  @insertion_stmatrix
                else
                  $logger.error "Maximum matrix tag doesn't match with any matrix"
                  exit 1
                end

      loop do
        until cur_mat[i, j][:jump].nil?
          case pre_mat
          when :mat
            ali_prf << pss[i-1]
            ali_seq << aas[j-1]
          when :del
            ali_prf << pss[i-1]
            ali_seq << '-'
          when :ins
            ali_prf << StructuralProfilePosition.new('-'*str_cnt)
            ali_seq << aas[j-1]
          end

          jm, mi, jj  = cur_mat[i, j][:jump].split('-')
          i, j        = Integer(mi), Integer(jj)
          cur_mat     = case jm
                        when 'M'
                          $logger.debug "%-15s : %-4s : %s" % ["MAT[#{i}, #{j}]", "JUMP", "#{ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
                          pre_mat = :mat
                          @match_stmatrix
                        when 'D'
                          $logger.debug "%-15s : %-4s : %s" % ["DEL[#{i}, #{j}]", "JUMP", "#{ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
                          pre_mat = :del
                          @deletion_stmatrix
                        when 'I'
                          $logger.debug "%-15s : %-4s : %s" % ["INS[#{i}, #{j}]", "JUMP", "#{ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
                          pre_mat = :ins
                          @insertion_stmatrix
                        else
                          $logger.warn "Something wrong in jumping step"
                          exit 1
                        end
        end

        point = cur_mat[i, j][:point]
        case point
        when NONE
          $logger.debug "FINISH TRACEBACK"
          break
        when DIAG
          ali_prf << pss[i-1]
          ali_seq << aas[j-1]
          i -= 1
          j -= 1
          $logger.debug "%-15s : %-4s : %s" % ["MAT[#{i}, #{j}]", "DIAG", "#{ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
        when LEFT
          ali_prf << pss[i-1]
          ali_seq << '-'
          i -= 1
          $logger.debug "%-15s : %-4s : %s" % ["DEL[#{i}, #{j}]", "LEFT", "#{ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
        when UP
          ali_prf << StructuralProfilePosition.new('-'*str_cnt)
          ali_seq << aas[j-1]
          j -= 1
          $logger.debug "%-15s : %-4s : %s" % ["INS[#{i}, #{j}]", "UP", "#{ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
        else
          $logger.error "Something wrong with pointing stage at i: #{i}, j: #{j}"
          exit 1
        end
      end
      [ali_prf.reverse!, ali_seq.reverse!]
    end

    def calculate_z_score(iter=100)
      scores = NArray.int(100)
      (0...iter).each do |i|
        aligner   = ProfileSequenceAligner.new(@structural_profile, @sequence.shuffle)
        scores[i] = aligner.local_alignment_with_affine_gap_penalty.raw_score
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
      pss     = @structural_profile.positions
      aas     = @sequence.amino_acids
      str_cnt = @structural_profile.sequences.size
      ali_prf = []
      ali_seq = []
      i = @structural_profile.length
      j = @sequence.length

      loop do
        case @point[j][i]
        when NONE
          break
        when DIAG
          ali_prf << pss[i-1]
          ali_seq << aas[j-1]
          i -= 1
          j -= 1
          $logger.debug "%-15s : %-4s : %s" % ["MAT[#{i}, #{j}]", "DIAG", "#{ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
        when LEFT
          ali_prf << pss[i-1]
          ali_seq << '-'
          i -= 1
          $logger.debug "%-15s : %-4s : %s" % ["MAT[#{i}, #{j}]", "INS", "#{ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
        when UP
          ali_prf << StructuralProfilePosition.new('-'*str_cnt)
          ali_seq << aas[j-1]
          j -= 1
          $logger.debug "%-15s : %-4s : %s" % ["MAT[#{i}, #{j}]", "DEL", "#{ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
        else
          $logger.error "Something wrong with pointing stage at i: #{i}, j: #{j}"
          exit 1
        end
      end

      @aligned_structural_profile_positions, @aligned_amino_acids = ali_prf.reverse!, ali_seq.reverse!
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
          VALUE ali_prf = rb_ary_new();
          VALUE ali_seq = rb_ary_new();

          long i = NUM2LONG(rb_funcall(stp, rb_intern("length"), 0));
          long j = NUM2LONG(rb_funcall(seq, rb_intern("length"), 0));

          while(1) {
            int p = FIX2INT(rb_ary_entry(rb_ary_entry(pnt, j), i));
            if (p == 1) {
              // UP
              VALUE gap = rb_str_new2("-");
              for (int i=0; i<NUM2INT(str_cnt); i++) {
                rb_str_concat(gap, rb_str_new2("-"));
              }
              VALUE args[1];
              args[0] = gap;
              rb_ary_push(ali_prf, rb_class_new_instance(1, args, rb_path2class("StructuralProfilePosition")));
              rb_ary_push(ali_seq, rb_ary_entry(aas, j-1));
              j -= 1;
            } else if (p == 2) {
              // LEFT
              rb_ary_push(ali_prf, rb_ary_entry(pss, i-1));
              rb_ary_push(ali_seq, rb_str_new2("-"));
              i -= 1;
            } else if (p == 3) {
              // DIAG
              rb_ary_push(ali_prf, rb_ary_entry(pss, i-1));
              rb_ary_push(ali_seq, rb_ary_entry(aas, j-1));
              i -= 1;
              j -= 1;
            } else if (p == 0) {
              // NONE
              break;
            } else {
              rb_fatal("Something wrong with point matrix");
            }
          }

          rb_iv_set(self, "@aligned_structural_profile_positions", rb_ary_reverse(ali_prf));
          rb_iv_set(self, "@aligned_amino_acids", rb_ary_reverse(ali_seq));

          return Qnil;
        }
      EOCPP
    end

    def calculate_z_score(iter=100)
      scores = NArray.int(100)
      (0...iter).each do |i|
        aligner   = ProfileSequenceAligner.new(@structural_profile, @sequence.shuffle)
        scores[i] = aligner.global_alignment_with_linear_gap_penalty.raw_score
      end
      (@raw_score - scores.mean) / scores.stddev
    end

  end

  class ProfileSequenceGlobalAlignmentAffineGap < ProfileSequenceAlignment

    attr_reader :match_score, :match_point, :match_jump,
                :deletion_score, :deletion_point, :deletion_jump,
                :insertion_score, :insertion_point, :insertion_jump

    def initialize(prf, seq,
                   mat_score, mat_point, mat_jump,
                   del_score, del_point, del_jump,
                   ins_score, ins_point, ins_jump)
      super(prf, seq)
      @match_score      = mat_score
      @match_point      = mat_point
      @match_jump       = mat_jump
      @deletion_score   = del_score
      @deletion_point   = del_point
      @deletion_jump    = del_jump
      @insertion_score  = ins_score
      @insertion_point  = ins_point
      @insertion_jump   = ins_jump
      @raw_score        = [
        @match_score[-1][-1],
        @deletion_score[-1][-1],
        @insertion_score[-1][-1]
      ].max
      traceback
    end

    def traceback
      pss     = @structural_profile.positions
      aas     = @sequence.amino_acids
      str_cnt = @structural_profile.sequences.size
      ali_prf = []
      ali_seq = []
      pre_mat = nil
      pi, si  = [pss.length, aas.length]
      cur_score, cur_point, cur_jump =
        case @raw_score
        when @match_score[-1][-1]
          pre_mat = :mat
          [@match_score, @match_point, @match_jump]
        when @deletion_score[-1][-1]
          pre_mat = :del
          [@deletion_score, @deletion_point, @deletion_jump]
        when @insertion_score[-1][-1]
          pre_mat = :ins
          [@insertion_score, @insertion_point, @insertion_jump]
        else
          $logger.error "Raw score doesn't match with any matrix"
          exit 1
        end

      $logger.debug "START TRACEBACK"

      loop do
        # jumping between matrices
        until cur_jump[si][pi].nil?
          case pre_mat
          when :mat
            ali_prf << pss[pi-1]
            ali_seq << aas[si-1]
          when :del
            ali_prf << pss[pi-1]
            ali_seq << '-'
          when :ins
            ali_prf << StructuralProfilePosition.new('-'*str_cnt)
            ali_seq << aas[si-1]
          end

          mat, si, pi = cur_jump[si][pi]
          cur_score, cur_point, cur_jump =
            case mat
            when 'M'
              $logger.debug "%-15s : %-4s : %s" %
              ["MAT[#{si}][#{pi}]", "JUMP", "#{ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
              pre_mat = :mat
              [@match_score, @match_point, @match_jump]
            when 'D'
              $logger.debug "%-15s : %-4s : %s" %
              ["DEL[#{si}][#{pi}]", "JUMP", "#{ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
              pre_mat = :del
              [@deletion_score, @deletion_point, @deletion_jump]
            when 'I'
              $logger.debug "%-15s : %-4s : %s" %
              ["INS[#{si}][#{pi}]", "JUMP", "#{ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
              pre_mat = :ins
              [@insertion_score, @insertion_point, @insertion_jump]
            else
              $logger.warn "Something wrong in jumping step"
              exit 1
            end
        end

        # following pointers in a matrix
        point = cur_point[si][pi]

        case point
        when NONE
          $logger.debug "FINISH TRACEBACK"
          break
        when DIAG
          ali_prf << pss[pi-1]
          ali_seq << aas[si-1]
          pi -= 1
          si -= 1
          $logger.debug "%-15s : %-4s : %s" %
          ["MAT[#{si}][#{pi}]", "DIAG", "#{ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
        when LEFT
          ali_prf << pss[pi-1]
          ali_seq << '-'
          pi -= 1
          $logger.debug "%-15s : %-4s : %s" %
          ["DEL[#{si}][#{pi}]", "LEFT", "#{ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
        when UP
          ali_prf << StructuralProfilePosition.new('-'*str_cnt)
          ali_seq << aas[si-1]
          si -= 1
          $logger.debug "%-15s : %-4s : %s" %
          ["INS[#{si}][#{pi}]", "UP", "#{ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
        else
          $logger.error "Something wrong with pointing stage at si: #{si}, pi: #{pi}"
          exit 1
        end
      end

      @aligned_structural_profile_positions, @aligned_amino_acids = ali_prf.reverse!, ali_seq.reverse!
    end

    def calculate_z_score(iter=100)
      scores = NArray.int(100)
      (0...iter).each do |i|
        aligner   = ProfileSequenceAligner.new(@structural_profile, @sequence.shuffle)
        scores[i] = aligner.global_alignment_with_affine_gap_penalty.raw_score
      end
      (@raw_score - scores.mean) / scores.stddev
    end

  end
end
