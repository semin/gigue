module Gigue
  class ProfileSequenceAligner

    attr_reader :structural_profile, :sequence

    def initialize(prf, seq)
      @structural_profile = prf
      @sequence           = seq
    end

    def local_alignment_with_linear_gap_penalty(gap_del=100, gap_ins=100)
      pss   = @structural_profile.positions
      aas   = @sequence.amino_acids
      plen  = @structural_profile.length
      slen  = @sequence.length
      max_s = 0
      max_i = nil
      max_j = nil
      stmat = NArray.object(plen+1, slen+1)

      # initialize score and trace matrix
      (0..plen).each do |n|
        (0..slen).each do |m|
          stmat[n, m] = { :score => 0, :point => nil }
          if    (n == 0 && m == 0)
            stmat[n, m][:point] = NONE
          elsif (n >  0 && m == 0)
            stmat[n, m][:point] = LEFT
          elsif (n == 0 && m >  0)
            stmat[n, m][:point] = UP
          end
        end
      end

      # fill in score and trace matrix
      (1..plen).each do |n|
        (1..slen).each do |m|
          mat = stmat[n-1,  m-1][:score] + pss[n-1].mat_score(aas[m-1])
          del = stmat[n-1,    m][:score] - gap_del
          ins = stmat[n,    m-1][:score] - gap_ins
          max = [0, mat, del, ins].max

          stmat[n, m][:score] = max
          stmat[n, m][:point] = case max
                                when 0    then NONE
                                when mat  then DIAG
                                when del  then LEFT
                                when ins  then UP
                                end
          if max >= max_s
            max_i = n
            max_j = m
            max_s = max
          end
        end
      end

      ProfileSequenceLocalAlignmentLinearGap.new(@structural_profile, @sequence,
                                                 stmat, max_i, max_j)
    end

    def local_alignment_with_affine_gap_penalty
      pss   = @structural_profile.positions
      aas   = @sequence.amino_acids
      plen  = @structural_profile.length
      slen  = @sequence.length
      max_s = 0
      max_i = nil
      max_j = nil
      max_m = nil

      # create score, deletion, and insertion matrices
      mat = NArray.object(plen+1, slen+1)
      del = NArray.object(plen+1, slen+1)
      ins = NArray.object(plen+1, slen+1)

      # initialize score matrix and fill in the first row and column
      prev_gap_ins_ext = nil
      (0..plen).each do |pi|
        (0..slen).each do |si|
          mat[pi, si] = { :score => 0, :point => nil, :jump => nil }
          if    (pi == 0 && si == 0)
            mat[pi, si][:point] = NONE
          elsif (pi == 1 && si == 0)
            mat[pi, si][:score] = -pss[pi-1].gap_del_open
            mat[pi, si][:point] = LEFT
          elsif (pi >  1 && si == 0)
            mat[pi, si][:score] = mat[pi-1, si][:score] - pss[pi-1].gap_del_ext
            mat[pi, si][:point] = LEFT
          elsif (pi == 0 && si == 1)
            mat[pi, si][:score] = -pss[pi].gap_ins_open
            mat[pi, si][:point] = UP
            prev_gap_ins_ext    = pss[pi].gap_ins_ext
          elsif (pi == 0 && si >  1)
            mat[pi, si][:score] = mat[pi, si-1][:score] - prev_gap_ins_ext
            mat[pi, si][:point] = UP
          end
        end
      end

      # initialize deletion and insertion matrices
      infinity = 1/0.0
      (0..plen).each do |pi|
        (0..slen).each do |si|
          del[pi, si] = { :score => (pi == 0 || si == 0 ? -infinity : 0), :point => nil, :jump => nil }
          ins[pi, si] = { :score => (pi == 0 || si == 0 ? -infinity : 0), :point => nil, :jump => nil }
        end
      end

      # fill in match, deletion, and insertion matrices
      prev_gap_ins_ext = pss[0].gap_ins_ext
      (1..plen).each do |n|
        (1..slen).each do |m|
          pi, ai  = n-1, m-1

          # update match matrix
          mat_mat = mat[n-1, m-1][:score] + pss[pi].mat_score(aas[ai])
          mat_del = del[n-1, m-1][:score] + pss[pi].mat_score(aas[ai])
          mat_ins = ins[n-1, m-1][:score] + pss[pi].mat_score(aas[ai])
          mat_max = [0, mat_mat, mat_del, mat_ins].max
          mat[n, m][:score] = mat_max

          case mat_max
          when 0
            mat[n, m][:point] = NONE
            $logger.debug "%-10s %5s" % ["M-#{n}-#{m}", "NONE"]
          when mat_mat
            mat[n, m][:point] = DIAG
            $logger.debug "%-10s %5s %10s" % ["M-#{n}-#{m}", "POINT", "M-#{n-1}-#{m-1}"]
          when mat_del
            jmp = "D-#{n-1}-#{m-1}"
            mat[n, m][:jump] = jmp
            $logger.debug "%-10s %5s %10s" % ["M-#{n}-#{m}", "JUMP", jmp]
          when mat_ins
            jmp = "I-#{n-1}-#{m-1}"
            mat[n, m][:jump] = jmp
            $logger.debug "%-10s %5s %10s" % ["M-#{n}-#{m}", "JUMP", jmp]
          end

          # update deletion matrix
          del_mat = mat[n-1, m][:score] - pss[pi].gap_del_open
          del_del = del[n-1, m][:score] - pss[pi].gap_del_ext
          del_max = [0, del_mat, del_del].max
          del[n, m][:score] = del_max

          case del_max
          when 0
            del[n, m][:point] = NONE
            $logger.debug "%-10s %5s" % ["D-#{n}-#{m}", "NONE"]
          when del_mat
            jmp = "M-#{n-1}-#{m}"
            del[n, m][:jump] = jmp
            $logger.debug "%-10s %5s %10s" % ["D-#{n}-#{m}", "JUMP", jmp]
          when del_del
            del[n, m][:point] = LEFT
            $logger.debug "%-10s %5s %10s" % ["D-#{n}-#{m}", "POINT", "D-#{n-1}-#{m}"]
          end

          # update insertion matrix
          ins_mat = mat[n, m-1][:score] - pss[pi].gap_ins_open
          ins_ins = ins[n, m-1][:score] - prev_gap_ins_ext
          ins_max = [0, ins_mat, ins_ins].max
          ins[n, m][:score] = ins_max

          case ins_max
          when 0
            ins[n, m][:point] = NONE
            $logger.debug "%-10s %5s" % ["I-#{n}-#{m}", "NONE"]
          when ins_mat
            jmp = "M-#{n}-#{m-1}"
            ins[n, m][:jump] = jmp
            prev_gap_ins_ext = pss[pi].gap_ins_ext
            $logger.debug "%-10s %5s %10s" % ["I-#{n}-#{m}", "JUMP", jmp]
          when ins_ins
            ins[n, m][:point] = UP
            $logger.debug "%-10s %5s %10s" % ["I-#{n}-#{m}", "POINT", "D-#{n}-#{m-1}"]
          end

          # keep the record of a matrix and its indexes having maximum value so far
          max = [mat_max, del_max, ins_max].max
          if max >= max_s
            max_i = n
            max_j = m
            max_s = max
            max_m = case max
                    when mat_max then :mat
                    when del_max then :del
                    when ins_max then :ins
                    end
          end
        end
      end

      ProfileSequenceLocalAlignmentAffineGap.new(@structural_profile, @sequence,
                                                 mat, del, ins, max_i, max_j, max_m)
    end

    def global_alignment_with_linear_gap_penalty(gap_del=100, gap_ins=100)
      pss   = @structural_profile.positions
      aas   = @sequence.amino_acids
      plen  = @structural_profile.length
      slen  = @sequence.length

      # 1. create score and point matrices
      score = Array.new(slen+1)
      point = Array.new(slen+1)

      # 2. initialize score and point matrices
      (0..slen).each do |si|
        score[si] = Array.new(plen+1)
        point[si] = Array.new(plen+1)
        (0..plen).each do |pi|
          if    (pi == 0 && si == 0)
            score[si][pi] = 0
            point[si][pi] = NONE
          elsif (pi >  0 && si == 0)
            score[si][pi] = -gap_del*pi
            point[si][pi] = LEFT
          elsif (pi == 0 && si >  0)
            score[si][pi] = -gap_ins*si
            point[si][pi] = UP
          end
        end
      end

      # 3. fill in score and point matrices
      (1..slen).each do |m|
        (1..plen).each do |n|
          mat = score[m-1][n-1] + pss[n-1].mat_score(aas[m-1])
          del = score[m][n-1] - gap_del
          ins = score[m-1][n] - gap_ins
          max = [mat, del, ins].max

          score[m][n] = max
          point[m][n] = case max
                        when mat then DIAG
                        when del then LEFT
                        when ins then UP
                        end
        end
      end

      ProfileSequenceGlobalAlignmentLinearGap.new(@structural_profile, @sequence,
                                                  score, point)
    end

    inline(:C) do |builder|
      builder.add_compile_flags '-x c++', '-lstdc++'
      builder.c_raw %q{
        static VALUE global_alignment_with_linear_gap_penalty_cpp(int argc, VALUE *argv, VALUE self) {
          long gdel   = argv[0] == Qnil ? 100 : NUM2LONG(argv[0]);
          long gins   = argv[1] == Qnil ? 100 : NUM2LONG(argv[1]);
          VALUE stp   = rb_iv_get(self, "@structural_profile");
          VALUE seq   = rb_iv_get(self, "@sequence");
          VALUE pss   = rb_funcall(stp, rb_intern("positions"), 0);
          VALUE aas   = rb_funcall(seq, rb_intern("amino_acids"), 0);
          long plen   = NUM2LONG(rb_funcall(stp, rb_intern("length"), 0));
          long slen   = NUM2LONG(rb_funcall(seq, rb_intern("length"), 0));
          VALUE score = rb_ary_new2(slen+1);
          VALUE point = rb_ary_new2(slen+1);

          for (long si=0; si < slen+1; si++) {
            VALUE score_row = rb_ary_new2(plen+1);
            VALUE point_row = rb_ary_new2(plen+1);

            for (long pi=0; pi < plen+1; pi++) {
              if ((pi == 0) && (si == 0)) {
                rb_ary_store(score_row, pi, LONG2NUM(0));
                rb_ary_store(point_row, pi, INT2FIX(0));
              } else if (pi >  0 && si == 0) {
                rb_ary_store(score_row, pi, LONG2NUM(-gdel*pi));
                rb_ary_store(point_row, pi, INT2FIX(2));
              } else if (pi == 0 && si >  0) {
                rb_ary_store(score_row, pi, LONG2NUM(-gins*pi));
                rb_ary_store(point_row, pi, INT2FIX(1));
              }
            }
            rb_ary_store(score, si, score_row);
            rb_ary_store(point, si, point_row);
          }

          for (long m=1; m < slen+1; m++) {
            for (long n=1; n < plen+1; n++) {
              long mat = NUM2LONG(rb_ary_entry(rb_ary_entry(score, m-1), n-1)) +
                          NUM2LONG(rb_funcall(rb_ary_entry(pss, n-1), rb_intern("mat_score"), 1, rb_ary_entry(aas, m-1)));
              long del = NUM2LONG(rb_ary_entry(rb_ary_entry(score, m), n-1)) - gdel;
              long ins = NUM2LONG(rb_ary_entry(rb_ary_entry(score, m-1), n)) - gins;

              if (mat >= ins) {
                if (mat >= del) {
                  rb_ary_store(rb_ary_entry(score, m), n, LONG2NUM(mat));
                  rb_ary_store(rb_ary_entry(point, m), n, INT2FIX(3));
                } else {
                  rb_ary_store(rb_ary_entry(score, m), n, LONG2NUM(del));
                  rb_ary_store(rb_ary_entry(point, m), n, INT2FIX(2));
                }
              } else {
                if (ins > del) {
                  rb_ary_store(rb_ary_entry(score, m), n, LONG2NUM(ins));
                  rb_ary_store(rb_ary_entry(point, m), n, INT2FIX(1));
                } else {
                  rb_ary_store(rb_ary_entry(score, m), n, LONG2NUM(del));
                  rb_ary_store(rb_ary_entry(point, m), n, INT2FIX(2));
                }
              }
            }
          }
          VALUE args[4];
          args[0] = stp;
          args[1] = seq;
          args[2] = score;
          args[3] = point;

          VALUE psgl = rb_class_new_instance(4, args, rb_path2class("ProfileSequenceGlobalAlignmentLinearGap"));
          return psgl;
        }
      }
    end

    def global_alignment_with_affine_gap_penalty
      pss   = @structural_profile.positions
      aas   = @sequence.amino_acids
      plen  = @structural_profile.length
      slen  = @sequence.length

      # 1. create score, point, jump matrices for match, deletion, and insertion
      mat_score = Array.new(slen+1)
      mat_point = Array.new(slen+1)
      mat_jump  = Array.new(slen+1)

      del_score = Array.new(slen+1)
      del_point = Array.new(slen+1)
      del_jump  = Array.new(slen+1)

      ins_score = Array.new(slen+1)
      ins_point = Array.new(slen+1)
      ins_jump  = Array.new(slen+1)

      # 2. initialize match matrices
      prev_gap_ins_ext = nil

      (0..slen).each do |si|
        mat_score[si] = Array.new(plen+1)
        mat_point[si] = Array.new(plen+1)
        mat_jump[si]  = Array.new(plen+1)

        (0..plen).each do |pi|
          if    (pi == 0 && si == 0)
            mat_score[si][pi] = 0
            mat_point[si][pi] = NONE
          elsif (pi == 1 && si == 0)
            mat_score[si][pi] = -pss[pi-1].gap_del_open
            mat_point[si][pi] = LEFT
          elsif (pi >  1 && si == 0)
            mat_score[si][pi] = mat_score[si][pi-1] - pss[pi-1].gap_del_ext
            mat_point[si][pi] = LEFT
          elsif (pi == 0 && si == 1)
            mat_score[si][pi] = -pss[pi].gap_ins_open
            mat_point[si][pi] = UP
            prev_gap_ins_ext  = pss[pi].gap_ins_ext
          elsif (pi == 0 && si >  1)
            mat_score[si][pi] = mat_score[si-1][pi] - prev_gap_ins_ext
            mat_point[si][pi] = UP
          end
        end
      end

      # 3. initialize deletion and insertion matrices
      infinity = 1/0.0

      (0..slen).each do |si|
        del_score[si] = Array.new(plen+1)
        del_point[si] = Array.new(plen+1)
        del_jump[si]  = Array.new(plen+1)

        ins_score[si] = Array.new(plen+1)
        ins_point[si] = Array.new(plen+1)
        ins_jump[si]  = Array.new(plen+1)

        (0..plen).each do |pi|
          del_score[si][pi] = -infinity if (pi == 0 || si == 0)
          ins_score[si][pi] = -infinity if (pi == 0 || si == 0)
        end
      end

      # 4. fill in match, deletion, and insertion matrices
      prev_gap_ins_ext = pss[0].gap_ins_ext

      (1..slen).each do |m|
        (1..plen).each do |n|
          pi, ai  = n-1, m-1

          mat_mat = mat_score[m-1][n-1] + pss[pi].mat_score(aas[ai])
          mat_del = del_score[m-1][n-1] + pss[pi].mat_score(aas[ai])
          mat_ins = ins_score[m-1][n-1] + pss[pi].mat_score(aas[ai])
          mat_max = [mat_mat, mat_del, mat_ins].max
          mat_score[m][n] = mat_max

          case mat_max
          when mat_mat
            mat_point[m][n] = DIAG
            $logger.debug "%-10s %5s %10s" % ["M[#{m}][#{n}]", "POINT", "M[#{m-1}][#{n-1}]"]
          when mat_del
            mat_jump[m][n] = ['D', m-1, n-1]
            $logger.debug "%-10s %5s %10s" % ["M[#{m}][#{n}]", "JUMP", "D[#{m-1}][#{n-1}]"]
          when mat_ins
            mat_jump[m][n] = ['I', m-1, n-1]
            $logger.debug "%-10s %5s %10s" % ["M[#{m}][#{n}]", "JUMP", "I[#{m-1}][#{n-1}]"]
          end

          del_mat = mat_score[m][n-1] - pss[pi].gap_del_open
          del_del = del_score[m][n-1] - pss[pi].gap_del_ext
          del_max = [del_mat, del_del].max
          del_score[m][n] = del_max

          case del_max
          when del_mat
            del_jump[m][n] = ['M', m, n-1]
            $logger.debug "%-10s %5s %10s" % ["D[#{m}][#{n}]", "JUMP", "M[#{m}][#{n-1}]"]
          when del_del
            del_point[m][n] = LEFT
            $logger.debug "%-10s %5s %10s" % ["D[#{m}][#{n}]", "POINT", "D[#{m}][#{n-1}]"]
          end

          ins_mat = mat_score[m-1][n] - pss[pi].gap_ins_open
          ins_ins = ins_score[m-1][n] - prev_gap_ins_ext
          ins_max = [ins_mat, ins_ins].max
          ins_score[m][n] = ins_max

          case ins_max
          when ins_mat
            ins_jump[m][n] = ['M', m-1, n]
            prev_gap_ins_ext = pss[pi].gap_ins_ext
            $logger.debug "%-10s %5s %10s" % ["I[#{m}][#{n}]", "JUMP", "M[#{m-1}][#{n}]"]
          when ins_ins
            ins_point[m][n] = UP
            $logger.debug "%-10s %5s %10s" % ["I[#{m}][#{n}]", "POINT", "I[#{m-1}][#{n}]"]
          end
        end
      end

      ProfileSequenceGlobalAlignmentAffineGap.new(@structural_profile, @sequence,
                                                  mat_score, mat_point, mat_jump,
                                                  del_score, del_point, del_jump,
                                                  ins_score, ins_point, ins_jump)
    end

  end
end
