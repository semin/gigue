module Gigue
  class ProfileSequenceAligner

    attr_reader :structural_profile, :sequence, :algorithm

    def initialize(prf, seq, algo=nil)
      @structural_profile = prf
      @sequence           = seq
      @str_seq_ratio      = prf.length / Float(seq.length)
    end

    def local_alignment_linear_gap(options={})
      @algorithm = :local
      begin
        local_alignment_linear_gap_cpp(options)
      rescue
        local_alignment_linear_gap_rb(options)
      end
    end

    def local_alignment_affine_gap
      @algorithm = :local
      begin
        local_alignment_affine_gap_cpp
      rescue
        local_alignment_affine_gap_rb
      end
    end

    def global_alignment_linear_gap
      if    (@str_seq_ratio > 1.5)     then @algorithm = :glolocstr
      elsif (@str_seq_ratio < 0.6667)  then @algorithm = :glolocseq
      else @algorithm = :global
      end
      begin
        global_alignment_linear_gap_cpp(options={})
      rescue
        global_alignment_linear_gap_rb(options={})
      end
    end

    def global_alignment_affine_gap
      if    (@str_seq_ratio > 1.5)     then @algorithm = :glolocstr
      elsif (@str_seq_ratio < 0.6667)  then @algorithm = :glolocseq
      else @algorithm = :global
      end
      begin
        global_alignment_affine_gap_cpp
      rescue
        global_alignment_affine_gap_rb
      end
    end

    inline(:C) do |builder|
      builder.add_compile_flags '-x c++', '-lstdc++'
      builder.c_raw <<-EOCPP
        static VALUE local_alignment_linear_gap_cpp(int argc, VALUE *argv, VALUE self) {
          VALUE options = argv[0] == Qnil ? rb_hash_new() : argv[0];
          VALUE opts    = rb_hash_new();

          rb_hash_aset(opts, ID2SYM(rb_intern("gap_del")), INT2NUM(100));
          rb_hash_aset(opts, ID2SYM(rb_intern("gap_ins")), INT2NUM(100));
          rb_funcall(opts, rb_intern("merge!"), 1, options);

          long gdel   = NUM2LONG(rb_hash_aref(opts, ID2SYM(rb_intern("gap_del"))));
          long gins   = NUM2LONG(rb_hash_aref(opts, ID2SYM(rb_intern("gap_ins"))));
          VALUE stp   = rb_iv_get(self, "@structural_profile");
          VALUE seq   = rb_iv_get(self, "@sequence");
          VALUE pss   = rb_funcall(stp, rb_intern("positions"), 0);
          VALUE aas   = rb_funcall(seq, rb_intern("amino_acids"), 0);
          long plen   = NUM2LONG(rb_funcall(stp, rb_intern("length"), 0));
          long slen   = NUM2LONG(rb_funcall(seq, rb_intern("length"), 0));
          long max_s  = 0;
          long max_m  = 0;
          long max_n  = 0;

          // 1. create score and point matrices
          VALUE score = rb_ary_new2(slen+1);
          VALUE point = rb_ary_new2(slen+1);

          // 2. initialize score and point matrices
          for (long m = 0; m < slen+1; m++) {
            VALUE score_row = rb_ary_new2(plen+1);
            VALUE point_row = rb_ary_new2(plen+1);

            rb_ary_store(score, m, score_row);
            rb_ary_store(point, m, point_row);

            for (long n = 0; n < plen+1; n++) {
              if ((n == 0) && (m == 0)) {
                // NONE
                rb_ary_store(score_row, n, INT2NUM(0));
                rb_ary_store(point_row, n, INT2FIX(0));
              } else if ((n > 0) && (m == 0)) {
                // LEFT
                rb_ary_store(score_row, n, LONG2NUM(-gdel*n));
                rb_ary_store(point_row, n, INT2FIX(2));
              } else if ((n == 0) && (m > 0)) {
                // UP
                rb_ary_store(score_row, n, LONG2NUM(-gins*m));
                rb_ary_store(point_row, n, INT2FIX(1));
              }
            }
          }

          // 3. fill in score and point matrices
          for (long m = 1; m < slen+1; m++) {
            for (long n = 1; n < plen+1; n++) {
              long mat = NUM2LONG(rb_ary_entry(rb_ary_entry(score, m-1), n-1)) + NUM2LONG(rb_funcall(rb_ary_entry(pss, n-1), rb_intern("mat_score"), 1, rb_ary_entry(aas, m-1)));
              long del = NUM2LONG(rb_ary_entry(rb_ary_entry(score, m), n-1)) - gdel;
              long ins = NUM2LONG(rb_ary_entry(rb_ary_entry(score, m-1), n)) - gins;
              long max = 0;

              if (mat >= ins) {
                if (mat >= del) {
                  if (mat > 0) {
                    max = mat;
                    rb_ary_store(rb_ary_entry(score, m), n, LONG2NUM(max));
                    rb_ary_store(rb_ary_entry(point, m), n, INT2FIX(3));
                  } else {
                    rb_ary_store(rb_ary_entry(score, m), n, LONG2NUM(max));
                    rb_ary_store(rb_ary_entry(point, m), n, INT2FIX(0));
                  }
                } else {
                  if (del > 0) {
                    max = del;
                    rb_ary_store(rb_ary_entry(score, m), n, LONG2NUM(max));
                    rb_ary_store(rb_ary_entry(point, m), n, INT2FIX(2));
                  } else {
                    rb_ary_store(rb_ary_entry(score, m), n, LONG2NUM(max));
                    rb_ary_store(rb_ary_entry(point, m), n, INT2FIX(0));
                  }
                }
              } else {
                if (ins > del) {
                  if (ins > 0) {
                    max = ins;
                    rb_ary_store(rb_ary_entry(score, m), n, LONG2NUM(max));
                    rb_ary_store(rb_ary_entry(point, m), n, INT2FIX(1));
                  } else {
                    rb_ary_store(rb_ary_entry(score, m), n, LONG2NUM(max));
                    rb_ary_store(rb_ary_entry(point, m), n, INT2FIX(0));
                  }
                } else {
                  if (del > 0) {
                    max = del;
                    rb_ary_store(rb_ary_entry(score, m), n, LONG2NUM(max));
                    rb_ary_store(rb_ary_entry(point, m), n, INT2FIX(2));
                  } else {
                    rb_ary_store(rb_ary_entry(score, m), n, LONG2NUM(max));
                    rb_ary_store(rb_ary_entry(point, m), n, INT2FIX(0));
                  }
                }
              }

              if (max >= max_s) {
                max_m = m;
                max_n = n;
                max_s = max;
              }
            }
          }

          VALUE args[6];
          args[0] = stp;
          args[1] = seq;
          args[2] = score;
          args[3] = point;
          args[4] = LONG2NUM(max_m);
          args[5] = LONG2NUM(max_n);

          return rb_class_new_instance(6, args, rb_path2class("Gigue::ProfileSequenceLocalAlignmentLinearGap"));
        }
      EOCPP
    end

    def local_alignment_linear_gap_rb(options={})
      opts = {
        :gap_del => 100,
        :gap_ins => 100
      }.merge!(options)

      gdel  = opts[:gap_del]
      gins  = opts[:gap_ins]
      pss   = @structural_profile.positions
      plen  = @structural_profile.length
      aas   = @sequence.amino_acids
      slen  = @sequence.length
      max_s = 0
      max_m = nil
      max_n = nil

      # 1. create score and point matrices
      score = Array.new(slen+1)
      point = Array.new(slen+1)

      # 2. initialize score and point matrices
      (0..slen).each do |m|
        score[m] = Array.new(plen+1)
        point[m] = Array.new(plen+1)

        (0..plen).each do |n|
          if    (n == 0 && m == 0)
            score[m][n] = 0
            point[m][n] = NONE
          elsif (n >  0 && m == 0)
            score[m][n] = -gdel*n
            point[m][n] = LEFT
          elsif (n == 0 && m >  0)
            score[m][n] = -gins*m
            point[m][n] = UP
          end
        end
      end

      # 3. fill in score and point matrices
      (1..slen).each do |m|
        (1..plen).each do |n|
          mat = score[m-1][n-1] + pss[n-1].mat_score(aas[m-1])
          del = score[m][n-1] - gdel
          ins = score[m-1][n] - gins
          max = [0, mat, del, ins].max
          score[m][n] = max
          point[m][n] = case max
                        when 0    then NONE
                        when mat  then DIAG
                        when del  then LEFT
                        when ins  then UP
                        end
          if max >= max_s
            max_m = m
            max_n = n
            max_s = max
          end
        end
      end

      ProfileSequenceLocalAlignmentLinearGap.new(@structural_profile, @sequence,
                                                 score, point, max_m, max_n)
    end

    inline(:C) do |builder|
      builder.add_compile_flags '-x c++', '-lstdc++'
      builder.include '<limits>'
      builder.c_raw %q{
        static VALUE local_alignment_affine_gap_cpp(int argc, VALUE *argv, VALUE self) {
          VALUE seq   = rb_iv_get(self, "@sequence");
          VALUE stp   = rb_iv_get(self, "@structural_profile");
          VALUE aas   = rb_funcall(seq, rb_intern("amino_acids"), 0);
          VALUE pss   = rb_funcall(stp, rb_intern("positions"), 0);
          long plen   = NUM2LONG(rb_funcall(stp, rb_intern("length"), 0));
          long slen   = NUM2LONG(rb_funcall(seq, rb_intern("length"), 0));
          long max_s  = 0;
          long max_m  = 0;
          long max_n  = 0;

          // 1. create score, point, jump matrices for match, deletion, and insertion
          VALUE mat_score = rb_ary_new2(slen+1);
          VALUE mat_point = rb_ary_new2(slen+1);
          VALUE mat_jump  = rb_ary_new2(slen+1);

          VALUE del_score = rb_ary_new2(slen+1);
          VALUE del_point = rb_ary_new2(slen+1);
          VALUE del_jump  = rb_ary_new2(slen+1);

          VALUE ins_score = rb_ary_new2(slen+1);
          VALUE ins_point = rb_ary_new2(slen+1);
          VALUE ins_jump  = rb_ary_new2(slen+1);

          // 2. initialize match matrices
          long prev_gap_ins_ext = 0;

          for (long m = 0; m < slen+1; m++) {
            VALUE mat_score_row = rb_ary_new2(plen+1);
            VALUE mat_point_row = rb_ary_new2(plen+1);
            VALUE mat_jump_row  = rb_ary_new2(plen+1);

            rb_ary_store(mat_score, m, mat_score_row);
            rb_ary_store(mat_point, m, mat_point_row);
            rb_ary_store(mat_jump,  m, mat_jump_row);

            for (long n = 0; n < plen+1; n++) {
              if ((m == 0) && (n == 0)) {
                // NONE
                rb_ary_store(mat_score_row, n, LONG2NUM(0));
                rb_ary_store(mat_point_row, n, INT2FIX(0));
              } else if ((m == 0) && (n == 1)) {
                // LEFT
                VALUE ps  = rb_ary_entry(pss, n-1);
                long gap  = NUM2LONG(rb_funcall(ps, rb_intern("gap_del_open"), 0));
                rb_ary_store(mat_score_row, n, LONG2NUM(-gap));
                rb_ary_store(mat_point_row, n, INT2FIX(2));
              } else if ((m == 0) && (n > 1)) {
                // LEFT
                long prv  = NUM2LONG(rb_ary_entry(rb_ary_entry(mat_score, m), n-1));
                VALUE ps  = rb_ary_entry(pss, n-1);
                long gap  = NUM2LONG(rb_funcall(ps, rb_intern("gap_del_ext"), 0));
                rb_ary_store(mat_score_row, n, LONG2NUM(prv-gap));
                rb_ary_store(mat_point_row, n, INT2FIX(2));
              } else if ((m == 1) && (n == 0)) {
                // UP
                VALUE ps  = rb_ary_entry(pss, n);
                long gap  = NUM2LONG(rb_funcall(ps, rb_intern("gap_ins_open"), 0));
                rb_ary_store(mat_score_row, n, LONG2NUM(-gap));
                rb_ary_store(mat_point_row, n, INT2FIX(1));
                prev_gap_ins_ext = NUM2LONG(rb_funcall(ps, rb_intern("gap_ins_ext"), 0));
              } else if ((m > 1) && (n == 0)) {
                // UP
                long prv = NUM2LONG(rb_ary_entry(rb_ary_entry(mat_score, m-1), n));
                rb_ary_store(mat_score_row, n, LONG2NUM(prv-prev_gap_ins_ext));
                rb_ary_store(mat_point_row, n, INT2FIX(1));
              }
            }
          }

          // 3. initialize deletion and insertion matrices
          long n_infinity = std::numeric_limits<long>::min() / 2;

          for (long m = 0; m < slen+1; m++) {
            VALUE del_score_row = rb_ary_new2(plen+1);
            VALUE del_point_row = rb_ary_new2(plen+1);
            VALUE del_jump_row  = rb_ary_new2(plen+1);

            VALUE ins_score_row = rb_ary_new2(plen+1);
            VALUE ins_point_row = rb_ary_new2(plen+1);
            VALUE ins_jump_row  = rb_ary_new2(plen+1);

            rb_ary_store(del_score, m, del_score_row);
            rb_ary_store(del_point, m, del_point_row);
            rb_ary_store(del_jump,  m, del_jump_row);

            rb_ary_store(ins_score, m, ins_score_row);
            rb_ary_store(ins_point, m, ins_point_row);
            rb_ary_store(ins_jump,  m, ins_jump_row);

            // fill in  score matrices[0][0] with -infinity
            for (long n = 0; n < plen+1; n++) {
              if ((m == 0) || (n == 0)) {
                rb_ary_store(del_score_row, n, LONG2NUM(n_infinity));
                rb_ary_store(ins_score_row, n, LONG2NUM(n_infinity));
              }
            }
          }

          // 4. fill in match, deletion, and insertion matrices
          prev_gap_ins_ext = NUM2LONG(rb_funcall(rb_ary_entry(pss, 0), rb_intern("gap_ins_ext"), 0));

          for (long m = 1; m < slen+1; m++) {
            for (long n = 1; n < plen+1; n++) {
              VALUE aa = rb_ary_entry(aas, m-1);
              VALUE ps = rb_ary_entry(pss, n-1);

              long cur_mat_score = NUM2LONG(rb_funcall(ps, rb_intern("mat_score"), 1, aa));
              long prv_mat_score = NUM2LONG(rb_ary_entry(rb_ary_entry(mat_score, m-1), n-1));
              long prv_del_score = NUM2LONG(rb_ary_entry(rb_ary_entry(del_score, m-1), n-1));
              long prv_ins_score = NUM2LONG(rb_ary_entry(rb_ary_entry(ins_score, m-1), n-1));

              long mat_mat = prv_mat_score + cur_mat_score;
              long mat_del = prv_del_score + cur_mat_score;
              long mat_ins = prv_ins_score + cur_mat_score;
              long mat_max = 0;

              if (mat_mat >= mat_ins) {
                if (mat_mat >= mat_del) {
                  if (mat_mat > 0) {
                    // POINT DIAG
                    mat_max = mat_mat;
                    rb_ary_store(rb_ary_entry(mat_score, m), n, LONG2NUM(mat_max));
                    rb_ary_store(rb_ary_entry(mat_point, m), n, INT2FIX(3));
                  } else {
                    // NONE
                    rb_ary_store(rb_ary_entry(mat_score, m), n, LONG2NUM(0));
                    rb_ary_store(rb_ary_entry(mat_point, m), n, INT2FIX(0));
                  }
                } else {
                  if (mat_del > 0) {
                    // JUMP to DEL
                    mat_max = mat_del;
                    rb_ary_store(rb_ary_entry(mat_score, m), n, LONG2NUM(mat_max));
                    rb_ary_store(rb_ary_entry(mat_jump, m), n, rb_ary_new3(3, rb_str_new2("D"), LONG2NUM(m-1), LONG2NUM(n-1)));
                  } else {
                    // NONE
                    rb_ary_store(rb_ary_entry(mat_score, m), n, LONG2NUM(0));
                    rb_ary_store(rb_ary_entry(mat_point, m), n, INT2FIX(0));
                  }
                }
              } else {
                if (mat_ins > mat_del) {
                  if (mat_ins > 0) {
                    // JUMP to INS
                    mat_max = mat_ins;
                    rb_ary_store(rb_ary_entry(mat_score, m), n, LONG2NUM(mat_max));
                    rb_ary_store(rb_ary_entry(mat_jump, m), n, rb_ary_new3(3, rb_str_new2("I"), LONG2NUM(m-1), LONG2NUM(n-1)));
                  } else {
                    // NONE
                    rb_ary_store(rb_ary_entry(mat_score, m), n, LONG2NUM(0));
                    rb_ary_store(rb_ary_entry(mat_point, m), n, INT2FIX(0));
                  }
                } else {
                  if (mat_del > 0) {
                    // JUMP to DEL
                    mat_max = mat_del;
                    rb_ary_store(rb_ary_entry(mat_score, m), n, LONG2NUM(mat_max));
                    rb_ary_store(rb_ary_entry(mat_jump, m), n, rb_ary_new3(3, rb_str_new2("D"), LONG2NUM(m-1), LONG2NUM(n-1)));
                  } else {
                    // NONE
                    rb_ary_store(rb_ary_entry(mat_score, m), n, LONG2NUM(0));
                    rb_ary_store(rb_ary_entry(mat_point, m), n, INT2FIX(0));
                  }
                }
              }

              prv_mat_score     = NUM2LONG(rb_ary_entry(rb_ary_entry(mat_score, m), n-1));
              prv_del_score     = NUM2LONG(rb_ary_entry(rb_ary_entry(del_score, m), n-1));
              long gap_del_open = NUM2LONG(rb_funcall(ps, rb_intern("gap_del_open"), 0));
              long gap_del_ext  = NUM2LONG(rb_funcall(ps, rb_intern("gap_del_ext"), 0));

              long del_mat = prv_mat_score - gap_del_open;
              long del_del = prv_del_score - gap_del_ext;
              long del_max = 0;

              if (del_mat >= del_del) {
                if (del_mat > 0) {
                  // JUMP to MAT
                  del_max = del_mat;
                  rb_ary_store(rb_ary_entry(del_score, m), n, LONG2NUM(del_max));
                  rb_ary_store(rb_ary_entry(del_jump, m), n, rb_ary_new3(3, rb_str_new2("M"), LONG2NUM(m), LONG2NUM(n-1)));
                } else {
                  // NONE
                  rb_ary_store(rb_ary_entry(del_score, m), n, LONG2NUM(0));
                  rb_ary_store(rb_ary_entry(del_point, m), n, INT2FIX(0));
                }
              } else {
                if (del_del > 0) {
                  // POINT LEFT
                  del_max = del_del;
                  rb_ary_store(rb_ary_entry(del_score, m), n, LONG2NUM(del_max));
                  rb_ary_store(rb_ary_entry(del_point, m), n, INT2FIX(2));
                } else {
                  // NONE
                  rb_ary_store(rb_ary_entry(del_score, m), n, LONG2NUM(0));
                  rb_ary_store(rb_ary_entry(del_point, m), n, INT2FIX(0));
                }
              }

              prv_mat_score     = NUM2LONG(rb_ary_entry(rb_ary_entry(mat_score, m-1), n));
              prv_ins_score     = NUM2LONG(rb_ary_entry(rb_ary_entry(ins_score, m-1), n));
              long gap_ins_open = NUM2LONG(rb_funcall(ps, rb_intern("gap_ins_open"), 0));

              long ins_mat = prv_mat_score - gap_ins_open;
              long ins_ins = prv_ins_score - prev_gap_ins_ext;
              long ins_max = 0;

              if (ins_mat >= ins_ins) {
                if (ins_mat > 0) {
                  // JUMP to MAT
                  ins_max = ins_mat;
                  rb_ary_store(rb_ary_entry(ins_score, m), n, LONG2NUM(ins_max));
                  rb_ary_store(rb_ary_entry(ins_jump, m), n, rb_ary_new3(3, rb_str_new2("M"), LONG2NUM(m-1), LONG2NUM(n)));
                } else {
                  // NONE
                  rb_ary_store(rb_ary_entry(ins_score, m), n, LONG2NUM(0));
                  rb_ary_store(rb_ary_entry(ins_point, m), n, INT2FIX(0));
                }
              } else {
                if (ins_ins > 0) {
                  // POINT UP
                  ins_max = ins_ins;
                  rb_ary_store(rb_ary_entry(ins_score, m), n, LONG2NUM(ins_max));
                  rb_ary_store(rb_ary_entry(ins_point, m), n, INT2FIX(1));
                  prev_gap_ins_ext = NUM2LONG(rb_funcall(ps, rb_intern("gap_ins_ext"), 0));
                } else {
                  // NONE
                  rb_ary_store(rb_ary_entry(ins_score, m), n, LONG2NUM(0));
                  rb_ary_store(rb_ary_entry(ins_point, m), n, INT2FIX(0));
                }
              }

              long max = 0;

              if (mat_max >= ins_max) {
                if (mat_max >= del_max) {
                  max = mat_max;
                } else {
                  max = del_max;
                }
              } else {
                if (ins_max > del_max) {
                  max = ins_max;
                } else {
                  max = del_max;
                }
              }

              if (max >= max_s) {
                max_m = m;
                max_n = n;
                max_s = max;
              }
            }
          }

          VALUE args[13];
          args[0]   = stp;
          args[1]   = seq;
          args[2]   = mat_score;
          args[3]   = mat_point;
          args[4]   = mat_jump;
          args[5]   = del_score;
          args[6]   = del_point;
          args[7]   = del_jump;
          args[8]   = ins_score;
          args[9]   = ins_point;
          args[10]  = ins_jump;
          args[11]  = LONG2NUM(max_m);
          args[12]  = LONG2NUM(max_n);

          return rb_class_new_instance(13, args, rb_path2class("Gigue::ProfileSequenceLocalAlignmentAffineGap"));
        }
      }
    end

    def local_alignment_affine_gap_rb
      pss     = @structural_profile.positions
      plen    = @structural_profile.length
      aas     = @sequence.amino_acids
      slen    = @sequence.length
      max_s   = 0
      max_m   = nil
      max_n   = nil

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

      (0..slen).each do |m|
        mat_score[m] = Array.new(plen+1)
        mat_point[m] = Array.new(plen+1)
        mat_jump[m]  = Array.new(plen+1)

        (0..plen).each do |n|
          if    (n == 0 && m == 0)
            mat_score[m][n] = 0
            mat_point[m][n] = NONE
          elsif (n == 1 && m == 0)
            mat_score[m][n] = -pss[n-1].gap_del_open
            mat_point[m][n] = LEFT
          elsif (n >  1 && m == 0)
            mat_score[m][n] = mat_score[m][n-1] - pss[n-1].gap_del_ext
            mat_point[m][n] = LEFT
          elsif (n == 0 && m == 1)
            mat_score[m][n] = -pss[n].gap_ins_open
            mat_point[m][n] = UP
            prev_gap_ins_ext  = pss[n].gap_ins_ext
          elsif (n == 0 && m >  1)
            mat_score[m][n] = mat_score[m-1][n] - prev_gap_ins_ext
            mat_point[m][n] = UP
          end
        end
      end

      # 3. initialize deletion and insertion matrices
      infinity = 1/0.0

      (0..slen).each do |m|
        del_score[m] = Array.new(plen+1)
        del_point[m] = Array.new(plen+1)
        del_jump[m]  = Array.new(plen+1)

        ins_score[m] = Array.new(plen+1)
        ins_point[m] = Array.new(plen+1)
        ins_jump[m]  = Array.new(plen+1)

        (0..plen).each do |n|
          del_score[m][n] = -infinity if (m == 0 || n == 0)
          ins_score[m][n] = -infinity if (m == 0 || n == 0)
        end
      end

      # 4. fill in match, deletion, and insertion matrices
      log_fmt = "%-12s : %5s : %12s"
      prev_gap_ins_ext = pss[0].gap_ins_ext

      (1..slen).each do |m|
        (1..plen).each do |n|
          mat_mat = mat_score[m-1][n-1] + pss[n-1].mat_score(aas[m-1])
          mat_del = del_score[m-1][n-1] + pss[n-1].mat_score(aas[m-1])
          mat_ins = ins_score[m-1][n-1] + pss[n-1].mat_score(aas[m-1])
          mat_max = [0, mat_mat, mat_del, mat_ins].max
          mat_score[m][n] = mat_max
          log_mat = "M[#{m}][#{n}]"

          case mat_max
          when 0
            mat_point[m][n] = NONE
            $logger.debug log_fmt % [log_mat, "POINT", "NONE"]
          when mat_mat
            mat_point[m][n] = DIAG
            $logger.debug log_fmt % [log_mat, "POINT", "M[#{m-1}][#{n-1}]"]
          when mat_del
            mat_jump[m][n] = ['D', m-1, n-1]
            $logger.debug log_fmt % [log_mat, "JUMP", "D[#{m-1}][#{n-1}]"]
          when mat_ins
            mat_jump[m][n] = ['I', m-1, n-1]
            $logger.debug log_fmt % [log_mat, "JUMP", "I[#{m-1}][#{n-1}]"]
          end

          del_mat = mat_score[m][n-1] - pss[n-1].gap_del_open
          del_del = del_score[m][n-1] - pss[n-1].gap_del_ext
          del_max = [0, del_mat, del_del].max
          del_score[m][n] = del_max
          log_mat = "D[#{m}][#{n}]"

          case del_max
          when 0
            del_point[m][n] = NONE
            $logger.debug log_fmt % [log_mat, "POINT", "NONE"]
          when del_mat
            del_jump[m][n] = ['M', m, n-1]
            $logger.debug log_fmt % [log_mat, "JUMP", "M[#{m}][#{n-1}]"]
          when del_del
            del_point[m][n] = LEFT
            $logger.debug log_fmt % [log_mat, "POINT", "D[#{m}][#{n-1}]"]
          end

          ins_mat = mat_score[m-1][n] - pss[n-1].gap_ins_open
          ins_ins = ins_score[m-1][n] - prev_gap_ins_ext
          ins_max = [0, ins_mat, ins_ins].max
          ins_score[m][n] = ins_max
          log_mat = "I[#{m}][#{n}]"

          case ins_max
          when 0
            ins_point[m][n] = NONE
            $logger.debug log_fmt % [log_mat, "POINT", "NONE"]
          when ins_mat
            ins_jump[m][n] = ['M', m-1, n]
            prev_gap_ins_ext = pss[n-1].gap_ins_ext
            $logger.debug log_fmt % [log_mat, "JUMP", "M[#{m-1}][#{n}]"]
          when ins_ins
            ins_point[m][n] = UP
            $logger.debug log_fmt % [log_mat, "POINT", "I[#{m-1}][#{n}]"]
          end

          # keep the record of a matrix and its indexes having maximum value so far
          max = [mat_max, del_max, ins_max].max
          if max >= max_s
            max_m   = m
            max_n   = n
            max_s   = max
          end
        end
      end

      ProfileSequenceLocalAlignmentAffineGap.new(@structural_profile, @sequence,
                                                  mat_score, mat_point, mat_jump,
                                                  del_score, del_point, del_jump,
                                                  ins_score, ins_point, ins_jump,
                                                  max_m, max_n)
    end

    inline(:C) do |builder|
      builder.add_compile_flags '-x c++', '-lstdc++'
      builder.c_raw <<-EOCPP
        static VALUE global_alignment_linear_gap_cpp(int argc, VALUE *argv, VALUE self) {
          VALUE options = argv[0] == Qnil ? rb_hash_new() : argv[0];
          VALUE opts    = rb_hash_new();

          rb_hash_aset(opts, ID2SYM(rb_intern("gap_del")), INT2NUM(100));
          rb_hash_aset(opts, ID2SYM(rb_intern("gap_ins")), INT2NUM(100));
          rb_funcall(opts, rb_intern("merge!"), 1, options);

          long gdel = NUM2LONG(rb_hash_aref(opts, ID2SYM(rb_intern("gap_del"))));
          long gins = NUM2LONG(rb_hash_aref(opts, ID2SYM(rb_intern("gap_ins"))));
          VALUE alg = rb_iv_get(self, "@algorithm");
          VALUE stp = rb_iv_get(self, "@structural_profile");
          VALUE seq = rb_iv_get(self, "@sequence");
          VALUE pss = rb_funcall(stp, rb_intern("positions"), 0);
          VALUE aas = rb_funcall(seq, rb_intern("amino_acids"), 0);
          long plen = NUM2LONG(rb_funcall(stp, rb_intern("length"), 0));
          long slen = NUM2LONG(rb_funcall(seq, rb_intern("length"), 0));


          // 1. create score and point matrices
          VALUE score = rb_ary_new2(slen+1);
          VALUE point = rb_ary_new2(slen+1);

          // 2. initialize score and point matrices
          for (long m = 0; m < slen+1; m++) {
            VALUE score_row = rb_ary_new2(plen+1);
            VALUE point_row = rb_ary_new2(plen+1);

            rb_ary_store(score, m, score_row);
            rb_ary_store(point, m, point_row);

            for (long n = 0; n < plen+1; n++) {
              if ((n == 0) && (m == 0)) {
                // NONE
                rb_ary_store(score_row, n, LONG2NUM(0));
                rb_ary_store(point_row, n, INT2FIX(0));
              } else if (n >  0 && m == 0) {
                // LEFT
                if (alg == ID2SYM(rb_intern("glolocstr"))) {
                  // when GLOCAL and str : seq > 1.5 : 1 -> no terminal deletion penalty
                  rb_ary_store(score_row, n, LONG2NUM(0));
                  rb_ary_store(point_row, n, INT2FIX(2));
                } else  {
                  rb_ary_store(score_row, n, LONG2NUM(-gdel*n));
                  rb_ary_store(point_row, n, INT2FIX(2));
                }
              } else if (n == 0 && m >  0) {
                // UP
                if (alg == ID2SYM(rb_intern("glolocseq"))) {
                  // when GLOCAL and str : seq < 1 : 1.5 -> no terminal insertion penalty
                  rb_ary_store(score_row, n, LONG2NUM(0));
                  rb_ary_store(point_row, n, INT2FIX(1));
                } else { 
                  rb_ary_store(score_row, n, LONG2NUM(-gins*m));
                  rb_ary_store(point_row, n, INT2FIX(1));
                }
              }
            }
          }

          // 3. fill in score and point matrices
          for (long m = 1; m < slen+1; m++) {
            VALUE aa = rb_ary_entry(aas, m-1);
            for (long n = 1; n < plen+1; n++) {
              VALUE ps = rb_ary_entry(pss, n-1);

              long cur = NUM2LONG(rb_funcall(ps, rb_intern("mat_score"), 1, aa));
              long mat = NUM2LONG(rb_ary_entry(rb_ary_entry(score, m-1), n-1)) + cur;
              long del = 0;

              if ((m == slen) && (alg == ID2SYM(rb_intern("glolocstr")))) {
                // when GLOCAL and str : seq > 1.5 : 1 -> no deletion penalty
                del = NUM2LONG(rb_ary_entry(rb_ary_entry(score, m), n-1));
              } else {
                del = NUM2LONG(rb_ary_entry(rb_ary_entry(score, m), n-1)) - gdel;
              }

              long ins = 0;

              if ((n == plen) && (alg == ID2SYM(rb_intern("glolocseq")))) {
                // when GLOCAL and str : seq > 1.5 : 1 -> no deletion penalty
                ins = NUM2LONG(rb_ary_entry(rb_ary_entry(score, m-1), n));
              } else {
                ins = NUM2LONG(rb_ary_entry(rb_ary_entry(score, m-1), n)) - gins;
              }

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

          return rb_class_new_instance(4, args, rb_path2class("Gigue::ProfileSequenceGlobalAlignmentLinearGap"));
        }
      EOCPP
    end

    def global_alignment_linear_gap_rb(options={})
      opts = {
        :gap_del => 100,
        :gap_ins => 100
      }.merge!(options)

      gdel  = opts[:gap_del]
      gins  = opts[:gap_ins]
      pss   = @structural_profile.positions
      plen  = @structural_profile.length
      aas   = @sequence.amino_acids
      slen  = @sequence.length

      # 1. create score and point matrices
      score = Array.new(slen+1)
      point = Array.new(slen+1)

      # 2. initialize score and point matrices
      (0..slen).each do |m|
        score[m] = Array.new(plen+1)
        point[m] = Array.new(plen+1)

        (0..plen).each do |n|
          if    (n == 0 && m == 0)
            score[m][n] = 0
            point[m][n] = NONE
          elsif (n >  0 && m == 0)
            score[m][n] = -gdel*n
            point[m][n] = LEFT
          elsif (n == 0 && m >  0)
            score[m][n] = -gins*m
            point[m][n] = UP
          end
        end
      end

      # 3. fill in score and point matrices
      (1..slen).each do |m|
        (1..plen).each do |n|
          mat = score[m-1][n-1] + pss[n-1].mat_score(aas[m-1])
          del = score[m][n-1] - gdel
          ins = score[m-1][n] - gins
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
      builder.include '<limits>'
      builder.c_raw %q{
        static VALUE global_alignment_affine_gap_cpp(int argc, VALUE *argv, VALUE self) {
          VALUE alg = rb_iv_get(self, "@algorithm");
          VALUE seq = rb_iv_get(self, "@sequence");
          VALUE stp = rb_iv_get(self, "@structural_profile");
          VALUE aas = rb_funcall(seq, rb_intern("amino_acids"), 0);
          VALUE pss = rb_funcall(stp, rb_intern("positions"), 0);
          long plen = NUM2LONG(rb_funcall(stp, rb_intern("length"), 0));
          long slen = NUM2LONG(rb_funcall(seq, rb_intern("length"), 0));

          // 1. create score, point, jump matrices for match, deletion, and insertion
          VALUE mat_score = rb_ary_new2(slen+1);
          VALUE mat_point = rb_ary_new2(slen+1);
          VALUE mat_jump  = rb_ary_new2(slen+1);

          VALUE del_score = rb_ary_new2(slen+1);
          VALUE del_point = rb_ary_new2(slen+1);
          VALUE del_jump  = rb_ary_new2(slen+1);

          VALUE ins_score = rb_ary_new2(slen+1);
          VALUE ins_point = rb_ary_new2(slen+1);
          VALUE ins_jump  = rb_ary_new2(slen+1);

          // 2. initialize match matrices
          long prev_gap_ins_ext = 0;

          for (long m = 0; m < slen+1; m++) {
            VALUE mat_score_row = rb_ary_new2(plen+1);
            VALUE mat_point_row = rb_ary_new2(plen+1);
            VALUE mat_jump_row  = rb_ary_new2(plen+1);

            rb_ary_store(mat_score, m, mat_score_row);
            rb_ary_store(mat_point, m, mat_point_row);
            rb_ary_store(mat_jump,  m, mat_jump_row);

            for (long n = 0; n < plen+1; n++) {
              if ((m == 0) && (n == 0)) {
                // NONE
                rb_ary_store(mat_score_row, n, LONG2NUM(0));
                rb_ary_store(mat_point_row, n, INT2FIX(0));
              } else if ((m == 0) && (n == 1)) {
                // LEFT
                if (alg == ID2SYM(rb_intern("glolocstr"))) {
                  rb_ary_store(mat_score_row, n, LONG2NUM(0));
                  rb_ary_store(mat_point_row, n, INT2FIX(2));
                } else {
                  VALUE ps  = rb_ary_entry(pss, n-1);
                  long gap  = NUM2LONG(rb_funcall(ps, rb_intern("gap_del_open"), 0));
                  rb_ary_store(mat_score_row, n, LONG2NUM(-gap));
                  rb_ary_store(mat_point_row, n, INT2FIX(2));
                }
              } else if ((m == 0) && (n > 1)) {
                // LEFT
                if (alg == ID2SYM(rb_intern("glolocstr"))) {
                  rb_ary_store(mat_score_row, n, LONG2NUM(0));
                  rb_ary_store(mat_point_row, n, INT2FIX(2));
                } else {
                  long prv  = NUM2LONG(rb_ary_entry(rb_ary_entry(mat_score, m), n-1));
                  VALUE ps  = rb_ary_entry(pss, n-1);
                  long gap  = NUM2LONG(rb_funcall(ps, rb_intern("gap_del_ext"), 0));
                  rb_ary_store(mat_score_row, n, LONG2NUM(prv-gap));
                  rb_ary_store(mat_point_row, n, INT2FIX(2));
                }
              } else if ((m == 1) && (n == 0)) {
                // UP
                if (alg == ID2SYM(rb_intern("glolocseq"))) {
                  rb_ary_store(mat_score_row, n, LONG2NUM(0));
                  rb_ary_store(mat_point_row, n, INT2FIX(1));
                } else {
                  VALUE ps  = rb_ary_entry(pss, n);
                  long gap  = NUM2LONG(rb_funcall(ps, rb_intern("gap_ins_open"), 0));
                  rb_ary_store(mat_score_row, n, LONG2NUM(-gap));
                  rb_ary_store(mat_point_row, n, INT2FIX(1));
                  prev_gap_ins_ext = NUM2LONG(rb_funcall(ps, rb_intern("gap_ins_ext"), 0));
                }
              } else if ((m > 1) && (n == 0)) {
                // UP
                if (alg == ID2SYM(rb_intern("glolocseq"))) {
                  rb_ary_store(mat_score_row, n, LONG2NUM(0));
                  rb_ary_store(mat_point_row, n, INT2FIX(1));
                } else {
                  long prv = NUM2LONG(rb_ary_entry(rb_ary_entry(mat_score, m-1), n));
                  rb_ary_store(mat_score_row, n, LONG2NUM(prv-prev_gap_ins_ext));
                  rb_ary_store(mat_point_row, n, INT2FIX(1));
                }
              }
            }
          }

          // 3. initialize deletion and insertion matrices
          long n_infinity = std::numeric_limits<long>::min() / 2;

          for (long m = 0; m < slen+1; m++) {
            VALUE del_score_row = rb_ary_new2(plen+1);
            VALUE del_point_row = rb_ary_new2(plen+1);
            VALUE del_jump_row  = rb_ary_new2(plen+1);

            VALUE ins_score_row = rb_ary_new2(plen+1);
            VALUE ins_point_row = rb_ary_new2(plen+1);
            VALUE ins_jump_row  = rb_ary_new2(plen+1);

            rb_ary_store(del_score, m, del_score_row);
            rb_ary_store(del_point, m, del_point_row);
            rb_ary_store(del_jump,  m, del_jump_row);

            rb_ary_store(ins_score, m, ins_score_row);
            rb_ary_store(ins_point, m, ins_point_row);
            rb_ary_store(ins_jump,  m, ins_jump_row);

            // fill in  score matrices[0][0] with -infinity
            for (long n = 0; n < plen+1; n++) {
              if ((m == 0) || (n == 0)) {
                rb_ary_store(del_score_row, n, LONG2NUM(n_infinity));
                rb_ary_store(ins_score_row, n, LONG2NUM(n_infinity));
              }
            }
          }

          // 4. fill in match, deletion, and insertion matrices
          prev_gap_ins_ext = NUM2LONG(rb_funcall(rb_ary_entry(pss, 0), rb_intern("gap_ins_ext"), 0));

          for (long m = 1; m < slen+1; m++) {
            VALUE aa = rb_ary_entry(aas, m-1);
            for (long n = 1; n < plen+1; n++) {
              VALUE ps = rb_ary_entry(pss, n-1);

              long cur_mat_score = NUM2LONG(rb_funcall(ps, rb_intern("mat_score"), 1, aa));
              long prv_mat_score = NUM2LONG(rb_ary_entry(rb_ary_entry(mat_score, m-1), n-1));
              long prv_del_score = NUM2LONG(rb_ary_entry(rb_ary_entry(del_score, m-1), n-1));
              long prv_ins_score = NUM2LONG(rb_ary_entry(rb_ary_entry(ins_score, m-1), n-1));

              long mat_mat = prv_mat_score + cur_mat_score;
              long mat_del = prv_del_score + cur_mat_score;
              long mat_ins = prv_ins_score + cur_mat_score;

              if (mat_mat >= mat_ins) {
                if (mat_mat >= mat_del) {
                  // POINT DIAG
                  rb_ary_store(rb_ary_entry(mat_score, m), n, LONG2NUM(mat_mat));
                  rb_ary_store(rb_ary_entry(mat_point, m), n, INT2FIX(3));
                } else {
                  // JUMP to DEL
                  rb_ary_store(rb_ary_entry(mat_score, m), n, LONG2NUM(mat_del));
                  rb_ary_store(rb_ary_entry(mat_jump, m), n, rb_ary_new3(3, rb_str_new2("D"), LONG2NUM(m-1), LONG2NUM(n-1)));
                }
              } else {
                if (mat_ins > mat_del) {
                  // JUMP to INS
                  rb_ary_store(rb_ary_entry(mat_score, m), n, LONG2NUM(mat_ins));
                  rb_ary_store(rb_ary_entry(mat_jump, m), n, rb_ary_new3(3, rb_str_new2("I"), LONG2NUM(m-1), LONG2NUM(n-1)));
                } else {
                  // JUMP to DEL
                  rb_ary_store(rb_ary_entry(mat_score, m), n, LONG2NUM(mat_del));
                  rb_ary_store(rb_ary_entry(mat_jump, m), n, rb_ary_new3(3, rb_str_new2("D"), LONG2NUM(m-1), LONG2NUM(n-1)));
                }
              }

              prv_mat_score     = NUM2LONG(rb_ary_entry(rb_ary_entry(mat_score, m), n-1));
              prv_del_score     = NUM2LONG(rb_ary_entry(rb_ary_entry(del_score, m), n-1));
              long gap_del_open = 0;
              long gap_del_ext  = 0;

              if ((m != slen) || (alg != ID2SYM(rb_intern("glolocstr")))) {
                gap_del_open = NUM2LONG(rb_funcall(ps, rb_intern("gap_del_open"), 0));
                gap_del_ext  = NUM2LONG(rb_funcall(ps, rb_intern("gap_del_ext"), 0));
              }

              long del_mat = prv_mat_score - gap_del_open;
              long del_del = prv_del_score - gap_del_ext;

              if (del_mat >= del_del) {
                // JUMP to MAT
                rb_ary_store(rb_ary_entry(del_score, m), n, LONG2NUM(del_mat));
                rb_ary_store(rb_ary_entry(del_jump, m), n, rb_ary_new3(3, rb_str_new2("M"), LONG2NUM(m), LONG2NUM(n-1)));
              } else {
                // POINT LEFT
                rb_ary_store(rb_ary_entry(del_score, m), n, LONG2NUM(del_del));
                rb_ary_store(rb_ary_entry(del_point, m), n, INT2FIX(2));
              }

              prv_mat_score     = NUM2LONG(rb_ary_entry(rb_ary_entry(mat_score, m-1), n));
              prv_ins_score     = NUM2LONG(rb_ary_entry(rb_ary_entry(ins_score, m-1), n));
              long gap_ins_open = 0;

              if ((n != plen) || (alg != ID2SYM(rb_intern("glolocseq")))) {
                gap_ins_open = NUM2LONG(rb_funcall(ps, rb_intern("gap_ins_open"), 0));
              } else {
                prev_gap_ins_ext = 0;
              }

              long ins_mat = prv_mat_score - gap_ins_open;
              long ins_ins = prv_ins_score - prev_gap_ins_ext;

              if (ins_mat >= ins_ins) {
                // JUMP to MAT
                rb_ary_store(rb_ary_entry(ins_score, m), n, LONG2NUM(ins_mat));
                rb_ary_store(rb_ary_entry(ins_jump, m), n, rb_ary_new3(3, rb_str_new2("M"), LONG2NUM(m-1), LONG2NUM(n)));
                prev_gap_ins_ext = NUM2LONG(rb_funcall(ps, rb_intern("gap_ins_ext"), 0));
              } else {
                // POINT UP
                rb_ary_store(rb_ary_entry(ins_score, m), n, LONG2NUM(ins_ins));
                rb_ary_store(rb_ary_entry(ins_point, m), n, INT2FIX(1));
              }
            }
          }

          VALUE args[11];
          args[0]   = stp;
          args[1]   = seq;
          args[2]   = mat_score;
          args[3]   = mat_point;
          args[4]   = mat_jump;
          args[5]   = del_score;
          args[6]   = del_point;
          args[7]   = del_jump;
          args[8]   = ins_score;
          args[9]   = ins_point;
          args[10]  = ins_jump;

          return rb_class_new_instance(11, args, rb_path2class("Gigue::ProfileSequenceGlobalAlignmentAffineGap"));
        }
      }
    end

    def global_alignment_affine_gap_rb
      pss   = @structural_profile.positions
      plen  = @structural_profile.length
      aas   = @sequence.amino_acids
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

      (0..slen).each do |m|
        mat_score[m] = Array.new(plen+1)
        mat_point[m] = Array.new(plen+1)
        mat_jump[m]  = Array.new(plen+1)

        (0..plen).each do |n|
          if    (n == 0 && m == 0)
            mat_score[m][n] = 0
            mat_point[m][n] = NONE
          elsif (n == 1 && m == 0)
            mat_score[m][n] = -pss[n-1].gap_del_open
            mat_point[m][n] = LEFT
          elsif (n >  1 && m == 0)
            mat_score[m][n] = mat_score[m][n-1] - pss[n-1].gap_del_ext
            mat_point[m][n] = LEFT
          elsif (n == 0 && m == 1)
            mat_score[m][n] = -pss[n].gap_ins_open
            mat_point[m][n] = UP
            prev_gap_ins_ext  = pss[n].gap_ins_ext
          elsif (n == 0 && m >  1)
            mat_score[m][n] = mat_score[m-1][n] - prev_gap_ins_ext
            mat_point[m][n] = UP
          end
        end
      end

      # 3. initialize deletion and insertion matrices
      infinity = 1/0.0

      (0..slen).each do |m|
        del_score[m] = Array.new(plen+1)
        del_point[m] = Array.new(plen+1)
        del_jump[m]  = Array.new(plen+1)

        ins_score[m] = Array.new(plen+1)
        ins_point[m] = Array.new(plen+1)
        ins_jump[m]  = Array.new(plen+1)

        (0..plen).each do |n|
          del_score[m][n] = -infinity if (m == 0 || n == 0)
          ins_score[m][n] = -infinity if (m == 0 || n == 0)
        end
      end

      # 4. fill in match, deletion, and insertion matrices
      log_fmt = "%-12s : %5s : %12s"
      prev_gap_ins_ext = pss[0].gap_ins_ext

      (1..slen).each do |m|
        (1..plen).each do |n|
          mat_mat = mat_score[m-1][n-1] + pss[n-1].mat_score(aas[m-1])
          mat_del = del_score[m-1][n-1] + pss[n-1].mat_score(aas[m-1])
          mat_ins = ins_score[m-1][n-1] + pss[n-1].mat_score(aas[m-1])
          mat_max = [mat_mat, mat_del, mat_ins].max
          mat_score[m][n] = mat_max
          log_mat = "M[#{m}][#{n}]"

          case mat_max
          when mat_mat
            mat_point[m][n] = DIAG
            $logger.debug log_fmt % [log_mat, "POINT", "M[#{m-1}][#{n-1}]"]
          when mat_del
            mat_jump[m][n] = ['D', m-1, n-1]
            $logger.debug log_fmt % [log_mat, "JUMP", "D[#{m-1}][#{n-1}]"]
          when mat_ins
            mat_jump[m][n] = ['I', m-1, n-1]
            $logger.debug log_fmt % [log_mat, "JUMP", "I[#{m-1}][#{n-1}]"]
          end

          del_mat = mat_score[m][n-1] - pss[n-1].gap_del_open
          del_del = del_score[m][n-1] - pss[n-1].gap_del_ext
          del_max = [del_mat, del_del].max
          del_score[m][n] = del_max
          log_mat = "D[#{m}][#{n}]"

          case del_max
          when del_mat
            del_jump[m][n] = ['M', m, n-1]
            $logger.debug log_fmt % [log_mat, "JUMP", "M[#{m}][#{n-1}]"]
          when del_del
            del_point[m][n] = LEFT
            $logger.debug log_fmt % [log_mat, "POINT", "D[#{m}][#{n-1}]"]
          end

          ins_mat = mat_score[m-1][n] - pss[n-1].gap_ins_open
          ins_ins = ins_score[m-1][n] - prev_gap_ins_ext
          ins_max = [ins_mat, ins_ins].max
          ins_score[m][n] = ins_max
          log_mat = "I[#{m}][#{n}]"

          case ins_max
          when ins_mat
            ins_jump[m][n] = ['M', m-1, n]
            prev_gap_ins_ext = pss[n-1].gap_ins_ext
            $logger.debug log_fmt % [log_mat, "JUMP", "M[#{m-1}][#{n}]"]
          when ins_ins
            ins_point[m][n] = UP
            $logger.debug log_fmt % [log_mat, "POINT", "I[#{m-1}][#{n}]"]
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
