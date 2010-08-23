module Gigue
  class SequenceSequenceAligner

    attr_reader :sequence1, :sequence2

    def initialize(seq1, seq2)
      @sequence1 = seq1
      @sequence2 = seq2
    end

    def local_alignment_linear_gap(options={})
      begin
        local_alignment_linear_gap_cpp(options)
      rescue
        local_alignment_linear_gap_rb(options)
      end
    end

    def local_alignment_affine_gap
      begin
        local_alignment_affine_gap_cpp
      rescue
        local_alignment_affine_gap_rb
      end
    end

    def global_alignment_linear_gap
      begin
        global_alignment_linear_gap_cpp(options={})
      rescue
        global_alignment_linear_gap_rb(options={})
      end
    end

    def global_alignment_affine_gap
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
          VALUE klass   = rb_const_get(rb_cObject, rb_intern("Gigue"));
          VALUE klass2  = rb_const_get(klass, rb_intern("SubstitutionTable"));
          VALUE subst   = rb_funcall(klass2, rb_intern("blosum62"), 0);

          VALUE options = argv[0] == Qnil ? rb_hash_new() : argv[0];
          VALUE opts    = rb_hash_new();

          rb_hash_aset(opts, ID2SYM(rb_intern("subst")), subst);
          rb_hash_aset(opts, ID2SYM(rb_intern("gap")), INT2NUM(8));
          rb_funcall(opts, rb_intern("merge!"), 1, options);

          long gap    = NUM2LONG(rb_hash_aref(opts, ID2SYM(rb_intern("gap"))));
          VALUE seq1  = rb_iv_get(self, "@sequence1");
          VALUE seq2  = rb_iv_get(self, "@sequence2");
          VALUE aas1  = rb_funcall(seq1, rb_intern("amino_acids"), 0);
          VALUE aas2  = rb_funcall(seq2, rb_intern("amino_acids"), 0);
          long slen1  = NUM2LONG(rb_funcall(seq1, rb_intern("length"), 0));
          long slen2  = NUM2LONG(rb_funcall(seq2, rb_intern("length"), 0));
          long max_s  = 0;
          long max_m  = 0;
          long max_n  = 0;

          // 1. create score and point matrices
          VALUE score = rb_ary_new2(slen2+1);
          VALUE point = rb_ary_new2(slen2+1);

          // 2. initialize score and point matrices
          for (long m = 0; m < slen2+1; m++) {
            VALUE score_row = rb_ary_new2(slen1+1);
            VALUE point_row = rb_ary_new2(slen1+1);

            rb_ary_store(score, m, score_row);
            rb_ary_store(point, m, point_row);

            for (long n = 0; n < slen1+1; n++) {
              if ((n == 0) && (m == 0)) {
                // NONE
                rb_ary_store(score_row, n, INT2NUM(0));
                rb_ary_store(point_row, n, INT2FIX(0));
              } else if ((n > 0) && (m == 0)) {
                // LEFT
                rb_ary_store(score_row, n, LONG2NUM(-gap*n));
                rb_ary_store(point_row, n, INT2FIX(2));
              } else if ((n == 0) && (m > 0)) {
                // UP
                rb_ary_store(score_row, n, LONG2NUM(-gap*m));
                rb_ary_store(point_row, n, INT2FIX(1));
              }
            }
          }

          // 3. fill in score and point matrices
          for (long m = 1; m < slen2+1; m++) {
            VALUE aa2 = rb_ary_entry(aas2, m-1);
            for (long n = 1; n < slen1+1; n++) {
              VALUE aa1 = rb_ary_entry(aas1, n-1);
              
              long cur_mat_score = NUM2LONG(rb_funcall(subst, rb_intern("score"), 2, aa1, aa2));
              long mat = NUM2LONG(rb_ary_entry(rb_ary_entry(score, m-1), n-1)) + cur_mat_score;
              long del = NUM2LONG(rb_ary_entry(rb_ary_entry(score, m), n-1)) - gap;
              long ins = NUM2LONG(rb_ary_entry(rb_ary_entry(score, m-1), n)) - gap;
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
          args[0] = seq1;
          args[1] = seq2;
          args[2] = score;
          args[3] = point;
          args[4] = LONG2NUM(max_m);
          args[5] = LONG2NUM(max_n);

          return rb_class_new_instance(6, args, rb_path2class("SequenceSequenceLocalAlignmentLinearGap"));
        }
      EOCPP
    end

    def local_alignment_linear_gap_rb(options={})
      opts = {
        :subst  => SubstitutionTable::blosum62,
        :gap    => 8
      }.merge!(options)

      subst = opts[:subst]
      gap   = opts[:gap]
      aas1  = @sequence1.amino_acids
      slen1 = @sequence1.length
      aas2  = @sequence2.amino_acids
      slen2 = @sequence2.length
      max_s = 0
      max_m = nil
      max_n = nil

      # 1. create score and point matrices
      score = Array.new(slen2+1)
      point = Array.new(slen2+1)

      # 2. initialize score and point matrices
      (0..slen2).each do |m|
        score[m] = Array.new(slen1+1)
        point[m] = Array.new(slen1+1)

        (0..slen1).each do |n|
          if    (n == 0 && m == 0)
            score[m][n] = 0
            point[m][n] = NONE
          elsif (n >  0 && m == 0)
            score[m][n] = -gap*n
            point[m][n] = LEFT
          elsif (n == 0 && m >  0)
            score[m][n] = -gap*m
            point[m][n] = UP
          end
        end
      end

      # 3. fill in score and point matrices
      (1..slen2).each do |m|
        (1..slen1).each do |n|
          mat = score[m-1][n-1] + subst.score(aas1[n-1], aas2[m-1])
          del = score[m][n-1] - gap
          ins = score[m-1][n] - gap
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

      SequenceSequenceLocalAlignmentLinearGap.new(@sequence1, @sequence2,
                                                  score, point, max_m, max_n)
    end

    inline(:C) do |builder|
      builder.add_compile_flags '-x c++', '-lstdc++'
      builder.include '<limits>'
      builder.c_raw %q{
        static VALUE local_alignment_affine_gap_cpp(int argc, VALUE *argv, VALUE self) {
          VALUE klass   = rb_const_get(rb_cObject, rb_intern("Gigue"));
          VALUE klass2  = rb_const_get(klass, rb_intern("SubstitutionTable"));
          VALUE subst   = rb_funcall(klass2, rb_intern("blosum62"), 0);

          VALUE options = argv[0] == Qnil ? rb_hash_new() : argv[0];
          VALUE opts    = rb_hash_new();

          rb_hash_aset(opts, ID2SYM(rb_intern("subst")), subst);
          rb_hash_aset(opts, ID2SYM(rb_intern("gopen")), DBL2NUM(10.0));
          rb_hash_aset(opts, ID2SYM(rb_intern("gext")), DBL2NUM(0.5));
          rb_funcall(opts, rb_intern("merge!"), 1, options);

          double gopen  = NUM2DBL(rb_hash_aref(opts, ID2SYM(rb_intern("gopen"))));
          double gext   = NUM2DBL(rb_hash_aref(opts, ID2SYM(rb_intern("gext"))));

          VALUE seq1  = rb_iv_get(self, "@sequence1");
          VALUE seq2  = rb_iv_get(self, "@sequence2");
          VALUE aas1  = rb_funcall(seq1, rb_intern("amino_acids"), 0);
          VALUE aas2  = rb_funcall(seq2, rb_intern("amino_acids"), 0);
          long slen1  = NUM2LONG(rb_funcall(seq1, rb_intern("length"), 0));
          long slen2  = NUM2LONG(rb_funcall(seq2, rb_intern("length"), 0));
          double max_s  = 0;
          long max_m    = 0;
          long max_n    = 0;

          // 1. create score, point, jump matrices for match, deletion, and insertion
          VALUE mat_score = rb_ary_new2(slen2+1);
          VALUE mat_point = rb_ary_new2(slen2+1);
          VALUE mat_jump  = rb_ary_new2(slen2+1);

          VALUE del_score = rb_ary_new2(slen2+1);
          VALUE del_point = rb_ary_new2(slen2+1);
          VALUE del_jump  = rb_ary_new2(slen2+1);

          VALUE ins_score = rb_ary_new2(slen2+1);
          VALUE ins_point = rb_ary_new2(slen2+1);
          VALUE ins_jump  = rb_ary_new2(slen2+1);

          // 2. initialize match matrices
          for (long m = 0; m < slen2+1; m++) {
            VALUE mat_score_row = rb_ary_new2(slen1+1);
            VALUE mat_point_row = rb_ary_new2(slen1+1);
            VALUE mat_jump_row  = rb_ary_new2(slen1+1);

            rb_ary_store(mat_score, m, mat_score_row);
            rb_ary_store(mat_point, m, mat_point_row);
            rb_ary_store(mat_jump,  m, mat_jump_row);

            for (long n = 0; n < slen1+1; n++) {
              if ((m == 0) && (n == 0)) {
                // NONE
                rb_ary_store(mat_score_row, n, DBL2NUM(0.0));
                rb_ary_store(mat_point_row, n, INT2FIX(0));
              } else if ((m == 0) && (n == 1)) {
                // LEFT
                rb_ary_store(mat_score_row, n, DBL2NUM(-gopen));
                rb_ary_store(mat_point_row, n, INT2FIX(2));
              } else if ((m == 0) && (n > 1)) {
                // LEFT
                double prv = NUM2DBL(rb_ary_entry(rb_ary_entry(mat_score, m), n-1));
                rb_ary_store(mat_score_row, n, DBL2NUM(prv-gext));
                rb_ary_store(mat_point_row, n, INT2FIX(2));
              } else if ((m == 1) && (n == 0)) {
                // UP
                rb_ary_store(mat_score_row, n, DBL2NUM(-gopen));
                rb_ary_store(mat_point_row, n, INT2FIX(1));
              } else if ((m > 1) && (n == 0)) {
                // UP
                double prv = NUM2DBL(rb_ary_entry(rb_ary_entry(mat_score, m-1), n));
                rb_ary_store(mat_score_row, n, DBL2NUM(prv-gext));
                rb_ary_store(mat_point_row, n, INT2FIX(1));
              }
            }
          }

          // 3. initialize deletion and insertion matrices
          double n_infinity = -std::numeric_limits<double>::infinity();

          for (long m = 0; m < slen2+1; m++) {
            VALUE del_score_row = rb_ary_new2(slen1+1);
            VALUE del_point_row = rb_ary_new2(slen1+1);
            VALUE del_jump_row  = rb_ary_new2(slen1+1);

            VALUE ins_score_row = rb_ary_new2(slen1+1);
            VALUE ins_point_row = rb_ary_new2(slen1+1);
            VALUE ins_jump_row  = rb_ary_new2(slen1+1);

            rb_ary_store(del_score, m, del_score_row);
            rb_ary_store(del_point, m, del_point_row);
            rb_ary_store(del_jump,  m, del_jump_row);

            rb_ary_store(ins_score, m, ins_score_row);
            rb_ary_store(ins_point, m, ins_point_row);
            rb_ary_store(ins_jump,  m, ins_jump_row);

            // fill in  score matrices[0][0] with -infinity
            for (long n = 0; n < slen1+1; n++) {
              if ((m == 0) || (n == 0)) {
                rb_ary_store(del_score_row, n, DBL2NUM(n_infinity));
                rb_ary_store(ins_score_row, n, DBL2NUM(n_infinity));
              }
            }
          }

          // 4. fill in match, deletion, and insertion matrices
          for (long m = 1; m < slen2+1; m++) {
            VALUE aa2 = rb_ary_entry(aas2, m-1);
            for (long n = 1; n < slen1+1; n++) {
              VALUE aa1 = rb_ary_entry(aas1, n-1);

              double cur_mat_score = NUM2DBL(rb_funcall(subst, rb_intern("score"), 2, aa2, aa1));
              double prv_mat_score = NUM2DBL(rb_ary_entry(rb_ary_entry(mat_score, m-1), n-1));
              double prv_del_score = NUM2DBL(rb_ary_entry(rb_ary_entry(del_score, m-1), n-1));
              double prv_ins_score = NUM2DBL(rb_ary_entry(rb_ary_entry(ins_score, m-1), n-1));

              double mat_mat = prv_mat_score + cur_mat_score;
              double mat_del = prv_del_score + cur_mat_score;
              double mat_ins = prv_ins_score + cur_mat_score;
              double mat_max = 0;

              if (mat_mat >= mat_ins) {
                if (mat_mat >= mat_del) {
                  if (mat_mat > 0) {
                    // POINT DIAG
                    mat_max = mat_mat;
                    rb_ary_store(rb_ary_entry(mat_score, m), n, DBL2NUM(mat_max));
                    rb_ary_store(rb_ary_entry(mat_point, m), n, INT2FIX(3));
                  } else {
                    // NONE
                    rb_ary_store(rb_ary_entry(mat_score, m), n, DBL2NUM(0.0));
                    rb_ary_store(rb_ary_entry(mat_point, m), n, INT2FIX(0));
                  }
                } else {
                  if (mat_del > 0) {
                    // JUMP to DEL
                    mat_max = mat_del;
                    rb_ary_store(rb_ary_entry(mat_score, m), n, DBL2NUM(mat_max));
                    rb_ary_store(rb_ary_entry(mat_jump, m), n, rb_ary_new3(3, rb_str_new2("D"), LONG2NUM(m-1), LONG2NUM(n-1)));
                  } else {
                    // NONE
                    rb_ary_store(rb_ary_entry(mat_score, m), n, DBL2NUM(0.0));
                    rb_ary_store(rb_ary_entry(mat_point, m), n, INT2FIX(0));
                  }
                }
              } else {
                if (mat_ins > mat_del) {
                  if (mat_ins > 0) {
                    // JUMP to INS
                    mat_max = mat_ins;
                    rb_ary_store(rb_ary_entry(mat_score, m), n, DBL2NUM(mat_max));
                    rb_ary_store(rb_ary_entry(mat_jump, m), n, rb_ary_new3(3, rb_str_new2("I"), LONG2NUM(m-1), LONG2NUM(n-1)));
                  } else {
                    // NONE
                    rb_ary_store(rb_ary_entry(mat_score, m), n, DBL2NUM(0.0));
                    rb_ary_store(rb_ary_entry(mat_point, m), n, INT2FIX(0));
                  }
                } else {
                  if (mat_del > 0) {
                    // JUMP to DEL
                    mat_max = mat_del;
                    rb_ary_store(rb_ary_entry(mat_score, m), n, DBL2NUM(mat_max));
                    rb_ary_store(rb_ary_entry(mat_jump, m), n, rb_ary_new3(3, rb_str_new2("D"), LONG2NUM(m-1), LONG2NUM(n-1)));
                  } else {
                    // NONE
                    rb_ary_store(rb_ary_entry(mat_score, m), n, DBL2NUM(0.0));
                    rb_ary_store(rb_ary_entry(mat_point, m), n, INT2FIX(0));
                  }
                }
              }

              prv_mat_score   = NUM2DBL(rb_ary_entry(rb_ary_entry(mat_score, m), n-1));
              prv_del_score   = NUM2DBL(rb_ary_entry(rb_ary_entry(del_score, m), n-1));
              double del_mat  = prv_mat_score - gopen;
              double del_del  = prv_del_score - gext;;
              double del_max  = 0;

              if (del_mat >= del_del) {
                if (del_mat > 0) {
                  // JUMP to MAT
                  del_max = del_mat;
                  rb_ary_store(rb_ary_entry(del_score, m), n, DBL2NUM(del_max));
                  rb_ary_store(rb_ary_entry(del_jump, m), n, rb_ary_new3(3, rb_str_new2("M"), LONG2NUM(m), LONG2NUM(n-1)));
                } else {
                  // NONE
                  rb_ary_store(rb_ary_entry(del_score, m), n, DBL2NUM(0));
                  rb_ary_store(rb_ary_entry(del_point, m), n, INT2FIX(0));
                }
              } else {
                if (del_del > 0) {
                  // POINT LEFT
                  del_max = del_del;
                  rb_ary_store(rb_ary_entry(del_score, m), n, DBL2NUM(del_max));
                  rb_ary_store(rb_ary_entry(del_point, m), n, INT2FIX(2));
                } else {
                  // NONE
                  rb_ary_store(rb_ary_entry(del_score, m), n, DBL2NUM(0.0));
                  rb_ary_store(rb_ary_entry(del_point, m), n, INT2FIX(0));
                }
              }

              prv_mat_score   = NUM2DBL(rb_ary_entry(rb_ary_entry(mat_score, m-1), n));
              prv_ins_score   = NUM2DBL(rb_ary_entry(rb_ary_entry(ins_score, m-1), n));
              double ins_mat  = prv_mat_score - gopen;
              double ins_ins  = prv_ins_score - gext;;
              double ins_max  = 0;

              if (ins_mat >= ins_ins) {
                if (ins_mat > 0) {
                  // JUMP to MAT
                  ins_max = ins_mat;
                  rb_ary_store(rb_ary_entry(ins_score, m), n, DBL2NUM(ins_max));
                  rb_ary_store(rb_ary_entry(ins_jump, m), n, rb_ary_new3(3, rb_str_new2("M"), LONG2NUM(m-1), LONG2NUM(n)));
                } else {
                  // NONE
                  rb_ary_store(rb_ary_entry(ins_score, m), n, DBL2NUM(0.0));
                  rb_ary_store(rb_ary_entry(ins_point, m), n, INT2FIX(0));
                }
              } else {
                if (ins_ins > 0) {
                  // POINT UP
                  ins_max = ins_ins;
                  rb_ary_store(rb_ary_entry(ins_score, m), n, DBL2NUM(ins_max));
                  rb_ary_store(rb_ary_entry(ins_point, m), n, INT2FIX(1));
                } else {
                  // NONE
                  rb_ary_store(rb_ary_entry(ins_score, m), n, DBL2NUM(0.0));
                  rb_ary_store(rb_ary_entry(ins_point, m), n, INT2FIX(0));
                }
              }

              double max = 0;

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
          args[0]   = seq1;
          args[1]   = seq2;
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

          return rb_class_new_instance(13, args, rb_path2class("SequenceSequenceLocalAlignmentAffineGap"));
        }
      }
    end

    def local_alignment_affine_gap_rb(options={})
      opts = {
        :subst  => SubstitutionTable::blosum62,
        :gopen  => 10,
        :gext   => 0.5
      }.merge!(options)

      subst = opts[:subst]
      gopen = opts[:gopen]
      gext  = opts[:gext]
      aas1  = @sequence1.amino_acids
      slen1 = @sequence1.length
      aas2  = @sequence2.amino_acids
      slen2 = @sequence2.length
      max_s = 0
      max_m = nil
      max_n = nil

      # 1. create score, point, jump matrices for match, deletion, and insertion
      mat_score = Array.new(slen2+1)
      mat_point = Array.new(slen2+1)
      mat_jump  = Array.new(slen2+1)

      del_score = Array.new(slen2+1)
      del_point = Array.new(slen2+1)
      del_jump  = Array.new(slen2+1)

      ins_score = Array.new(slen2+1)
      ins_point = Array.new(slen2+1)
      ins_jump  = Array.new(slen2+1)

      # 2. initialize match matrices
      (0..slen2).each do |m|
        mat_score[m] = Array.new(slen1+1)
        mat_point[m] = Array.new(slen1+1)
        mat_jump[m]  = Array.new(slen1+1)

        (0..slen1).each do |n|
          if    (n == 0 && m == 0)
            mat_score[m][n] = 0
            mat_point[m][n] = NONE
          elsif (n == 1 && m == 0)
            mat_score[m][n] = -gopen
            mat_point[m][n] = LEFT
          elsif (n >  1 && m == 0)
            mat_score[m][n] = mat_score[m][n-1] - gext
            mat_point[m][n] = LEFT
          elsif (n == 0 && m == 1)
            mat_score[m][n] = -gopen
            mat_point[m][n] = UP
          elsif (n == 0 && m >  1)
            mat_score[m][n] = mat_score[m-1][n] - gext
            mat_point[m][n] = UP
          end
        end
      end

      # 3. initialize deletion and insertion matrices
      infinity = 1/0.0

      (0..slen2).each do |m|
        del_score[m] = Array.new(slen1+1)
        del_point[m] = Array.new(slen1+1)
        del_jump[m]  = Array.new(slen1+1)

        ins_score[m] = Array.new(slen1+1)
        ins_point[m] = Array.new(slen1+1)
        ins_jump[m]  = Array.new(slen1+1)

        (0..slen1).each do |n|
          del_score[m][n] = -infinity if (m == 0 || n == 0)
          ins_score[m][n] = -infinity if (m == 0 || n == 0)
        end
      end

      # 4. fill in match, deletion, and insertion matrices
      log_fmt = "%-12s : %5s : %12s"

      (1..slen2).each do |m|
        (1..slen1).each do |n|
          mat_mat = mat_score[m-1][n-1] + subst.score(aas1[n-1], aas2[m-1])
          mat_del = del_score[m-1][n-1] + subst.score(aas1[n-1], aas2[m-1])
          mat_ins = ins_score[m-1][n-1] + subst.score(aas1[n-1], aas2[m-1])
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

          del_mat = mat_score[m][n-1] - gopen
          del_del = del_score[m][n-1] - gext
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

          ins_mat = mat_score[m-1][n] - gopen
          ins_ins = ins_score[m-1][n] - gext
          ins_max = [0, ins_mat, ins_ins].max
          ins_score[m][n] = ins_max
          log_mat = "I[#{m}][#{n}]"

          case ins_max
          when 0
            ins_point[m][n] = NONE
            $logger.debug log_fmt % [log_mat, "POINT", "NONE"]
          when ins_mat
            ins_jump[m][n] = ['M', m-1, n]
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

      SequenceSequenceLocalAlignmentAffineGap.new(@sequence1, @sequence2,
                                                  mat_score, mat_point, mat_jump,
                                                  del_score, del_point, del_jump,
                                                  ins_score, ins_point, ins_jump,
                                                  max_m, max_n)
    end

    inline(:C) do |builder|
      builder.add_compile_flags '-x c++', '-lstdc++'
      builder.c_raw <<-EOCPP
        static VALUE global_alignment_linear_gap_cpp(int argc, VALUE *argv, VALUE self) {
          VALUE klass   = rb_const_get(rb_cObject, rb_intern("Gigue"));
          VALUE klass2  = rb_const_get(klass, rb_intern("SubstitutionTable"));
          VALUE subst   = rb_funcall(klass2, rb_intern("blosum62"), 0);

          VALUE options = argv[0] == Qnil ? rb_hash_new() : argv[0];
          VALUE opts    = rb_hash_new();

          rb_hash_aset(opts, ID2SYM(rb_intern("subst")), subst);
          rb_hash_aset(opts, ID2SYM(rb_intern("gap")), INT2NUM(8));
          rb_funcall(opts, rb_intern("merge!"), 1, options);

          long gap    = NUM2LONG(rb_hash_aref(opts, ID2SYM(rb_intern("gap"))));
          VALUE seq1  = rb_iv_get(self, "@sequence1");
          VALUE seq2  = rb_iv_get(self, "@sequence2");
          VALUE aas1  = rb_funcall(seq1, rb_intern("amino_acids"), 0);
          VALUE aas2  = rb_funcall(seq2, rb_intern("amino_acids"), 0);
          long slen1  = NUM2LONG(rb_funcall(seq1, rb_intern("length"), 0));
          long slen2  = NUM2LONG(rb_funcall(seq2, rb_intern("length"), 0));
          long max_s  = 0;
          long max_m  = 0;
          long max_n  = 0;

          // 1. create score and point matrices
          VALUE score = rb_ary_new2(slen2+1);
          VALUE point = rb_ary_new2(slen2+1);

          // 2. initialize score and point matrices
          for (long m = 0; m < slen2+1; m++) {
            VALUE score_row = rb_ary_new2(slen1+1);
            VALUE point_row = rb_ary_new2(slen1+1);

            rb_ary_store(score, m, score_row);
            rb_ary_store(point, m, point_row);

            for (long n = 0; n < slen1+1; n++) {
              if ((n == 0) && (m == 0)) {
                // NONE
                rb_ary_store(score_row, n, INT2NUM(0));
                rb_ary_store(point_row, n, INT2FIX(0));
              } else if ((n > 0) && (m == 0)) {
                // LEFT
                rb_ary_store(score_row, n, LONG2NUM(-gap*n));
                rb_ary_store(point_row, n, INT2FIX(2));
              } else if ((n == 0) && (m > 0)) {
                // UP
                rb_ary_store(score_row, n, LONG2NUM(-gap*m));
                rb_ary_store(point_row, n, INT2FIX(1));
              }
            }
          }

          // 3. fill in score and point matrices
          for (long m = 1; m < slen2+1; m++) {
            VALUE aa2 = rb_ary_entry(aas2, m-1);
            for (long n = 1; n < slen1+1; n++) {
              VALUE aa1 = rb_ary_entry(aas1, n-1);

              long cur_mat_score = NUM2LONG(rb_funcall(subst, rb_intern("score"), 2, aa1, aa2));
              long mat = NUM2LONG(rb_ary_entry(rb_ary_entry(score, m-1), n-1)) + cur_mat_score;
              long del = NUM2LONG(rb_ary_entry(rb_ary_entry(score, m), n-1)) - gap;
              long ins = NUM2LONG(rb_ary_entry(rb_ary_entry(score, m-1), n)) - gap;

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
          args[0] = seq1;
          args[1] = seq2;
          args[2] = score;
          args[3] = point;

          return rb_class_new_instance(4, args, rb_path2class("SequenceSequenceGlobalAlignmentLinearGap"));
        }
      EOCPP
    end

    def global_alignment_linear_gap_rb(options={})
      opts = {
        :subst  => SubstitutionTable::blosum62,
        :gap    => 8,
      }.merge!(options)

      subst = opts[:subst]
      gap   = opts[:gap]
      aas1  = @sequence1.amino_acids
      slen1 = @sequence1.length
      aas2  = @sequence2.amino_acids
      slen2 = @sequence2.length

      # 1. create score and point matrices
      score = Array.new(slen2+1)
      point = Array.new(slen2+1)

      # 2. initialize score and point matrices
      (0..slen2).each do |m|
        score[m] = Array.new(slen1+1)
        point[m] = Array.new(slen1+1)

        (0..slen1).each do |n|
          if    (n == 0 && m == 0)
            score[m][n] = 0
            point[m][n] = NONE
          elsif (n >  0 && m == 0)
            score[m][n] = -gap*n
            point[m][n] = LEFT
          elsif (n == 0 && m >  0)
            score[m][n] = -gap*m
            point[m][n] = UP
          end
        end
      end

      # 3. fill in score and point matrices
      (1..slen2).each do |m|
        (1..slen1).each do |n|
          mat = score[m-1][n-1] + subst.score(aas1[n-1], aas2[m-1])
          del = score[m][n-1] - gap
          ins = score[m-1][n] - gap
          max = [mat, del, ins].max

          score[m][n] = max
          point[m][n] = case max
                        when mat then DIAG
                        when del then LEFT
                        when ins then UP
                        end
        end
      end

      SequenceSequenceGlobalAlignmentLinearGap.new(@sequence1, @sequence2,
                                                  score, point)
    end

    inline(:C) do |builder|
      builder.add_compile_flags '-x c++', '-lstdc++'
      builder.include '<limits>'
      builder.c_raw %q{
        static VALUE global_alignment_affine_gap_cpp(int argc, VALUE *argv, VALUE self) {
          VALUE klass   = rb_const_get(rb_cObject, rb_intern("Gigue"));
          VALUE klass2  = rb_const_get(klass, rb_intern("SubstitutionTable"));
          VALUE subst   = rb_funcall(klass2, rb_intern("blosum62"), 0);
          double gopen  = 10.0;
          double gext   = 0.5;

          VALUE seq2  = rb_iv_get(self, "@sequence2");
          VALUE seq1  = rb_iv_get(self, "@sequence1");
          VALUE aas2  = rb_funcall(seq2, rb_intern("amino_acids"), 0);
          VALUE aas1  = rb_funcall(seq1, rb_intern("amino_acids"), 0);
          long slen1  = RARRAY_LEN(aas1);
          long slen2  = RARRAY_LEN(aas2);

          // 1. create score, point, jump matrices for match, deletion, and insertion
          VALUE mat_score = rb_ary_new2(slen2+1);
          VALUE mat_point = rb_ary_new2(slen2+1);
          VALUE mat_jump  = rb_ary_new2(slen2+1);

          VALUE del_score = rb_ary_new2(slen2+1);
          VALUE del_point = rb_ary_new2(slen2+1);
          VALUE del_jump  = rb_ary_new2(slen2+1);

          VALUE ins_score = rb_ary_new2(slen2+1);
          VALUE ins_point = rb_ary_new2(slen2+1);
          VALUE ins_jump  = rb_ary_new2(slen2+1);

          // 2. initialize match matrices

          for (long m = 0; m < slen2+1; m++) {
            VALUE mat_score_row = rb_ary_new2(slen1+1);
            VALUE mat_point_row = rb_ary_new2(slen1+1);
            VALUE mat_jump_row  = rb_ary_new2(slen1+1);

            rb_ary_store(mat_score, m, mat_score_row);
            rb_ary_store(mat_point, m, mat_point_row);
            rb_ary_store(mat_jump,  m, mat_jump_row);

            for (long n = 0; n < slen1+1; n++) {
              if ((m == 0) && (n == 0)) {
                // NONE
                rb_ary_store(mat_score_row, n, DBL2NUM(0.0));
                rb_ary_store(mat_point_row, n, INT2FIX(0));
              } else if ((m == 0) && (n == 1)) {
                // LEFT
                rb_ary_store(mat_score_row, n, DBL2NUM(-gopen));
                rb_ary_store(mat_point_row, n, INT2FIX(2));
              } else if ((m == 0) && (n > 1)) {
                // LEFT
                double prv = NUM2DBL(rb_ary_entry(rb_ary_entry(mat_score, m), n-1));
                rb_ary_store(mat_score_row, n, DBL2NUM(prv-gext));
                rb_ary_store(mat_point_row, n, INT2FIX(2));
              } else if ((m == 1) && (n == 0)) {
                // UP
                rb_ary_store(mat_score_row, n, DBL2NUM(-gopen));
                rb_ary_store(mat_point_row, n, INT2FIX(1));
              } else if ((m > 1) && (n == 0)) {
                // UP
                double prv = NUM2DBL(rb_ary_entry(rb_ary_entry(mat_score, m-1), n));
                rb_ary_store(mat_score_row, n, DBL2NUM(prv-gext));
                rb_ary_store(mat_point_row, n, INT2FIX(1));
              }
            }
          }

          // 3. initialize deletion and insertion matrices
          double n_infinity = -std::numeric_limits<double>::infinity();

          for (long m = 0; m < slen2+1; m++) {
            VALUE del_score_row = rb_ary_new2(slen1+1);
            VALUE del_point_row = rb_ary_new2(slen1+1);
            VALUE del_jump_row  = rb_ary_new2(slen1+1);

            VALUE ins_score_row = rb_ary_new2(slen1+1);
            VALUE ins_point_row = rb_ary_new2(slen1+1);
            VALUE ins_jump_row  = rb_ary_new2(slen1+1);

            rb_ary_store(del_score, m, del_score_row);
            rb_ary_store(del_point, m, del_point_row);
            rb_ary_store(del_jump,  m, del_jump_row);

            rb_ary_store(ins_score, m, ins_score_row);
            rb_ary_store(ins_point, m, ins_point_row);
            rb_ary_store(ins_jump,  m, ins_jump_row);

            // fill in  score matrices[0][0] with -infinity
            for (long n = 0; n < slen1+1; n++) {
              if ((m == 0) || (n == 0)) {
                rb_ary_store(del_score_row, n, DBL2NUM(n_infinity));
                rb_ary_store(ins_score_row, n, DBL2NUM(n_infinity));
              }
            }
          }

          // 4. fill in match, deletion, and insertion matrices
          for (long m = 1; m < slen2+1; m++) {
            VALUE aa2 = rb_ary_entry(aas2, m-1);
            for (long n = 1; n < slen1+1; n++) {
              VALUE aa1 = rb_ary_entry(aas1, n-1);

              double cur_mat_score = NUM2DBL(rb_funcall(subst, rb_intern("score"), 2, aa2, aa1));
              double prv_mat_score = NUM2DBL(rb_ary_entry(rb_ary_entry(mat_score, m-1), n-1));
              double prv_del_score = NUM2DBL(rb_ary_entry(rb_ary_entry(del_score, m-1), n-1));
              double prv_ins_score = NUM2DBL(rb_ary_entry(rb_ary_entry(ins_score, m-1), n-1));

              double mat_mat = prv_mat_score + cur_mat_score;
              double mat_del = prv_del_score + cur_mat_score;
              double mat_ins = prv_ins_score + cur_mat_score;

              if (mat_mat >= mat_ins) {
                if (mat_mat >= mat_del) {
                  // POINT DIAG
                  rb_ary_store(rb_ary_entry(mat_score, m), n, DBL2NUM(mat_mat));
                  rb_ary_store(rb_ary_entry(mat_point, m), n, INT2FIX(3));
                } else {
                  // JUMP to DEL
                  rb_ary_store(rb_ary_entry(mat_score, m), n, DBL2NUM(mat_del));
                  rb_ary_store(rb_ary_entry(mat_jump, m), n, rb_ary_new3(3, rb_str_new2("D"), LONG2NUM(m-1), LONG2NUM(n-1)));
                }
              } else {
                if (mat_ins > mat_del) {
                  // JUMP to INS
                  rb_ary_store(rb_ary_entry(mat_score, m), n, DBL2NUM(mat_ins));
                  rb_ary_store(rb_ary_entry(mat_jump, m), n, rb_ary_new3(3, rb_str_new2("I"), LONG2NUM(m-1), LONG2NUM(n-1)));
                } else {
                  // JUMP to DEL
                  rb_ary_store(rb_ary_entry(mat_score, m), n, DBL2NUM(mat_del));
                  rb_ary_store(rb_ary_entry(mat_jump, m), n, rb_ary_new3(3, rb_str_new2("D"), LONG2NUM(m-1), LONG2NUM(n-1)));
                }
              }

              prv_mat_score   = NUM2DBL(rb_ary_entry(rb_ary_entry(mat_score, m), n-1));
              prv_del_score   = NUM2DBL(rb_ary_entry(rb_ary_entry(del_score, m), n-1));
              double del_mat  = prv_mat_score - gopen;
              double del_del  = prv_del_score - gext;

              if (del_mat >= del_del) {
                // JUMP to MAT
                rb_ary_store(rb_ary_entry(del_score, m), n, DBL2NUM(del_mat));
                rb_ary_store(rb_ary_entry(del_jump, m), n, rb_ary_new3(3, rb_str_new2("M"), LONG2NUM(m), LONG2NUM(n-1)));
              } else {
                // POINT LEFT
                rb_ary_store(rb_ary_entry(del_score, m), n, DBL2NUM(del_del));
                rb_ary_store(rb_ary_entry(del_point, m), n, INT2FIX(2));
              }

              prv_mat_score   = NUM2DBL(rb_ary_entry(rb_ary_entry(mat_score, m-1), n));
              prv_ins_score   = NUM2DBL(rb_ary_entry(rb_ary_entry(ins_score, m-1), n));
              double ins_mat  = prv_mat_score - gopen;
              double ins_ins  = prv_ins_score - gext;

              if (ins_mat >= ins_ins) {
                // JUMP to MAT
                rb_ary_store(rb_ary_entry(ins_score, m), n, DBL2NUM(ins_mat));
                rb_ary_store(rb_ary_entry(ins_jump, m), n, rb_ary_new3(3, rb_str_new2("M"), LONG2NUM(m-1), LONG2NUM(n)));
              } else {
                // POINT UP
                rb_ary_store(rb_ary_entry(ins_score, m), n, DBL2NUM(ins_ins));
                rb_ary_store(rb_ary_entry(ins_point, m), n, INT2FIX(1));
              }
            }
          }

          VALUE args[11];
          args[0]   = seq1;
          args[1]   = seq2;
          args[2]   = mat_score;
          args[3]   = mat_point;
          args[4]   = mat_jump;
          args[5]   = del_score;
          args[6]   = del_point;
          args[7]   = del_jump;
          args[8]   = ins_score;
          args[9]   = ins_point;
          args[10]  = ins_jump;

          return rb_class_new_instance(11, args, rb_path2class("SequenceSequenceGlobalAlignmentAffineGap"));
        }
      }
    end

    def global_alignment_affine_gap_rb
      subst = SubstitutionTable::blosum62
      gopen = 10.0
      gext  = 0.5
      aas1  = @sequence1.amino_acids
      slen1 = @sequence1.length
      aas2  = @sequence2.amino_acids
      slen2 = @sequence2.length

      # 1. create score, point, jump matrices for match, deletion, and insertion
      mat_score = Array.new(slen2+1)
      mat_point = Array.new(slen2+1)
      mat_jump  = Array.new(slen2+1)

      del_score = Array.new(slen2+1)
      del_point = Array.new(slen2+1)
      del_jump  = Array.new(slen2+1)

      ins_score = Array.new(slen2+1)
      ins_point = Array.new(slen2+1)
      ins_jump  = Array.new(slen2+1)

      # 2. initialize match matrices
      (0..slen2).each do |m|
        mat_score[m] = Array.new(slen1+1)
        mat_point[m] = Array.new(slen1+1)
        mat_jump[m]  = Array.new(slen1+1)

        (0..slen1).each do |n|
          if    (n == 0 && m == 0)
            mat_score[m][n] = 0
            mat_point[m][n] = NONE
          elsif (n == 1 && m == 0)
            mat_score[m][n] = -gopen
            mat_point[m][n] = LEFT
          elsif (n >  1 && m == 0)
            mat_score[m][n] = mat_score[m][n-1] - gext
            mat_point[m][n] = LEFT
          elsif (n == 0 && m == 1)
            mat_score[m][n] = -gopen
            mat_point[m][n] = UP
          elsif (n == 0 && m >  1)
            mat_score[m][n] = mat_score[m-1][n] - gext
            mat_point[m][n] = UP
          end
        end
      end

      # 3. initialize deletion and insertion matrices
      infinity = 1/0.0

      (0..slen2).each do |m|
        del_score[m] = Array.new(slen1+1)
        del_point[m] = Array.new(slen1+1)
        del_jump[m]  = Array.new(slen1+1)

        ins_score[m] = Array.new(slen1+1)
        ins_point[m] = Array.new(slen1+1)
        ins_jump[m]  = Array.new(slen1+1)

        (0..slen1).each do |n|
          del_score[m][n] = -infinity if (m == 0 || n == 0)
          ins_score[m][n] = -infinity if (m == 0 || n == 0)
        end
      end

      # 4. fill in match, deletion, and insertion matrices
      log_fmt = "%-12s : %5s : %12s"

      (1..slen2).each do |m|
        (1..slen1).each do |n|
          mat_mat = mat_score[m-1][n-1] + subst.score(aas1[n-1], aas2[m-1])
          mat_del = del_score[m-1][n-1] + subst.score(aas1[n-1], aas2[m-1])
          mat_ins = ins_score[m-1][n-1] + subst.score(aas1[n-1], aas2[m-1])
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

          del_mat = mat_score[m][n-1] - gopen
          del_del = del_score[m][n-1] - gext
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

          ins_mat = mat_score[m-1][n] - gopen
          ins_ins = ins_score[m-1][n] - gext
          ins_max = [ins_mat, ins_ins].max
          ins_score[m][n] = ins_max
          log_mat = "I[#{m}][#{n}]"

          case ins_max
          when ins_mat
            ins_jump[m][n] = ['M', m-1, n]
            $logger.debug log_fmt % [log_mat, "JUMP", "M[#{m-1}][#{n}]"]
          when ins_ins
            ins_point[m][n] = UP
            $logger.debug log_fmt % [log_mat, "POINT", "I[#{m-1}][#{n}]"]
          end
        end
      end

      SequenceSequenceGlobalAlignmentAffineGap.new(@sequence1, @sequence2,
                                                   mat_score, mat_point, mat_jump,
                                                   del_score, del_point, del_jump,
                                                   ins_score, ins_point, ins_jump)
    end

  end
end
