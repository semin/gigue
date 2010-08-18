module Gigue
  class Sequence

    attr_reader :amino_acids, :code, :description

    def initialize(data, code=nil, desc=nil)
      @code         = code
      @description  = desc
      @amino_acids  = data.split('')
    end

    def data
      @amino_acids.join('')
    end

    def length
      @amino_acids.length
    end

    def [](index)
      @amino_acids[index]
    end

    def []=(index, value)
      @amino_acids[index] = value
    end

    def gapless
      gapless_data = @amino_acids.reject { |a| a == '-' }.join('')
      self.class.new(gapless_data, @code, @description)
    end

    def gapless!
      @amino_acids.delete('-')
      self
    end

    def shuffle
      self.class.new(@amino_acids.shuffle.join(''), @code, @description)
    end

    def shuffle!
      @amino_acids.shuffle!
      self
    end

    inline do |builder|
      builder.include "<stdlib.h>"
      builder.c_raw <<-EOC
        static VALUE shuffle_c(int argc, VALUE *argv, VALUE self) {
          VALUE aas1 = rb_iv_get(self, "@amino_acids");
          VALUE code = rb_iv_get(self, "@code");
          VALUE desc = rb_iv_get(self, "@description");
          long len1 = RARRAY_LEN(aas1);
          VALUE aas2 = rb_ary_new2(len1);
          long i, r;

          for(i = 0; i < len1; i++) {
            r = rand() % (len1 - i) + i;
            rb_ary_store(aas2, r, rb_ary_entry(aas1, i));
            rb_ary_store(aas2, i, rb_ary_entry(aas1, r));
          }

          VALUE data = rb_funcall(aas2, rb_intern("join"), 1, rb_str_new2(""));
          VALUE args[3];

          args[0] = data;
          args[1] = code;
          args[2] = desc;

          return rb_class_new_instance(3, args, rb_path2class("Sequence"));
        }
      EOC
    end

  end
end
