module Gigue
  class Sequence

    def self.calculate_equal_weights(seqs)
      seqs.each { |s| s.weight = 1.0 / seqs.size }
    end

    def self.calculate_va_weights(seqs)
      begin
        self.calculate_va_weights_cpp(seqs)
      rescue
        self.calculate_va_weights_rb(seqs)
      end
    end

    def self.calculate_va_weights(seqs)
      tot   = 0.0
      dists = Array.new(seqs.size)

      # calculate all by all percentage dissimilarity
      (0...(seqs.size-1)).each do |i|
        dists[i] = Array.new(seqs.size)
        ((i+1)...seqs.size).each do |j|
          d = 1 - seqs[i].pid(seqs[j])
          tot += d
          dists[i][j] = d
        end
      end

      # calculate VA weights
      (0...seqs.size).each do |i|
        sum = 0.0
        (0...seqs.size).each do |j|
          if (i < j)
            sum += dists[i][j]
          elsif (i > j)
            sum += dists[j][i]
          end
        end
        w = (sum / 2.0) / tot
        seqs[i].weight = w
      end
    end

    inline(:C) do |builder|
      builder.add_compile_flags '-x c++', '-lstdc++'
      builder.c_singleton <<-EOCPP
        static VALUE calculate_va_weights_cpp(VALUE seqs) {
          VALUE *seqs_ptr = RARRAY_PTR(seqs);
          long seq_cnt    = RARRAY_LEN(seqs);
          double tot      = 0.0;
          double dists[seq_cnt][seq_cnt];

          for (long i = 0; i < seq_cnt-1; i++) {
            for (long j = i+1; j < seq_cnt; j++) {
              double d = 1.0 - NUM2DBL(rb_funcall(seqs_ptr[i], rb_intern("pid"), 1, seqs_ptr[j])); 
              tot += d;
              dists[i][j] = d;
            }
          }

          for (long i = 0; i < seq_cnt; i++) {
            double sum = 0.0;
            for (long j = 0; j < seq_cnt; j++) {
              if (i < j) {
                sum += dists[i][j];
              } else if (i > j) {
                sum += dists[j][i];
              }
            }
            double w = (sum / 2.0) / tot;
            rb_funcall(seqs_ptr[i], rb_intern("weight="), 1, DBL2NUM(w));
          }
          return Qnil;
        }
      EOCPP
    end

    def self.calculate_blosum_weights(seqs, weight=0.6)
      begin
        self.calculate_blosum_weights_cpp(seqs, weight)
      rescue
        self.calculate_blosum_weights_rb(seqs, weight)
      end
    end

    def self.calculate_blosum_weights_rb(seqs, weight=0.6)
      clusters = seqs.map { |s| [s] }
      begin
        continue = false
        0.upto(clusters.size-2) do |i|
          indexes = []
          (i+1).upto(clusters.size-1) do |j|
            found = false
            clusters[i].each do |s1|
              clusters[j].each do |s2|
                if s1.pid(s2) >= weight
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

      seq_cnt = Float(seqs.size)

      clusters.each do |cluster|
        cluster.each do |seq|
          weight = cluster.size / seq_cnt
          seq.weight = weight
        end
      end
    end


    attr_accessor :data, :code, :description,
                  :weight, :environments, :gap_cnt,
                  :ungapped_data, :ungapped_length, :ungapped_environments


    def initialize(data, code=nil, desc=nil)
      @data         = data.gsub('X','').gsub('Z','Q')
      @code         = code
      @description  = desc
    end

    def amino_acids
      @data.split('')
    end

    def length
      @data.length
    end

    def [](index)
      @data[index]
    end

    def []=(index, value)
      @data[index] = value
    end

    def gapless
      self.class.new(@data.gsub('-', ''), @code, @description)
    end

    def gapless!
      @data.gsub!('-', '')
      self
    end

    def reverse
      self.class.new(@data.reverse, @code, @description)
    end

    def reverse!
      @data.reverse!
      self
    end

    def shuffle
      self.class.new(@data.split('').shuffle.join(''), @code, @description)
    end

    def shuffle!
      @data = @data.split('').shuffle.join('')
      self
    end

    def pid(other)
      begin
        pid_cpp(other)
      rescue
        pid_rb(other)
      end
    end

    def pid_rb(other)
      aas1  = amino_acids
      aas2  = other.amino_acids
      cols  = aas1.zip(aas2)
      gap   = '-'
      align = 0.0 # no. of aligned columns
      ident = 0.0 # no. of identical columns
      intgp = 0.0 # no. of internal gaps
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
      Float(ident) / (align + intgp)
    end

    inline(:C) do |builder|
      builder.add_compile_flags '-x c++', '-lstdc++'
      builder.c <<-EOCPP
        static VALUE pid_cpp(VALUE other) {
          VALUE seq1 = rb_iv_get(self, "@data");
          VALUE seq2 = rb_iv_get(other, "@data");
          char *seq1_ptr = RSTRING_PTR(seq1);
          char *seq2_ptr = RSTRING_PTR(seq2);
          double align = 0.0;
          double ident = 0.0;
          double intgp = 0.0;

          for (long i = 0; i < RSTRING_LEN(seq1); i++) {
            if ((seq1_ptr[i] != '-') && (seq2_ptr[i] != '-')) {
              align += 1;
              if (seq1_ptr[i] == seq2_ptr[i]) {
                ident += 1;
              }
            } else if (((seq1_ptr[i] == '-') && (seq2_ptr[i] != '-')) ||
                       ((seq1_ptr[i] != '-') && (seq2_ptr[i] == '-'))) {
              intgp += 1;
            }
          }
          return DBL2NUM(ident / (align + intgp));
        }
      EOCPP
    end

    def to_flatfile(options={})
      opts = {
        :os     => STDOUT,
        :type   => :pir,
        :width  => 70
      }.merge!(options)

      out = opts[:os].is_a?(String) ? File.open(opts[:os], 'w') : opts[:os]

      out.puts opts[:type] == :pir ? ">P1;#{@code}" : ">#{@code}"
      out.puts "sequence" if opts[:type] == :pir

      aas = amino_acids << '*'
      out.puts aas.map_with_index { |a, ai|
        a + (ai > 0 && (ai+1) % opts[:width] == 0 ? "\n" : '')
      }.join('').chomp

      out.close if [File, String].include?(out.class)
    end
  end
end
