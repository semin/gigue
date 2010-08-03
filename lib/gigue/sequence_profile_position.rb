module Gigue
  class SequenceProfilePosition

    attr_reader :probe

    def initialize(probe)
      @probe    = probe
      @aa_freqs = Hash.new(0)
      @aa_probs = Hash.new(0.0)
      @rel_aa_probs = Hash.new(0.0)
      @probe.split('').each do |c|
        if AMINO_ACIDS.include?(c) || c == '-'
          @aa_freqs[c] += 1
        else
          $logger.warn "#{aa} is a unknown type of amino acid and ignored."
        end
      end
      sum     = @aa_freqs.values.sum.to_f
      aa_sum  = sum - @aa_freqs['-']
      @aa_freqs.each do |a, f|
        @aa_probs[a]      = f / sum
        @rel_aa_probs[a]  = f / aa_sum
      end
    end

    def frequency_of(aa)
      @aa_freqs[aa]
    end

    def probability_of(aa)
      @aa_probs[aa]
    end

    def relative_probability_of(aa)
      @rel_aa_probs[aa]
    end

  end
end
