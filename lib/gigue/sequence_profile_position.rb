module Gigue
  class SequenceProfilePosition

    attr_reader :probe

    def initialize(probe, aa_raw_frqs=nil, aa_frqs=nil, aa_prbs=nil, aa_rel_prbs=nil)
      @probe        = probe
      @aa_raw_frqs  = aa_raw_frqs
      @aa_frqs      = aa_frqs
      @aa_prbs      = aa_prbs
      @aa_rel_prbs  = aa_rel_prbs
    end

    def raw_frequency_of(aa)
      @aa_raw_frqs[aa]
    end

    # weighting applied
    def frequency_of(aa)
      @aa_frqs[aa]
    end

    # no weighting applied
    def raw_probability_of(aa)
      sum = Float(@aa_raw_frqs.values.sum)
      @aa_raw_frqs[aa] / sum
    end

    # no weighting applied and no gap considered
    def raw_relative_probability_of(aa)
      gap = @aa_raw_frqs['-']
      sum = Float(@aa_raw_frqs.values.sum - gap)
      @aa_raw_frqs[aa] / sum
    end

    # weighting applied and gap considered
    def probability_of(aa)
      @aa_prbs[aa]
    end

    # weighting applied and no gap considered
    def relative_probability_of(aa)
      @aa_rel_prbs[aa]
    end

  end
end
