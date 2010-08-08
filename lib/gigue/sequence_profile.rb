module Gigue
  class SequenceProfile

    attr_reader :msa, :positions

    def initialize(msa)
      @msa        = msa
      @positions  = @msa.columns.map { |c| SequenceProfilePosition.new(c.probe) }
    end

    def sequences
      @msa.sequences
    end

    def length
      @msa.length
    end

  end
end
