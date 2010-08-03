module Gigue
  class SequenceProfile

    attr_reader :msa, :length, :positions

    def initialize(msa)
      @msa        = msa
      @length     = @msa.length
      @positions  = @msa.columns.map { |c| SequenceProfilePosition.new(c.probe) }
    end

  end
end
