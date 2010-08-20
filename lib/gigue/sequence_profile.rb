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

    def shuffle
      self.class.new(@msa.shuffle)
    end

    def shuffle!
      @msa = @msa.shuffle
      @positions = @msa.columns.map { |c| SequenceProfilePosition.new(c.probe) }
    end

  end
end
