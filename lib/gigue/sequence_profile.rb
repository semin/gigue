module Gigue
  class SequenceProfile

    attr_reader :msa, :positions

    def initialize(msa, pss)
      @msa = msa
      @positions = pss
    end

    def sequences
      @msa.sequences
      #seqs = Array.new(depth, '')
      #(0...length).each do |pi|
        #(0...depth).each do |si|
          #seqs[si] += @positions[pi].probe[si]
        #end
      #end
      #seqs
    end

    # no. of sequences
    def depth
      @msa.depth
      #@positions[0].probe.size
    end

    def length
      @msa.length
      #@positions.length
    end

    def shuffle
      self.class.new(@msa, @positions.shuffle)
    end

    def shuffle!
      @positions.shuffle!
    end

  end
end
