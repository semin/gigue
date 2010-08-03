module ArrayExtensions

  unless method_defined?(:to_hash)
    def to_hash(other)
      #Hash[ *(0...self.size()).inject([]) { |arr, ix|
        #arr.push(self[ix], other[ix])
      #} ]
      Hash[*self.zip(other).flatten]
    end
  end

  unless method_defined?(:sum)
    def sum(identity = 0, &block)
      if block_given?
        map(&block).sum
      else
        inject{ |sum, element| sum + element } || identity
      end
    end
  end

end

Array.send :include, ArrayExtensions
