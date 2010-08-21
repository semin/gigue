module ArrayExtensions

  unless method_defined?(:to_hash)
    def to_hash(other)
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
