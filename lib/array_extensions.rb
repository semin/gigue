module ArrayExtensions

  unless method_defined?(:to_hash)
    def to_hash(other)
      Hash[ *(0...self.size()).inject([]) { |arr, ix|
        arr.push(self[ix], other[ix])
      } ]
    end
  end

end

Array.send :include, ArrayExtensions
