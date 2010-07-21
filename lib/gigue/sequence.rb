module Gigue
  class Sequence

    attr_reader :code, :data

    def initialize(code, data)
      @code = code
      @data = data
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

  end
end
