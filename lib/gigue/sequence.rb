module Gigue
  class Sequence

    attr_reader :code, :data, :description

    def initialize(data, code=nil, desc=nil)
      @code = code
      @data = data
      @description = desc
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

    def gapless
      Sequence.new(@data.gsub('-', ''), @code, @description)
    end

    def gapless!
      @data.gsub!('-', '')
      self
    end

  end
end
