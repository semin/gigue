module Gigue
  class Esst

    attr_accessor :type, :name, :no, :colnames, :rownames, :matrix

    def initialize(type, name, no, colnames = [], rownames = [], matrix = nil)
      @type     = type
      @name     = name
      @no       = no
      @colnames = colnames
      @rownames = rownames
      @matrix   = matrix
    end

    def scores_from(aa)
      i = colnames.index(aa)
      @matrix[i, 0..-1]
    end

    def scores_to(aa)
      j = rownames.index(aa)
      @matrix[0..-1, j]
    end

  end
end
