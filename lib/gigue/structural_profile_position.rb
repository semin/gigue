module Gigue
  class StructuralProfilePosition

    attr_reader :probe

    def initialize(probe, mat_score=nil, gap_score=nil)
      @probe      = probe
      @mat_score  = mat_score
      @gap_score  = gap_score
    end

    def mat_score(aa)
      @mat_score[aa]
    end

    def gap_score(scheme)
      @gap_score[scheme]
    end

    def gap_ins_open
      @gap_score['InsO']
    end

    def gap_ins_ext
      @gap_score['InsE']
    end

    def gap_del_open
      @gap_score['DelO']
    end

    def gap_del_ext
      @gap_score['DelE']
    end

  end
end
