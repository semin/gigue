module Gigue
  class ProfilePosition

    attr_reader :probe, :_mat_score, :_gap_score

    def initialize(probe, mat_score, gap_score)
      @probe      = probe
      @_mat_score  = mat_score
      @_gap_score  = gap_score
    end

  end

  class FugueProfilePosition < ProfilePosition

    attr_reader :_env_score

    def initialize(probe, mat_score, gap_score, env_score)
      super(probe, mat_score, gap_score)
      @_env_score = env_score
    end

    def mat_score(aa)
      @_mat_score[aa]
    end

    def gap_ins_open
      -@_gap_score['InsO']
    end

    def gap_ins_ext
      -@_gap_score['InsE']
    end

    def gap_del_open
      -@_gap_score['DelO']
    end

    def gap_del_ext
      -@_gap_score['DelE']
    end

  end
end
