module Gigue

  class ProfilePosition

    attr_reader :probe

    def initialize(probe, mat_score, gap_score)
      @probe      = probe
      @mat_score  = mat_score
      @gap_score  = gap_score
    end

  end

  class StructuralProfilePosition < ProfilePosition

    def initialize(probe, mat_score, gap_score)
      super(probe, mat_score, gap_score)
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

  class FugueProfilePosition < ProfilePosition

    def initialize(probe, mat_score, gap_score, env_score)
      super(probe, mat_score, gap_score)
      @env_score = env_score
    end

    def mat_score(aa)
      @mat_score[aa]
    end

    def gap_score(scheme)
      @gap_score[scheme]
    end

    def env_score(scheme)
      @env_score[scheme]
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