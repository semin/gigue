module Gigue
  class ProfileProfileAligner

    attr_reader :structural_profile, :sequence_profile

    def initialize(str_prf, seq_prf)
      @structural_profile = str_prf
      @sequence_profile   = seq_prf
    end

    def local_alignment_with_linear_gap_penalty(gap_del=100, gap_ins=100)
      str_pss   = @structural_profile.positions
      seq_pss   = @sequence_profile.positions
      str_plen  = @structural_profile.length
      seq_plen  = @sequence_profile.length
      max_s     = 0
      max_i     = nil
      max_j     = nil
      stmat     = NArray.object(str_plen+1, seq_plen+1)

      # initialize score and trace matrix
      (0..str_plen).each do |str_i|
        (0..seq_plen).each do |seq_i|
          stmat[str_i, seq_i] = { :score => 0, :point => nil }
          if    (str_i == 0 && seq_i == 0)
            stmat[str_i, seq_i][:point] = NONE
          elsif (str_i >  0 && seq_i == 0)
            stmat[str_i, seq_i][:point] = LEFT
          elsif (str_i == 0 && seq_i >  0)
            stmat[str_i, seq_i][:point] = UP
          end
        end
      end

      # fill in score and trace matrix
      (1..str_plen).each do |n|
        (1..seq_plen).each do |m|
          # calculate profile-profile match score
          mat_score = 0
          AMINO_ACIDS.split('').each do |aa|
            #w = seq_pss[seq_i].probability_of(aa)
            w = seq_pss[m-1].relative_probability_of(aa)
            s = str_pss[n-1].mat_score(aa == 'C' ? 'U' : aa)
            mat_score += (w * s)
          end
          mat_score = mat_score.round

          mat = stmat[n-1,  m-1][:score] + mat_score
          del = stmat[n-1,    m][:score] - gap_del
          ins = stmat[n,    m-1][:score] - gap_ins
          max = [0, mat, del, ins].max

          stmat[n, m][:score] = max
          stmat[n, m][:point] = case max
                                when 0    then NONE
                                when mat  then DIAG
                                when del  then LEFT
                                when ins  then UP
                                end
          if max >= max_s
            max_i = n
            max_j = m
            max_s = max
          end
        end
      end

      ProfileProfileLocalAlignmentLinearGap.new(@structural_profile, @sequence_profile,
                                                stmat, max_i, max_j)
    end

    def local_alignment_with_affine_gap_penalty
      str_pss   = @structural_profile.positions
      seq_pss   = @sequence_profile.positions
      str_plen  = @structural_profile.length
      seq_plen  = @sequence_profile.length
      max_s     = 0
      max_i     = nil
      max_j     = nil
      max_m     = nil

      # create score, deletion, and insertion matrices
      mat = NArray.object(str_plen+1, seq_plen+1)
      del = NArray.object(str_plen+1, seq_plen+1)
      ins = NArray.object(str_plen+1, seq_plen+1)

      # initialize score matrix and fill in the first row and column
      prev_gap_ins_ext = nil
      (0..str_plen).each do |str_i|
        (0..seq_plen).each do |seq_i|
          mat[str_i, seq_i] = { :score => 0, :point => nil, :jump => nil }
          if    (str_i == 0 && seq_i == 0)
            mat[str_i, seq_i][:point] = NONE
          elsif (str_i == 1 && seq_i == 0)
            mat[str_i, seq_i][:score] = -str_pss[str_i-1].gap_del_open
            mat[str_i, seq_i][:point] = LEFT
          elsif (str_i >  1 && seq_i == 0)
            mat[str_i, seq_i][:score] = mat[str_i-1, seq_i][:score] - str_pss[str_i-1].gap_del_ext
            mat[str_i, seq_i][:point] = LEFT
          elsif (str_i == 0 && seq_i == 1)
            mat[str_i, seq_i][:score] = -str_pss[str_i].gap_ins_open
            mat[str_i, seq_i][:point] = UP
            prev_gap_ins_ext    = str_pss[str_i].gap_ins_ext
          elsif (str_i == 0 && seq_i >  1)
            mat[str_i, seq_i][:score] = mat[str_i, seq_i-1][:score] - prev_gap_ins_ext
            mat[str_i, seq_i][:point] = UP
          end
        end
      end

      # initialize deletion and insertion matrices
      infinity = 1/0.0
      (0..str_plen).each do |str_i|
        (0..seq_plen).each do |seq_i|
          del[str_i, seq_i] = { :score => (str_i == 0 || seq_i == 0 ? -infinity : 0), :point => nil, :jump => nil }
          ins[str_i, seq_i] = { :score => (str_i == 0 || seq_i == 0 ? -infinity : 0), :point => nil, :jump => nil }
        end
      end

      # fill in match, deletion, and insertion matrices
      prev_gap_ins_ext = str_pss[0].gap_ins_ext
      (1..str_plen).each do |n|
        (1..seq_plen).each do |m|
          str_i, seq_i = n-1, m-1

          # calculate profile-profile match score
          mat_score = 0
          AMINO_ACIDS.split('').each do |aa|
            #w = seq_pss[seq_i].probability_of(aa)
            w = seq_pss[seq_i].relative_probability_of(aa)
            s = str_pss[str_i].mat_score(aa == 'C' ? 'U' : aa)
            mat_score += (w * s)
          end
          mat_score = mat_score.round

          # update match matrix
          mat_mat = mat[n-1, m-1][:score] + mat_score
          mat_del = del[n-1, m-1][:score] + mat_score
          mat_ins = ins[n-1, m-1][:score] + mat_score
          mat_max = [0, mat_mat, mat_del, mat_ins].max
          mat[n, m][:score] = mat_max

          case mat_max
          when 0
            mat[n, m][:point] = NONE
            $logger.debug "%-10s %5s" % ["M-#{n}-#{m}", "NONE"]
          when mat_mat
            mat[n, m][:point] = DIAG
            $logger.debug "%-10s %5s %10s" % ["M-#{n}-#{m}", "POINT", "M-#{n-1}-#{m-1}"]
          when mat_del
            jmp = "D-#{n-1}-#{m-1}"
            mat[n, m][:jump] = jmp
            $logger.debug "%-10s %5s %10s" % ["M-#{n}-#{m}", "JUMP", jmp]
          when mat_ins
            jmp = "I-#{n-1}-#{m-1}"
            mat[n, m][:jump] = jmp
            $logger.debug "%-10s %5s %10s" % ["M-#{n}-#{m}", "JUMP", jmp]
          end

          # update deletion matrix
          del_mat = mat[n-1, m][:score] - str_pss[str_i].gap_del_open
          del_del = del[n-1, m][:score] - str_pss[str_i].gap_del_ext
          del_max = [0, del_mat, del_del].max
          del[n, m][:score] = del_max

          case del_max
          when 0
            del[n, m][:point] = NONE
            $logger.debug "%-10s %5s" % ["D-#{n}-#{m}", "NONE"]
          when del_mat
            jmp = "M-#{n-1}-#{m}"
            del[n, m][:jump] = jmp
            $logger.debug "%-10s %5s %10s" % ["D-#{n}-#{m}", "JUMP", jmp]
          when del_del
            del[n, m][:point] = LEFT
            $logger.debug "%-10s %5s %10s" % ["D-#{n}-#{m}", "POINT", "D-#{n-1}-#{m}"]
          end

          # update insertion matrix
          ins_mat = mat[n, m-1][:score] - str_pss[str_i].gap_ins_open
          ins_ins = ins[n, m-1][:score] - prev_gap_ins_ext
          ins_max = [0, ins_mat, ins_ins].max
          ins[n, m][:score] = ins_max

          case ins_max
          when 0
            ins[n, m][:point] = NONE
            $logger.debug "%-10s %5s" % ["I-#{n}-#{m}", "NONE"]
          when ins_mat
            jmp = "M-#{n}-#{m-1}"
            ins[n, m][:jump] = jmp
            prev_gap_ins_ext = str_pss[str_i].gap_ins_ext
            $logger.debug "%-10s %5s %10s" % ["I-#{n}-#{m}", "JUMP", jmp]
          when ins_ins
            ins[n, m][:point] = UP
            $logger.debug "%-10s %5s %10s" % ["I-#{n}-#{m}", "POINT", "D-#{n}-#{m-1}"]
          end

          # keep the record of a matrix and its indexes having maximum value so far
          max = [mat_max, del_max, ins_max].max
          if max >= max_s
            max_i = n
            max_j = m
            max_s = max
            max_m = case max
                    when mat_max then :mat
                    when del_max then :del
                    when ins_max then :ins
                    end
          end
        end
      end

      ProfileProfileLocalAlignmentAffineGap.new(@structural_profile, @sequence_profile,
                                                mat, del, ins, max_i, max_j, max_m)
    end

    def global_alignment_with_linear_gap_penalty(gap_del=100, gap_ins=100)
      str_pss   = @structural_profile.positions
      seq_pss   = @sequence_profile.positions
      str_plen  = @structural_profile.length
      seq_plen  = @sequence_profile.length
      stmat     = NArray.object(str_plen+1, seq_plen+1)

      # initialize score and trace matrix
      (0..str_plen).each do |str_i|
        (0..seq_plen).each do |seq_i|
          stmat[str_i, seq_i] = { :score => nil, :point => nil }
          if    (str_i == 0 && seq_i == 0)
            stmat[str_i, seq_i] = { :score => 0, :point => NONE }
          elsif (str_i >  0 && seq_i == 0)
            stmat[str_i, seq_i] = { :score => -gap_del * str_i, :point => LEFT }
          elsif (str_i == 0 && seq_i >  0)
            stmat[str_i, seq_i] = { :score => -gap_ins * seq_i, :point => UP }
          end
        end
      end

      (1..str_plen).each do |n|
        (1..seq_plen).each do |m|
          # calculate profile-profile match score
          mat_score = 0
          AMINO_ACIDS.split('').each do |aa|
            #w = seq_pss[seq_i].probability_of(aa)
            w = seq_pss[m-1].relative_probability_of(aa)
            s = str_pss[n-1].mat_score(aa == 'C' ? 'U' : aa)
            mat_score += (w * s)
          end
          mat_score = mat_score.round

          mat = stmat[n-1,  m-1][:score] + mat_score
          del = stmat[n-1,    m][:score] - gap_del
          ins = stmat[n,    m-1][:score] - gap_ins
          max = [mat, del, ins].max

          stmat[n, m][:score] = max
          stmat[n, m][:point] = case max
                                when mat then DIAG
                                when del then LEFT
                                when ins then UP
                                end
        end
      end

      ProfileProfileGlobalAlignmentLinearGap.new(@structural_profile, @sequence_profile,
                                                 stmat)
    end

    def global_alignment_with_affine_gap_penalty
      str_pss   = @structural_profile.positions
      seq_pss   = @sequence_profile.positions
      str_plen  = @structural_profile.length
      seq_plen  = @sequence_profile.length

      # create score, deletion, and insertion matrices
      mat = NArray.object(str_plen+1, seq_plen+1)
      del = NArray.object(str_plen+1, seq_plen+1)
      ins = NArray.object(str_plen+1, seq_plen+1)

      # 1. initialize score matrix and fill in the first row and column
      prev_gap_ins_ext = nil

      # 1.1 Range.each version
      (0..str_plen).each do |str_i|
        (0..seq_plen).each do |seq_i|
          mat[str_i, seq_i] = { :score => nil, :point => nil, :jump => nil }
          if    (str_i == 0 && seq_i == 0)
            mat[str_i, seq_i][:score] = 0
            mat[str_i, seq_i][:point] = NONE
          elsif (str_i == 1 && seq_i == 0)
            mat[str_i, seq_i][:score] = -str_pss[str_i-1].gap_del_open
            mat[str_i, seq_i][:point] = LEFT
          elsif (str_i >  1 && seq_i == 0)
            mat[str_i, seq_i][:score] = mat[str_i-1, seq_i][:score] - str_pss[str_i-1].gap_del_ext
            mat[str_i, seq_i][:point] = LEFT
          elsif (str_i == 0 && seq_i == 1)
            mat[str_i, seq_i][:score] = -str_pss[str_i].gap_ins_open
            mat[str_i, seq_i][:point] = UP
            prev_gap_ins_ext    = str_pss[str_i].gap_ins_ext
          elsif (str_i == 0 && seq_i >  1)
            mat[str_i, seq_i][:score] = mat[str_i, seq_i-1][:score] - prev_gap_ins_ext
            mat[str_i, seq_i][:point] = UP
          end
        end
      end

      # 2. initialize deletion and insertion matrices
      infinity = 1/0.0

      # 2.1 Range.each version
      (0..str_plen).each do |str_i|
        (0..seq_plen).each do |seq_i|
          del[str_i, seq_i] = { :score => (str_i == 0 || seq_i == 0 ? -infinity : nil), :point => nil, :jump => nil }
          ins[str_i, seq_i] = { :score => (str_i == 0 || seq_i == 0 ? -infinity : nil), :point => nil, :jump => nil }
        end
      end

      # 3. fill in match, deletion, and insertion matrices
      prev_gap_ins_ext = str_pss[0].gap_ins_ext

      # 3.1 Range.each version
      (1..str_plen).each do |n|
        (1..seq_plen).each do |m|
          str_i, seq_i = n-1, m-1

          # calculate profile-profile match score
          mat_score = 0
          AMINO_ACIDS.split('').each do |aa|
            #w = seq_pss[seq_i].probability_of(aa)
            w = seq_pss[seq_i].relative_probability_of(aa)
            s = str_pss[str_i].mat_score(aa == 'C' ? 'U' : aa)
            mat_score += (w * s)
          end
          mat_score = mat_score.round

          mat_mat = mat[n-1, m-1][:score] + mat_score
          mat_del = del[n-1, m-1][:score] + mat_score
          mat_ins = ins[n-1, m-1][:score] + mat_score
          mat_max = [mat_mat, mat_del, mat_ins].max
          mat[n, m][:score] = mat_max

          case mat_max
          when mat_mat
            mat[n, m][:point] = DIAG
            $logger.debug "%-10s %5s %10s" % ["M-#{n}-#{m}", "POINT", "M-#{n-1}-#{m-1}"]
          when mat_del
            jmp = "D-#{n-1}-#{m-1}"
            mat[n, m][:jump] = jmp
            $logger.debug "%-10s %5s %10s" % ["M-#{n}-#{m}", "JUMP", jmp]
          when mat_ins
            jmp = "I-#{n-1}-#{m-1}"
            mat[n, m][:jump] = jmp
            $logger.debug "%-10s %5s %10s" % ["M-#{n}-#{m}", "JUMP", jmp]
          end

          del_mat = mat[n-1, m][:score] - str_pss[str_i].gap_del_open
          del_del = del[n-1, m][:score] - str_pss[str_i].gap_del_ext
          del_max = [del_mat, del_del].max
          del[n, m][:score] = del_max

          case del_max
          when del_mat
            jmp = "M-#{n-1}-#{m}"
            del[n, m][:jump] = jmp
            $logger.debug "%-10s %5s %10s" % ["D-#{n}-#{m}", "JUMP", jmp]
          when del_del
            del[n, m][:point] = LEFT
            $logger.debug "%-10s %5s %10s" % ["D-#{n}-#{m}", "POINT", "D-#{n-1}-#{m}"]
          end

          ins_mat = mat[n, m-1][:score] - str_pss[str_i].gap_ins_open
          ins_ins = ins[n, m-1][:score] - prev_gap_ins_ext
          ins_max = [ins_mat, ins_ins].max
          ins[n, m][:score] = ins_max

          case ins_max
          when ins_mat
            jmp = "M-#{n}-#{m-1}"
            ins[n, m][:jump] = jmp
            prev_gap_ins_ext = str_pss[str_i].gap_ins_ext
            $logger.debug "%-10s %5s %10s" % ["I-#{n}-#{m}", "JUMP", jmp]
          when ins_ins
            ins[n, m][:point] = UP
            $logger.debug "%-10s %5s %10s" % ["I-#{n}-#{m}", "POINT", "D-#{n}-#{m-1}"]
          end
        end
      end

      ProfileProfileGlobalAlignmentAffineGap.new(@structural_profile, @sequence_profile,
                                                 mat, del, ins)
    end
  end
end
