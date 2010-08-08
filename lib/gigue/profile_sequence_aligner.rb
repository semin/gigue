module Gigue
  class ProfileSequenceAligner

    attr_reader :profile, :sequence

    def initialize(prf, seq)
      @profile  = prf
      @sequence = seq
    end

    def local_alignment_with_linear_gap_penalty(gap_del = 100, gap_ins = 100)
      pss   = @profile.positions
      aas   = @sequence.amino_acids
      plen  = @profile.length
      slen  = @sequence.length
      max_s = 0
      max_i = nil
      max_j = nil
      stmat = NArray.object(plen+1, slen+1)

      # initialize score and trace matrix
      (0..plen).each do |n|
        (0..slen).each do |m|
          stmat[n, m] = { :score => 0, :point => nil }
          if    (n == 0 && m == 0)
            stmat[n, m][:point] = NONE
          elsif (n >  0 && m == 0)
            stmat[n, m][:point] = LEFT
          elsif (n == 0 && m >  0)
            stmat[n, m][:point] = UP
          end
        end
      end

      # fill in score and trace matrix
      (1..plen).each do |n|
        (1..slen).each do |m|
          mat = stmat[n-1,  m-1][:score] + pss[n-1].mat_score(aas[m-1])
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
      ProfileSequenceAlignmentLinearGap.new(@profile, @sequence, stmat, max_i, max_j)
    end

    def local_alignment_with_affine_gap_penalty
      pss   = @profile.positions
      aas   = @sequence.amino_acids
      plen  = @profile.length
      slen  = @sequence.length
      max_s = 0
      max_i = nil
      max_j = nil
      max_m = nil

      # create score, deletion, and insertion matrices
      mat = NArray.object(plen+1, slen+1)
      del = NArray.object(plen+1, slen+1)
      ins = NArray.object(plen+1, slen+1)

      # initialize score matrix and fill in the first row and column
      prev_gap_ins_ext = nil
      (0..plen).each do |pi|
        (0..slen).each do |si|
          mat[pi, si] = { :score => 0, :point => nil, :jump => nil }
          if    (pi == 0 && si == 0)
            mat[pi, si][:point] = NONE
          elsif (pi == 1 && si == 0)
            mat[pi, si][:score] = -pss[pi-1].gap_del_open
            mat[pi, si][:point] = LEFT
          elsif (pi >  1 && si == 0)
            mat[pi, si][:score] = mat[pi-1, si][:score] - pss[pi-1].gap_del_ext
            mat[pi, si][:point] = LEFT
          elsif (pi == 0 && si == 1)
            mat[pi, si][:score] = -pss[pi].gap_ins_open
            mat[pi, si][:point] = UP
            prev_gap_ins_ext    = pss[pi].gap_ins_ext
          elsif (pi == 0 && si >  1)
            mat[pi, si][:score] = mat[pi, si-1][:score] - prev_gap_ins_ext
            mat[pi, si][:point] = UP
          end
        end
      end

      # initialize deletion and insertion matrices
      infinity = 1/0.0
      (0..plen).each do |pi|
        (0..slen).each do |si|
          del[pi, si] = { :score => (pi == 0 || si == 0 ? -infinity : 0), :point => nil, :jump => nil }
          ins[pi, si] = { :score => (pi == 0 || si == 0 ? -infinity : 0), :point => nil, :jump => nil }
        end
      end

      # fill in match, deletion, and insertion matrices
      prev_gap_ins_ext = pss[0].gap_ins_ext
      (1..plen).each do |n|
        (1..slen).each do |m|
          pi, ai  = n-1, m-1

          # update match matrix
          mat_mat = mat[n-1, m-1][:score] + pss[pi].mat_score(aas[ai])
          mat_del = del[n-1, m-1][:score] + pss[pi].mat_score(aas[ai])
          mat_ins = ins[n-1, m-1][:score] + pss[pi].mat_score(aas[ai])
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
          del_mat = mat[n-1, m][:score] - pss[pi].gap_del_open
          del_del = del[n-1, m][:score] - pss[pi].gap_del_ext
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
          ins_mat = mat[n, m-1][:score] - pss[pi].gap_ins_open
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
            prev_gap_ins_ext = pss[pi].gap_ins_ext
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

      ProfileSequenceAlignmentAffineGap.new(@profile, @sequence,
                                            mat, del, ins, max_i, max_j, max_m)
    end

    def global_alignment_with_linear_gap_penalty(gap_del = 100, gap_ins = 100)
      pss   = @profile.positions
      aas   = @sequence.amino_acids
      plen  = @profile.length
      slen  = @sequence.length
      stmat = NArray.object(plen+1, slen+1) # matrix for score and traceback

      # initialize score and trace matrix
      (0..plen).each do |pi|
        (0..slen).each do |si|
          stmat[pi, si] = { :score => nil, :point => nil }
          if    (pi == 0 && si == 0)
            stmat[pi, si] = { :score => 0, :point => NONE }
          elsif (pi >  0 && si == 0)
            stmat[pi, si] = { :score => -gap_del * pi, :point => LEFT }
          elsif (pi == 0 && si >  0)
            stmat[pi, si] = { :score => -gap_ins * si, :point => UP }
          end
        end
      end

      (1..plen).each do |n|
        (1..slen).each do |m|
          mat = stmat[n-1,  m-1][:score] + pss[n-1].mat_score(aas[m-1])
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
      ProfileSequenceAlignmentLinearGap.new(@profile, @sequence, stmat)
    end

    def global_alignment_with_affine_gap_penalty
      pss   = @profile.positions
      aas   = @sequence.amino_acids
      plen  = @profile.length
      slen  = @sequence.length

      # create score, deletion, and insertion matrices
      mat = NArray.object(plen+1, slen+1)
      del = NArray.object(plen+1, slen+1)
      ins = NArray.object(plen+1, slen+1)

      # initialize score matrix and fill in the first row and column
      prev_gap_ins_ext = nil
      (0..plen).each do |pi|
        (0..slen).each do |si|
          mat[pi, si] = { :score => nil, :point => nil, :jump => nil }
          if    (pi == 0 && si == 0)
            mat[pi, si][:score] = 0
            mat[pi, si][:point] = NONE
          elsif (pi == 1 && si == 0)
            mat[pi, si][:score] = -pss[pi-1].gap_del_open
            mat[pi, si][:point] = LEFT
          elsif (pi >  1 && si == 0)
            mat[pi, si][:score] = mat[pi-1, si][:score] - pss[pi-1].gap_del_ext
            mat[pi, si][:point] = LEFT
          elsif (pi == 0 && si == 1)
            mat[pi, si][:score] = -pss[pi].gap_ins_open
            mat[pi, si][:point] = UP
            prev_gap_ins_ext    = pss[pi].gap_ins_ext
          elsif (pi == 0 && si >  1)
            mat[pi, si][:score] = mat[pi, si-1][:score] - prev_gap_ins_ext
            mat[pi, si][:point] = UP
          end
        end
      end

      # initialize deletion and insertion matrices
      infinity = 1/0.0
      (0..plen).each do |pi|
        (0..slen).each do |si|
          del[pi, si] = { :score => (pi == 0 || si == 0 ? -infinity : nil), :point => nil, :jump => nil }
          ins[pi, si] = { :score => (pi == 0 || si == 0 ? -infinity : nil), :point => nil, :jump => nil }
        end
      end

      # fill in match, deletion, and insertion matrices
      prev_gap_ins_ext = pss[0].gap_ins_ext
      (1..plen).each do |n|
        (1..slen).each do |m|
          pi, ai  = n-1, m-1

          mat_mat = mat[n-1, m-1][:score] + pss[pi].mat_score(aas[ai])
          mat_del = del[n-1, m-1][:score] + pss[pi].mat_score(aas[ai])
          mat_ins = ins[n-1, m-1][:score] + pss[pi].mat_score(aas[ai])
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

          del_mat = mat[n-1, m][:score] - pss[pi].gap_del_open
          del_del = del[n-1, m][:score] - pss[pi].gap_del_ext
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


          ins_mat = mat[n, m-1][:score] - pss[pi].gap_ins_open
          ins_ins = ins[n, m-1][:score] - prev_gap_ins_ext
          ins_max = [ins_mat, ins_ins].max
          ins[n, m][:score] = ins_max

          case ins_max
          when ins_mat
            jmp = "M-#{n}-#{m-1}"
            ins[n, m][:jump] = jmp
            prev_gap_ins_ext = pss[pi].gap_ins_ext
            $logger.debug "%-10s %5s %10s" % ["I-#{n}-#{m}", "JUMP", jmp]
          when ins_ins
            ins[n, m][:point] = UP
            $logger.debug "%-10s %5s %10s" % ["I-#{n}-#{m}", "POINT", "D-#{n}-#{m-1}"]
          end
        end
      end
      ProfileSequenceAlignmentAffineGap.new(@profile,
                                            @sequence,
                                            mat, del, ins)
    end

  end
end
