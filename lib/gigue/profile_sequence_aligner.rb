module Gigue
  class ProfileSequenceAligner

    def initialize(prf, seq, algo = :global)
      @profile  = prf
      @sequence = seq
    end

    def global_alignment
      pss     = @profile.positions
      aas     = @sequence.split('')
      prf_len = @profile.length
      seq_len = @sequence.length

      score = NMatrix.int(prf_len+1, seq_len+1).fill!(0)
      point = NMatrix.int(prf_len+1, seq_len+1).fill!(0)

      point[0, 0]           = NONE
      point[1..prf_len, 0]  = LEFT
      point[0, 1..seq_len]  = UP

      score[1, 0] = pss[0].gap_del_open
      (2..prf_len).each { |pi| score[pi, 0] = score[pi-1, 0] + pss[pi-1].gap_del_ext }
      score[0, 1] = pss[0].gap_ins_open
      (2..seq_len).each { |si| score[0, si] = score[0, si-1] + pss[0].gap_ins_ext }

      prev_stat = nil

      (0...aas.size).each do |ai|
        aa = aas[ai]
        (0...pss.size).each do |pi|
          ps = pss[pi]
          mat_score = score[pi,   ai] + ps.mat_score(aa)
          ins_score = score[pi+1, ai] + (prev_stat && prev_stat == :ins ? ps.gap_ins_ext : ps.gap_ins_open)
          del_score = score[pi, ai+1] + (prev_stat && prev_stat == :del ? ps.gap_del_ext : ps.gap_del_open)

          if (mat_score >= ins_score)
            if (mat_score >= del_score)
              score[pi+1, ai+1] = mat_score
              point[pi+1, ai+1] = DIAG
            else
              score[pi+1, ai+1] = del_score
              point[pi+1, ai+1] = LEFT
            end
          else
            if (ins_score > del_score)
              score[pi+1, ai+1] = ins_score
              point[pi+1, ai+1] = UP
            else
              score[pi+1, ai+1] = del_score
              point[pi+1, ai+1] = LEFT
            end
          end
        end
      end
      ProfileSequenceAlignment.new(@profile, @sequence, score, point)
    end

    def local_align(prf, seq)
    end

    def global_align(prf, seq)
    end
  end
end
