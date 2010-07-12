module Gigue
  class ProfileSequenceAlignment

    attr_reader :profile, :sequence, :score, :point

    def initialize(prf, seq, score, point)
      @profile  = prf
      @sequence = seq
      @score    = score
      @point    = point
    end

    def to_fasta(wrap = 70)
      pss     = @profile.positions
      aas     = @sequence.split('')
      ali_prf = []
      ali_seq = []
      i, j    = @point.shape

      loop do
        p = @point[i-1, j-1]
        break if p == NONE
        s = @score[i-1, j-1]
        if (p == DIAG)
          ali_prf << pss[i-2]
          ali_seq << aas[j-2]
          i += -1
          j += -1
        elsif (p == LEFT)
          ali_prf << pss[i-2]
          ali_seq << '-'
          i += -1
        elsif (p == UP)
          ali_prf << '-'
          ali_seq << aas[j-2]
          j += -1
        else
          warn "Something wrong!"
          exit 1
        end
      end

      ali_prf.reverse!
      ali_seq.reverse!

      # print aligned profile entries
      entries = @profile.entry_names
      entries.each_with_index do |entry, ei|
        puts ">#{entry}"
        puts ali_prf.map_with_index { |p, pi|
          ((p == '-' ? '-' : p.probe[ei]) +
           (pi > 0 && (pi+1) % wrap == 0 ? "\n" : ''))
        }.join('')
      end

      # print aligned query sequence
      puts ">test" # #{seq.name} should come
      puts ali_seq.map_with_index { |a, ai|
        a + (ai > 0 && (ai+1) % wrap == 0 ? "\n" : '')
      }.join('')
    end

  end
end
