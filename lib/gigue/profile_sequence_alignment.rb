module Gigue
  class ProfileSequenceAlignment

    attr_reader :profile, :sequence, :raw_score,
                :aligned_profile_positions,
                :aligned_amino_acids

    def initialize(profile, sequence)
      @profile  = profile
      @sequence = sequence
    end

    def to_flatfile(options = {})
      opts = {
        :os     => STDOUT,
        :type   => :pir,
        :width  => 70
      }.merge!(options)

      opts[:os] = File.open(opts[:os], 'w') if opts[:os].is_a? String

      # print aligned profile entries
      entries = @profile.entry_names
      entries.each_with_index do |entry, ei|
        opts[:os].puts ">#{entry}"
        opts[:os].puts "structure" if opts[:type] == :pir
        opts[:os].puts @aligned_profile_positions.map_with_index { |p, pi|
          ((p == '-' ? '-' : (p.probe[ei] == 'J' ? 'C' : p.probe[ei])) +
           (pi > 0 && (pi+1) % opts[:width] == 0 ? "\n" : ''))
        }.join('')
      end

      # print aligned query sequence
      opts[:os].puts ">#{sequence.code}"
      opts[:os].puts "sequence" if opts[:type] == :pir
      opts[:os].puts @aligned_amino_acids.map_with_index { |a, ai|
        a + (ai > 0 && (ai+1) % opts[:width] == 0 ? "\n" : '')
      }.join('')

      opts[:os].close if [File, String].include? opts[:os].class
    end

  end

  class ProfileSequenceAlignmentLinearGap < ProfileSequenceAlignment

    attr_reader :stmatrix

    def initialize(profile, sequence, stmatrix, i = nil, j = nil)
      super(profile, sequence)
      @stmatrix = stmatrix
      @aligned_profile_positions, @aligned_amino_acids = traceback(i, j)
    end

    def raw_score
      @stmatrix[-1,-1][:score]
    end

    def traceback(i, j)
      pss     = @profile.positions
      aas     = @sequence.amino_acids
      seq_cnt = @profile.seq_cnt
      ali_prf = []
      ali_seq = []

      if i.nil? || j.nil?
        i, j = @stmatrix.shape
        i -= 1
        j -= 1
        $logger.debug "No indexes provided, so use bottom right corner indexes (#{i}, #{j}) instead."
      end

      loop do
        case @stmatrix[i, j][:point]
        when NONE
          break
        when DIAG
          ali_prf << pss[i-1]
          ali_seq << aas[j-1]
          i -= 1
          j -= 1
          $logger.debug "%-15s : %-4s : %s" % ["MAT[#{i}, #{j}]", "DIAG", "#{ali_prf[-1] == '-' ? '-' * seq_cnt : ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
        when LEFT
          ali_prf << pss[i-1]
          ali_seq << '-'
          i -= 1
          $logger.debug "%-15s : %-4s : %s" % ["MAT[#{i}, #{j}]", "INS", "#{ali_prf[-1] == '-' ? '-' * seq_cnt : ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
        when UP
          ali_prf << '-'
          ali_seq << aas[j-1]
          j -= 1
          $logger.debug "%-15s : %-4s : %s" % ["MAT[#{i}, #{j}]", "DEL", "#{ali_prf[-1] == '-' ? '-' * seq_cnt : ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
        else
          $logger.error "Something wrong with pointing stage at i: #{i}, j: #{j}"
          exit 1
        end
      end
      [ali_prf.reverse!, ali_seq.reverse!]
    end

  end

  class ProfileSequenceAlignmentAffineGap < ProfileSequenceAlignment

    attr_reader :match_stmatrix, :insertion_stmatrix, :deletion_stmatrix

    def initialize(profile, sequence,
                   match_stmatrix, deletion_stmatrix, insertion_stmatrix,
                   i = nil, j = nil, max_mat = nil)
      super(profile, sequence)
      @match_stmatrix     = match_stmatrix
      @deletion_stmatrix  = deletion_stmatrix
      @insertion_stmatrix = insertion_stmatrix
      @raw_score          = if (!i.nil? && !j.nil? && !max_mat.nil?)
                              case max_mat
                              when :mat then @match_stmatrix[-1,-1][:score]
                              when :del then @deletion_stmatrix[-1,-1][:score]
                              when :ins then @insertion_stmatrix[-1,-1][:score]
                              end
                            else
                              [ @match_stmatrix[-1,-1][:score],
                                @deletion_stmatrix[-1,-1][:score],
                                @insertion_stmatrix[-1,-1][:score]
                              ].max
                            end
      @aligned_profile_positions, @aligned_amino_acids = traceback(i, j, max_mat)
    end

    def traceback(i, j, max_mat)
      pss     = @profile.positions
      aas     = @sequence.amino_acids
      seq_cnt = @profile.seq_cnt
      ali_prf = []
      ali_seq = []
      pre_mat = nil
      cur_mat = nil

      if i.nil? || j.nil?
        i, j = @match_stmatrix.shape
        i -= 1
        j -= 1
        $logger.debug "No indexes provided, so use bottom right corner indexes (#{i}, #{j}) instead."
      end

      cur_mat = if max_mat.nil?
                  case @raw_score
                  when @match_stmatrix[-1, -1][:score]
                    $logger.debug "START TRACEBACK FROM MAT MATRIX"
                    pre_mat = :mat
                    @match_stmatrix
                  when @deletion_stmatrix[-1, -1][:score]
                    $logger.debug "START TRACEBACK FROM DEL MATRIX"
                    pre_mat = :del
                    @deletion_stmatrix
                  when @insertion_stmatrix[-1, -1][:score]
                    $logger.debug "START TRACEBACK FROM INS MATRIX"
                    pre_mat = :ins
                    @insertion_stmatrix
                  else
                    $logger.error "Raw score doesn't match with any matrix"
                    exit 1
                  end
                else
                  case max_mat
                  when :mat
                    $logger.debug "START TRACEBACK FROM MAT MATRIX"
                    pre_mat = :mat
                    @match_stmatrix
                  when :del
                    $logger.debug "START TRACEBACK FROM DEL MATRIX"
                    pre_mat = :del
                    @deletion_stmatrix
                  when :ins
                    $logger.debug "START TRACEBACK FROM INS MATRIX"
                    pre_mat = :ins
                    @insertion_stmatrix
                  else
                    $logger.error "Maximum matrix tag doesn't match with any matrix"
                    exit 1
                  end
                end

      loop do
        until cur_mat[i, j][:jump].nil?
          case pre_mat
          when :mat
            ali_prf << pss[i-1]
            ali_seq << aas[j-1]
          when :del
            ali_prf << pss[i-1]
            ali_seq << '-'
          when :ins
            ali_prf << '-'
            ali_seq << aas[j-1]
          end

          jm, mi, jj  = cur_mat[i, j][:jump].split('-')
          i, j      = Integer(mi), Integer(jj)
          cur_mat = case jm
                    when 'M'
                      $logger.debug "%-15s : %-4s : %s" % ["MAT[#{i}, #{j}]", "JUMP", "#{ali_prf[-1] == '-' ? '-' * seq_cnt : ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
                      pre_mat = :mat
                      @match_stmatrix
                    when 'D'
                      $logger.debug "%-15s : %-4s : %s" % ["DEL[#{i}, #{j}]", "JUMP", "#{ali_prf[-1] == '-' ? '-' * seq_cnt : ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
                      pre_mat = :del
                      @deletion_stmatrix
                    when 'I'
                      $logger.debug "%-15s : %-4s : %s" % ["INS[#{i}, #{j}]", "JUMP", "#{ali_prf[-1] == '-' ? '-' * seq_cnt : ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
                      pre_mat = :ins
                      @insertion_stmatrix
                    else
                      $logger.warn "Something wrong in jumping step"
                      exit 1
                    end
        end

        point = cur_mat[i, j][:point]
        case point
        when NONE
          $logger.debug "FINISH TRACEBACK"
          break
        when DIAG
          ali_prf << pss[i-1]
          ali_seq << aas[j-1]
          i -= 1
          j -= 1
          $logger.debug "%-15s : %-4s : %s" % ["MAT[#{i}, #{j}]", "DIAG", "#{ali_prf[-1] == '-' ? '-' * seq_cnt : ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
        when LEFT
          ali_prf << pss[i-1]
          ali_seq << '-'
          i -= 1
          $logger.debug "%-15s : %-4s : %s" % ["DEL[#{i}, #{j}]", "LEFT", "#{ali_prf[-1] == '-' ? '-' * seq_cnt : ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
        when UP
          ali_prf << '-'
          ali_seq << aas[j-1]
          j -= 1
          $logger.debug "%-15s : %-4s : %s" % ["INS[#{i}, #{j}]", "UP", "#{ali_prf[-1] == '-' ? '-' * seq_cnt : ali_prf[-1].probe} <=> #{ali_seq[-1]}"]
        else
          $logger.error "Something wrong with pointing stage at i: #{i}, j: #{j}"
          exit 1
        end
      end
      [ali_prf.reverse!, ali_seq.reverse!]
    end

  end
end
