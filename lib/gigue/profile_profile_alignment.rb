module Gigue
  class ProfileProfileAlignment

    attr_reader :structural_profile, :sequence_profile, :raw_score,
                :aligned_structural_profile_positions,
                :aligned_sequence_profile_positions

    def initialize(str_prf, seq_prf)
      @structural_profile = str_prf
      @sequence_profile   = seq_prf
    end

    def to_flatfile(options={})
      opts = {
        :os     => STDOUT,
        :type   => :pir,
        :width  => 70
      }.merge!(options)

      opts[:os] = File.open(opts[:os], 'w') if opts[:os].is_a? String

      # print aligned structural profile sequences
      @structural_profile.sequences.each_with_index do |pseq, psi|
        opts[:os].puts ">#{pseq.code}"
        opts[:os].puts "structure" if opts[:type] == :pir
        opts[:os].puts @aligned_structural_profile_positions.map_with_index { |p, pi|
          ((p == '-' ? '-' : (p.probe[psi] == 'J' ? 'C' : p.probe[psi])) +
           (pi > 0 && (pi+1) % opts[:width] == 0 ? "\n" : ''))
        }.join('')
      end

      # print aligned sequence profile sequences
      @sequence_profile.msa.sequences.each_with_index do |pseq, psi|
        opts[:os].puts ">#{pseq.code}"
        opts[:os].puts "sequence" if opts[:type] == :pir
        opts[:os].puts @aligned_sequence_profile_positions.map_with_index { |p, pi|
          ((p == '-' ? '-' : (p.probe[psi] == 'J' ? 'C' : p.probe[psi])) +
           (pi > 0 && (pi+1) % opts[:width] == 0 ? "\n" : ''))
        }.join('')
      end

      opts[:os].close if [File, String].include? opts[:os].class
    end
  end

  class ProfileProfileAlignmentLinearGap < ProfileProfileAlignment

    attr_reader :stmatrix

    def initialize(str_prf, seq_prf, stmatrix, i = nil, j = nil)
      super(str_prf, seq_prf)
      @stmatrix = stmatrix
      @aligned_structural_profile_positions, @aligned_sequence_profile_positions = traceback(i, j)
    end

    def raw_score
      @stmatrix[-1,-1][:score]
    end

    def traceback(i, j)
      str_pss     = @structural_profile.positions
      seq_pss     = @sequence_profile.positions
      str_cnt     = @structural_profile.sequences.size
      seq_cnt     = @sequence_profile.msa.sequences.size
      ali_str_pss = []
      ali_seq_pss = []

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
          ali_str_pss << str_pss[i-1]
          ali_seq_pss << seq_pss[j-1]
          i -= 1
          j -= 1
          $logger.debug "%-15s : %-4s : %s" % ["MAT[#{i}, #{j}]", "DIAG", "#{ali_str_pss[-1].probe} <=> #{ali_seq_pss[-1].probe}"]
        when LEFT
          ali_str_pss << str_pss[i-1]
          ali_seq_pss << SequenceProfilePosition.new('-'*seq_cnt)
          i -= 1
          $logger.debug "%-15s : %-4s : %s" % ["MAT[#{i}, #{j}]", "INS", "#{ali_str_pss[-1].probe} <=> #{ali_seq_pss[-1].probe}"]
        when UP
          ali_str_pss << StructuralProfilePosition.new('-'*str_cnt)
          ali_seq_pss << seq_pss[j-1]
          j -= 1
          $logger.debug "%-15s : %-4s : %s" % ["MAT[#{i}, #{j}]", "DEL", "#{ali_str_pss[-1].probe} <=> #{ali_seq_pss[-1].probe}"]
        else
          $logger.error "Something wrong with pointing stage at i: #{i}, j: #{j}"
          exit 1
        end
      end
      [ali_str_pss.reverse!, ali_seq_pss.reverse!]
    end
  end

  class ProfileProfileAlignmentAffineGap < ProfileProfileAlignment

    attr_reader :match_stmatrix, :insertion_stmatrix, :deletion_stmatrix

    def initialize(str_prf, seq_prf,
                   match_stmatrix, deletion_stmatrix, insertion_stmatrix,
                   i = nil, j = nil, max_mat = nil)
      super(str_prf, seq_prf)
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
      @aligned_structural_profile_positions, @aligned_sequence_profile_positions = traceback(i, j, max_mat)
    end

    def traceback(i, j, max_mat)
      str_pss     = @structural_profile.positions
      seq_pss     = @sequence_profile.positions
      str_cnt     = @structural_profile.sequences.size
      seq_cnt     = @sequence_profile.msa.sequences.size
      ali_str_pss = []
      ali_seq_pss = []
      pre_mat     = nil
      cur_mat     = nil

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
            ali_str_pss << str_pss[i-1]
            ali_seq_pss << seq_pss[j-1]
          when :del
            ali_str_pss << str_pss[i-1]
            ali_seq_pss << SequenceProfilePosition.new('-'*seq_cnt)
          when :ins
            ali_str_pss << StructuralProfilePosition.new('-'*str_cnt)
            ali_seq_pss << seq_pss[j-1]
          end

          jm, mi, jj  = cur_mat[i, j][:jump].split('-')
          i, j        = Integer(mi), Integer(jj)
          cur_mat     = case jm
                        when 'M'
                          $logger.debug "%-15s : %-4s : %s" % ["MAT[#{i}, #{j}]", "JUMP", "#{ali_str_pss[-1].probe} <=> #{ali_seq_pss[-1].probe}"]
                          pre_mat = :mat
                          @match_stmatrix
                        when 'D'
                          $logger.debug "%-15s : %-4s : %s" % ["DEL[#{i}, #{j}]", "JUMP", "#{ali_str_pss[-1].probe} <=> #{ali_seq_pss[-1].probe}"]
                          pre_mat = :del
                          @deletion_stmatrix
                        when 'I'
                          $logger.debug "%-15s : %-4s : %s" % ["INS[#{i}, #{j}]", "JUMP", "#{ali_str_pss[-1].probe} <=> #{ali_seq_pss[-1].probe}"]
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
          ali_str_pss << str_pss[i-1]
          ali_seq_pss << seq_pss[j-1]
          i -= 1
          j -= 1
          $logger.debug "%-15s : %-4s : %s" % ["MAT[#{i}, #{j}]", "DIAG", "#{ali_str_pss[-1].probe} <=> #{ali_seq_pss[-1].probe}"]
        when LEFT
          ali_str_pss << str_pss[i-1]
          ali_seq_pss << SequenceProfilePosition.new('-'*seq_cnt)
          i -= 1
          $logger.debug "%-15s : %-4s : %s" % ["DEL[#{i}, #{j}]", "LEFT", "#{ali_str_pss[-1].probe} <=> #{ali_seq_pss[-1].probe}"]
        when UP
          ali_str_pss << StructuralProfilePosition.new('-'*str_cnt)
          ali_seq_pss << seq_pss[j-1]
          j -= 1
          $logger.debug "%-15s : %-4s : %s" % ["INS[#{i}, #{j}]", "UP", "#{ali_str_pss[-1].probe} <=> #{ali_seq_pss[-1].probe}"]
        else
          $logger.error "Something wrong with pointing stage at i: #{i}, j: #{j}"
          exit 1
        end
      end
      [ali_str_pss.reverse!, ali_seq_pss.reverse!]
    end

  end
end
