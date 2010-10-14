module Gigue
  class CLI

    def self.execute(stdout, arguments=[])
      # default options
      options = {
        :output     => STDOUT,
        :algorithm  => :auto,
        :weighting  => :va,
        :process    => 1,
        :purge      => 0.5,
        :verbose    => Logger::ERROR,
      }

      # option setting
      globalopts = OptionParser.new do |opts|
        opts.banner = <<-BANNER
GIGUE: A Sequence-Structure Homology Detection Method Using Environment-Specific Substitution Tables

Usage: #{File.basename($0)} [--version] [--help|-h] COMMAND [ARGS]

GIGUE commands are:
    build     build profile(s) from JOY template(s) using ESSTs
    search    search query sequence(s) or alignment against target profile(s)
    align     align query sequence(s) or alignment with target profile(s)

See 'gigue COMMAND [--help|-h]' for more information on a specific command.
        BANNER
        opts.on('--version', 'display program version information') {
          stdout.puts "GIGUE #{VERSION}"
          exit
        }
        opts.on('-h', '--help', 'show this help message') {
          stdout.puts opts.banner
          exit
        }
      end

      subopts = {
        'build' => OptionParser.new do |opts|
          opts.banner = <<-BANNER
GIGUE: A Sequence-Structure Homology Detection Method Using Environment-Specific Substitution Tables

Usage: #{File.basename($0)} build [options]

Options:
          BANNER
          opts.on('-t', '--template FILENAME',  String, 'set JOY template file') { |o| options[:joytem] = o }
          opts.on('-e', '--essts FILENAME',     String, 'set environment-specific substition tables file') { |o| options[:essts] = o }
          opts.on('-o', '--output FILENAME',    String, 'set an output profile file (default: STDOUT)') { |o| options[:output] = o }
          opts.on('-w', '--weighting INTEGER',  Integer,
                  'set weighting scheme',
                  '0    EqualWeighting  -- weighting each sequence equally',
                  '1    BlosumWeighting -- weighting scheme based on single linkage clustering',
                  '2    VAWeighting     -- Vingron and Argos weighting (default)') { |o|
            options[:weighting] = case o
                                  when 0 then :equal
                                  when 1 then :blosum
                                  when 2 then :va
                                  else
                                    $logger.error "[--verbose|-v] #{o} is not supported (use 'gigue build -h' for help)"
                                    exit
                                  end

          }
          opts.on('-v', '--verbose INTEGER', Integer,
                  'show detailed console output',
                  '0    Error level (default)',
                  '1    Warning level',
                  '2    Information level',
                  '3    Debugging level') { |o|
            options[:verbose] = case o
                                when 0 then Logger::ERROR
                                when 1 then Logger::WARN
                                when 2 then Logger::INFO
                                when 3 then Logger::DEBUG
                                else
                                  $logger.error "[--verbose|-v] #{o} is not supported (use 'gigue build -h' for help)"
                                  exit
                                end
            $logger.level = options[:verbose]
          }
          opts.on('-h', '--help', 'show this help message') {
            stdout.puts opts
            exit
          }
        end,

        'search' => OptionParser.new do |opts|
          opts.banner = <<-BANNER
GIGUE: A Sequence-Structure Homology Detection Method Using Environment-Specific Substitution Tables

Usage: #{File.basename($0)} search [options]

Options:
          BANNER
          opts.on('-p', '--profile FILENAME',   String, 'set target profile') { |o| options[:profile] = o }
          opts.on('-s', '--sequence FILENAME',  String, 'set query sequence(s)') { |o| options[:sequence] = o }
          opts.on('-o', '--output FILENAME',    String, 'set output file name (default: STDOUT)') { |o| options[:output] = o }
          opts.on('-t', '--toprank INTEGER',    Integer, 'output scoring information about top N HITs') { |o| options[:toprank] = o }
          opts.on('-z', '--zcutoff FLOAT',      Float, 'output scoring information about HITs with Z-scores > cutoff') { |o| options[:zcutoff] = o }
          opts.on('-c', '--process INTEGER',    Integer, 'set number of processes to use (default: 1)') { |o| options[:process] = o }
          opts.on('-a', '--algorithm INTEGER',  Integer,
                  'set alignment algorithm',
                  '0      Global/Glocal (default)',
                  '2      Local') { |o|
            options[:algorithm] = case o
                                  when 0 then :auto
                                  when 1 then :local
                                  else
                                    $logger.error "[--algorithm|-a] #{o} is not supported (use 'gigue search -h' for help)"
                                    exit
                                  end
          }
          opts.on('-v', '--verbose INTEGER', Integer,
                  'show detailed console output',
                  '0      Error level (default)',
                  '1      Warning level',
                  '2      Information level',
                  '3      Debugging level') { |o|
            options[:verbose] = case o
                                when 0 then Logger::ERROR
                                when 1 then Logger::WARN
                                when 2 then Logger::INFO
                                when 3 then Logger::DEBUG
                                else
                                  $logger.error "[--verbose|-v] #{o} is not supported (use 'gigue search -h' for help)"
                                  exit
                                end
            $logger.level = options[:verbose]
          }
          opts.on('-h', '--help', 'show this help message') {
            stdout.puts opts
            exit
          }
        end,

        'align' => OptionParser.new do |opts|
          opts.banner = <<-BANNER
GIGUE: A Sequence-Structure Homology Detection Method Using Environment-Specific Substitution Tables

Usage: #{File.basename($0)} align [options]

Options:
          BANNER
          opts.on('-p', '--profile FILENAME',   String, 'set target profile') { |o| options[:profile] = o }
          opts.on('-s', '--sequence FILENAME',  String, 'set query sequence(s)') { |o| options[:sequence] = o }
          opts.on('-a', '--alignment FILENAME', String, 'set query alignment') { |o| options[:alignment] = o }
          opts.on('-o', '--output FILENAME',    String, 'set output file name (default: STDOUT)') { |o| options[:output] = o }
          opts.on('-a', '--algorithm INTEGER',  Integer,
                  'set alignment algorithm',
                  '0      Global/Glocal (default)',
                  '2      Local') { |o|
            options[:algorithm] = case o
                                  when 0 then :auto
                                  when 1 then :local
                                  else
                                    $logger.error "[--algorithm|-a] #{o} is not supported (use 'gigue align -h' for help)"
                                    exit
                                  end
          }
          opts.on('-v', '--verbose INTEGER', Integer,
                  'show detailed console output',
                  '0      Error level (default)',
                  '1      Warning level',
                  '2      Information level',
                  '3      Debugging level') { |o|
            options[:verbose] = case o
                                when 0 then Logger::ERROR
                                when 1 then Logger::WARN
                                when 2 then Logger::INFO
                                when 3 then Logger::DEBUG
                                else
                                  $logger.error "[--verbose|-v] #{o} is not supported (use 'gigue build -h' for help)"
                                  exit
                                end
            $logger.level = options[:verbose]
          }
          opts.on('-h', '--help', 'show this help message') {
            stdout.puts opts
            exit
          }
        end
      }

      begin
        arguments   = globalopts.order!
        subcommand  = arguments.shift
      rescue => e
        $logger.error "#{e.message} (use -h for help)"
        exit
      end

      if subcommand.nil?
        stdout.puts globalopts.banner
        exit
      end

      if (subcommand != 'build') && (subcommand != 'search') && (subcommand != 'align')
        $logger.error "first argument must be either 'build', 'search', or 'align' (use -h for help)"
        exit
      end

      begin
        subopts[subcommand].order!
      rescue => e
        $logger.error "#{e.message} (use 'gigue #{subcommand} -h' for help)"
        exit
      end

      # Perform required task
      case subcommand
      when 'build'
        if options[:joytem].nil? && options[:essts].nil?
          stdout.puts subopts['build']
          exit
        end

        if options[:joytem].nil?
          $logger.error "JOY template file must be provided (use 'gigue build -h' for help)"
          exit
        end

        if options[:essts].nil?
          $logger.error "ESST file must be provided (use 'gigue build -h' for help)"
          exit
        end

        build_profile(options)
      when 'search'
        if options[:profile].nil? && options[:sequence].nil?
          stdout.puts subopts['search']
          exit
        end

        if options[:profile].nil?
          $logger.error "profile must be provided (use 'gigue search -h' for help)"
          exit
        end

        if options[:sequence].nil?
          $logger.error "sequence file must be provided (use 'gigue search -h' for help)"
          exit
        end

        search_profile(options)
      when 'align'
        if options[:profile].nil? && options[:sequence].nil?
          stdout.puts subopts['align']
          exit
        end

        if options[:profile].nil?
          $logger.error "profile must be provided (use 'gigue align -h' for help)"
          exit
        end

        if options[:sequence].nil?
          $logger.error "sequence file must be provided (use 'gigue align -h' for help)"
          exit
        end

        align_profile(options)
      end
    end

    def self.build_profile(opts)
      $logger.debug "Building a profile from #{opts[:joytem]} using #{opts[:essts]}"

      unless File.exist? opts[:joytem]
        $logger.error "cannot find JOY template file, #{opts[:joytem]}"
        exit
      end

      unless File.exists? opts[:essts]
        $logger.error "cannot find ESST file, #{opts[:essts]}"
        exit
      end

      prf = StructuralProfile.create_from_joy_tem_and_essts(opts[:joytem], opts[:essts], :weighting => opts[:weighting]);
      out = opts[:output].nil? ? STDOUT : opts[:output]
      prf.to_gig(out)
    end

    def self.search_profile(opts)
      $logger.debug "Searching #{opts[:profile]} against #{opts[:sequence]}"

      unless File.exist? opts[:profile]
        $logger.error "cannot find profile, #{opts[:profile]}"
        exit
      end

      unless File.exist? opts[:sequence]
        $logger.error "cannot find sequence file, #{opts[:sequence]}"
        exit
      end

      bsn   = File.basename(opts[:profile])
      prf   = FugueProfile.new(opts[:profile])
      ff    = Bio::FlatFile.auto(opts[:sequence])
      hits  = Parallel.map(ff.entries, :in_processes => opts[:process]) do |ent|
        seq = Sequence.new(ent.aaseq, ent.entry_id, ent.definition)
        psa = ProfileSequenceAligner.new(prf, seq)
        ali = if opts[:algorithm] == :local
                  psa.local_alignment_affine_gap
                else
                  psa.global_alignment_affine_gap
                end
        [ bsn, prf.length, seq.code, seq.length, ali.raw_score, ali.reverse_score, ali.calculate_z_score(75) ]
      end

      shits = hits.sort_by { |h| h[-1] }.reverse
      out   = out.nil? ? STDOUT : File.open(out, 'w')
      hdr   = "#%-15s %10s  %-15s %10s %10s %10s %10s" % %w[Prof ProfLen Seq SeqLen RawScore RevScore Z-score]

      out.puts hdr

      shits.each do |h|
        hit_fmt = " %-15s %10d  %-15s %10d %10d %10d %10.3f"
        out.puts hit_fmt % h
      end

      out.close if out.is_a? File
    end

    def self.align_profile(opts)
      $logger.debug "Aligning #{opts[:profile]} against #{opts[:sequence]}"

      unless File.exist? opts[:profile]
        $logger.error "cannot find profile, #{opts[:profile]}"
        exit
      end

      unless File.exist? opts[:sequence]
        $logger.error "cannot find sequence file, #{opts[:sequence]}"
        exit
      end

      bsn = File.basename(opts[:profile])
      prf = FugueProfile.new(opts[:profile])
      ent = Bio::FlatFile.auto(opts[:sequence]).next_entry
      seq = Sequence.new(ent.aaseq.to_s, ent.entry_id, ent.definition)
      psa = ProfileSequenceAligner.new(prf, seq)
      ali = if opts[:algorithm] == :local
              psa.local_alignment_affine_gap
            else
              psa.global_alignment_affine_gap
            end
      out = out.nil? ? STDOUT : File.open(out, 'w')
      ali.to_flatfile(:out => out)
    end

  end
end
