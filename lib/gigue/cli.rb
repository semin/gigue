module Gigue
  class CLI

    def self.execute(stdout, arguments=[])
      # default options
      $options = {
        :algorithm  => :auto,
        :weighting  => :va,
        :process    => 1,
        :purge      => 0.5
      }

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
          opts.on('-t', '--template FILENAME', String, 'set JOY template file') { |o|
            $options[:joytem] = o
          }
          opts.on('-e', '--essts FILENAME', String, 'set environment-specific substition tables file') { |o|
            $options[:essts] = o
          }
          opts.on('-w', '--weighting INTEGER', Integer,
                  'set weighting scheme',
                  '0    EqualWeighting  -- weighting each sequence equally',
                  '1    BlosumWeighting -- weighting scheme based on single linkage clustering',
                  '2    VAWeighting     -- Vingron and Argos weighting (default)') { |o|
            $options[:weighting] = case o
                                   when 0 then :equal
                                   when 1 then :blosum
                                   when 2 then :va
                                   else
                                     STDERR.puts "error: [--verbose|-v] #{o} is not supported (use 'gigue build -h' for help)"
                                     exit
                                   end

          }
          opts.on('-o', '--output FILENAME', String, 'set an output profile file (default: STDOUT)') { |o|
            $options[:output] = o
          }
          opts.on('-v', '--verbose INTEGER', Integer,
                  'show detailed console output',
                  '0    Error level (default)',
                  '1    Warning level',
                  '2    Information level',
                  '3    Debugging level') { |o|
            $options[:verbose] = case o
                                 when 0 then Logger::ERROR
                                 when 1 then Logger::WARN
                                 when 2 then Logger::INFO
                                 when 3 then Logger::DEBUG
                                 else
                                   STDERR.puts "error: [--verbose|-v] #{o} is not supported (use 'gigue build -h' for help)"
                                   exit
                                 end
            $logger.level = $options[:verbose]
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
          opts.on('-p', '--profile FILENAME', String, 'set target profile') { |o|
            $options[:profile] = o
          }
          opts.on('-s', '--sequence FILENAME', String, 'set query sequence(s)') { |o|
            $options[:sequence] = o
          }
          opts.on('-o', '--output FILENAME', String, 'set output file name (default: STDOUT)') { |o|
            $options[:output] = o
          }
          opts.on('-t', '--toprank INTEGER', Integer, 'output scoring information about top N HITs') { |o|
            $options[:toprank] = o
          }
          opts.on('-z', '--zcutoff FLOAT', Float, 'output scoring information about HITs with Z-scores > cutoff') { |o|
            $options[:zcutoff] = o
          }
          opts.on('-c', '--process INTEGER', Integer, 'set number of processes to use (default: 1)') { |o|
            $options[:process] = o
          }
          #opts.on('-r', '--purge FLOAT', Float, 'purge gap rich columns in a GIGUE profile (default: 0.5)') {
            #$options[:purge] = o
          #}
          opts.on('-a', '--algorithm INTEGER', Integer,
                  'set alignment algorithm',
                  '0      Auto (Global/Glocal) (default)',
                  '1      Global',
                  '2      Glocal',
                  '3      Local') { |o|
            $options[:algorithm] = case o
                                   when 0 then :auto
                                   when 1 then :global
                                   when 2 then :glocal
                                   when 3 then :local
                                   else
                                   STDERR.puts "error: [--algorithm|-a] #{o} is not supported (use 'gigue search -h' for help)"
                                   exit
                                 end
          }
          opts.on('-v', '--verbose INTEGER', Integer,
                  'show detailed console output',
                  '0      Error level (default)',
                  '1      Warning level',
                  '2      Information level',
                  '3      Debugging level') { |o|
            $options[:verbose] = case o
                                 when 0 then Logger::ERROR
                                 when 1 then Logger::WARN
                                 when 2 then Logger::INFO
                                 when 3 then Logger::DEBUG
                                 else
                                   STDERR.puts "error: [--verbose|-v] #{o} is not supported (use 'gigue search -h' for help)"
                                   exit
                                 end
            $logger.level = $options[:verbose]
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
          opts.on('-p', '--profile FILENAME', String, 'set target profile') { |o|
            $options[:profile] = o
          }
          opts.on('-s', '--sequence FILENAME', String, 'set query sequence(s)') { |o|
            $options[:sequence] = o
          }
          opts.on('-a', '--alignment FILENAME', String, 'set query alignment') { |o|
            $options[:alignment] = o
          }
          opts.on('-o', '--output FILENAME', String, 'set output file name (default: STDOUT)') { |o|
            $options[:output] = o
          }
          opts.on('-a', '--algorithm INTEGER', Integer,
                  'set alignment algorithm',
                  '0      Auto (Global/Glocal) (default)',
                  '1      Global',
                  '2      Glocal',
                  '3      Local') { |o|
            $options[:algorithm] = case o
                                   when 0 then :auto
                                   when 1 then :global
                                   when 2 then :glocal
                                   when 3 then :local
                                   else
                                   STDERR.puts "error: [--algorithm|-a] #{o} is not supported (use 'gigue search -h' for help)"
                                   exit
                                 end
          }
          opts.on('-v', '--verbose INTEGER', Integer,
                  'show detailed console output',
                  '0      Error level (default)',
                  '1      Warning level',
                  '2      Information level',
                  '3      Debugging level') { |o|
            $options[:verbose] = case o
                                 when 0 then Logger::ERROR
                                 when 1 then Logger::WARN
                                 when 2 then Logger::INFO
                                 when 3 then Logger::DEBUG
                                 else
                                   STDERR.puts "error: [--verbose|-v] #{o} is not supported (use 'gigue build -h' for help)"
                                   exit
                                 end
            $logger.level = $options[:verbose]
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
        STDERR.puts "error: #{e.message} (use -h for help)"
        exit
      end

      if subcommand.nil?
        STDERR.puts globalopts.banner
        exit
      end

      if (subcommand != 'build') && (subcommand != 'search') && (subcommand != 'align')
        STDERR.puts "error: first argument must be either 'build', 'search', or 'align' (use -h for help)"
        exit
      end

      begin
        subopts[subcommand].order!
      rescue => e
        STDERR.puts "error: #{e.message} (use 'gigue #{subcommand} -h' for help)"
        exit
      end

      # Perform required task
      case subcommand
      when 'build'
        if $options[:joytem].nil? && $options[:essts].nil?
          STDERR.puts subopts['build']
          exit
        end

        if $options[:joytem].nil?
          STDERR.puts "error: JOY template file must be provided (use 'gigue build -h' for help)"
          exit
        end

        if $options[:essts].nil?
          STDERR.puts "error: ESST file must be provided (use 'gigue build -h' for help)"
          exit
        end

        build_profile($options[:joytem], $options[:essts], $options[:weighting], $options[:output])
      when 'search'
        if $options[:profile].nil? && $options[:sequence].nil?
          STDERR.puts subopts['search']
          exit
        end

        if $options[:profile].nil?
          STDERR.puts "error: GIGUE profile must be provided (use 'gigue search -h' for help)"
          exit
        end

        if $options[:sequence].nil?
          STDERR.puts "error: sequence file must be provided (use 'gigue search -h' for help)"
          exit
        end

        search_profile($options[:profile], $options[:sequence], $options[:output])
      when 'align'
        if $options[:profile].nil? && $options[:sequence].nil?
          STDERR.puts subopts['align']
          exit
        end

        if $options[:profile].nil?
          STDERR.puts "error: GIGUE profile must be provided (use 'gigue align -h' for help)"
          exit
        end

        if $options[:sequence].nil?
          STDERR.puts "error: sequence file must be provided (use 'gigue align -h' for help)"
          exit
        end

        align_profile($options[:profile], $options[:sequence], $options[:output])
      end
    end

    def self.build_profile(tem, essts, wgt, out)
      $logger.debug "Building a GIGUE profile from #{tem} using #{essts}"

      unless File.exist? tem
        STDERR.puts "error: cannot find JOY template file, #{tem}"
        exit
      end

      unless File.exists? essts
        STDERR.puts "error: cannot find ESST file, #{essts}"
        exit
      end

      prf = StructuralProfile.create_from_joy_tem_and_essts(tem, essts, :weighting => wgt);
      out = out.nil? ? STDOUT : out
      prf.to_gig(out)
    end

    def self.search_profile(prf, db, out)
      $logger.debug "Searching #{prf} against #{db}"

      unless File.exist? prf
        STDERR.puts "error: cannot find GIGUE profile, #{prf}"
        exit
      end

      unless File.exist? db
        STDERR.puts "error: cannot find sequence file, #{db}"
        exit
      end

      stem  = File.basename(prf)
      prf   = FugueProfile.new(prf)
      ff    = Bio::FlatFile.auto(db)

      hits = Parallel.map(ff.entries, :in_processes => $options[:process]) do |ent|
        seq   = Sequence.new(ent.aaseq, ent.entry_id, ent.definition)
        psa   = ProfileSequenceAligner.new(prf, seq)
        ali   = psa.local_alignment_affine_gap
        cols  = [
          stem, prf.length, seq.code, seq.length,
          ali.raw_score, ali.reverse_score, ali.calculate_z_score(75)
        ]
      end

      shits = hits.sort_by { |h| h[-1] }.reverse
      out   = out.nil? ? STDOUT : File.open(out, 'w')
      hdr   = "#%-15s %10s  %-15s %10s %10s %10s %10s" %
                %w[Prof ProfLen Seq SeqLen RawScore RevScore Z-score]

      out.puts hdr

      shits.each do |h|
        hit_fmt = " %-15s %10d  %-15s %10d %10d %10d %10.3f"
        out.puts hit_fmt % h
      end

      out.close if out.is_a? File
    end

    def self.align_profile(prf, db, out)
      $logger.debug "Aligning #{prf} against #{db}"

      unless File.exist? prf
        STDERR.puts "error: cannot find GIGUE profile, #{prf}"
        exit
      end

      unless File.exist? db
        STDERR.puts "error: cannot find sequence file, #{db}"
        exit
      end

      stem  = File.basename(prf)
      prf   = FugueProfile.new(prf)
      ent   = Bio::FlatFile.auto(db).next_entry
      seq   = Sequence.new(ent.aaseq.to_s, ent.entry_id, ent.definition)
      psa   = ProfileSequenceAligner.new(prf, seq)
      ali   = psa.global_alignment_affine_gap
      out   = out.nil? ? STDOUT : File.open(out, 'w')

      ali.to_flatfile(:out => out)
    end

  end
end
