require 'optparse'

module Gigue
  class CLI

    def self.execute(stdout, arguments=[])
      $options    = {}
      globalopts  = OptionParser.new do |opts|
        opts.banner = <<-BANNER
GIGUE: A Sequence-Structure Homology Detection Method Using Environment-Specific Substitution Tables

Usage: #{File.basename($0)} [--version] [--help|-h] COMMAND [ARGS]

GIGUE commands are:
    build     build profile(s) from JOY template(s) using ESSTs
    search    search query sequence(s) or alignment against target profile(s)
    align     align query sequence(s) or alignment with target profile(s)

See 'gigue COMMAND [--help|-h]' for more information on a specific command.
        BANNER
        opts.on('-h', '--help', 'show this help message') {
          stdout.puts opts.banner
          exit
        }
        opts.on('--version', 'display program version information') {
          stdout.puts "GIGUE #{VERSION}"
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
          opts.on('-l', '--list FILENAME', String, 'set list file for multiple JOY templates') { |o|
            $options[:joytemlist] = o
          }
          opts.on('-s', '--essts FILENAME', String, 'set environment-specific substition tables file') { |o|
            $options[:essts] = o
          }
          opts.on('-m', '--master FLOAT', Float, 'set weight for a master sequence') { |o|
            $options[:master] = o
          }
          opts.on('-p', '--prob', 'output probabilities instead of log-odd scores') { |o|
            $options[:prob] = true
          }
          opts.on('-w', '--weight', Integer,
                  'set weighting scheme',
                  '0      EqualWeight     -- weight each sequence equally',
                  '1      BlosumWeight    -- weighting scheme based on single linkage clustering (default)',
                  '2      VAWeight        -- Vingron and Argos weight') { |o|
            $options[:weight] = o || 1
          }
          opts.on('-o', '--output FILENAME', String, 'set an output profile file') { |o|
            $options[:output] = o
          }
          opts.on('-h', '--help', 'show this help message') {
            stdout.puts opts
            exit
          }
          opts.on('-v', '--verbose INTEGER', Integer, 'show detailed console output', ' ') { |o|
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
        end,

        'search' => OptionParser.new do |opts|
          opts.banner = <<-BANNER
GIGUE: A Sequence-Structure Homology Detection Method Using Environment-Specific Substitution Tables

Usage: #{File.basename($0)} search [options]

Options:
          BANNER
          opts.on('-s', '--sequence FILENAME', String, 'load query sequence(s)') { |o|
            $options[:sequence] = o
          }
          opts.on('-a', '--alignment FILENAME', String, 'load query alignment') { |o|
            $options[:alignment] = o
          }
          opts.on('-p', '--profile FILENAME', String, 'load target profile') { |o|
            $options[:profile] = o
          }
          opts.on('-l', '--list FILENAME', String, 'load list of target profile(s)') { |o|
            $options[:list] = o
          }
          opts.on('-n', '--noz', 'do NOT calculate Z-score') {
            $options[:noz] = true
          }
          opts.on('-t', '--toprank INTEGER', Integer, 'output scoring information about top N HITs') { |o|
            $options[:toprank] = o
          }
          opts.on('-z', '--zcutoff FLOAT', Float, 'output scoring information about top N HITs') { |o|
            $options[:zcutoff] = o
          }
          opts.on('-v', '--verbose INTEGER', Integer, 'show detailed console output', ' ') { |o|
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
        end,

        'align' => OptionParser.new do |opts|
          opts.banner = <<-BANNER
GIGUE: A Sequence-Structure Homology Detection Method Using Environment-Specific Substitution Tables

Usage: #{File.basename($0)} align [options]

Options:
          BANNER
          opts.on('-s', '--sequence FILENAME', String, 'load query sequence(s)') { |o|
            $options[:sequence] = o
          }
          opts.on('-a', '--alignment FILENAME', String, 'load query alignment') { |o|
            $options[:alignment] = o
          }
          opts.on('-p', '--profile FILENAME', String, 'load target profile') { |o|
            $options[:profile] = o
          }
          opts.on('-l', '--list FILENAME', String, 'load list of target profile(s)') { |o|
            $options[:list] = o
          }
          opts.on('-v', '--verbose INTEGER', Integer, 'show detailed console output', ' ') { |o|
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
        end
      }

      begin
        arguments   = globalopts.order!
        subcommand  = arguments.shift
      rescue => e
        STDERR.puts "error: #{e.message} (use -h for help)"
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
        if $options[:joytem].nil? && $options[:joytemlist].nil?
          STDERR.puts "error: JOY template(s) must be provided to generate profile(s) (use 'gigue build -h' for help)"
          exit
        end
        if !$options[:joytem].nil? && !$options[:joytemlist].nil?
          STDERR.puts "error: both options [--template|-t] and [--list|-l] provided. (use 'gigue build -h' for help)"
          exit
        end
        if $options[:essts].nil?
          STDERR.puts "error: ESSTs must be provided to generate profile(s) (use 'gigue build -h' for help)"
          exit
        end
        build_profile
      when 'search'
        if $options[:sequence].nil? && $options[:alignment].nil?
          STDERR.puts "error: query seqeunce or alignment must be provided to search againt profile(s) (use 'gigue search -h' for help)"
          exit
        end
        if !$options[:sequence].nil? && !$options[:alignment].nil?
          STDERR.puts "error: both options [--sequence|-s] and [--alignment|-a] provided. (use 'gigue search -h' for help)"
          exit
        end
        if $options[:profile].nil?
          STDERR.puts "error: target profile(s) must be provided to search againt (use 'gigue search -h' for help)"
          exit
        end
        search_profile
      end
    end

    def self.build_profile
      $logger.debug 'building profile(s)...'
    end

    def self.search_profile
      $logger.debug 'searching profile(s)...'
    end

  end
end
