module Gigue
  class EnvironmentTypeClassDefinition

    attr_reader :file, :environments

    def initialize(file)
      @file         = file
      @environments = []

      IO.readlines(@file).each_with_index do |line, li|
        line.chomp!
        if    line =~ /^#/ || line =~ /^\s*$/
          next
        elsif (elems = line.split(';')).size == 5
          @environments << (env = OpenStruct.new)
          @environments[-1].name        = elems[0].strip
          @environments[-1].values      = elems[1].strip.split('')
          @environments[-1].labels      = elems[2].strip.split('')
          @environments[-1].constrained = elems[3].strip == 'T' ? true : false
          @environments[-1].silent      = elems[4].strip == 'T' ? true : false
        else
          $logger.error "Cannot parse line #{li+1}: #{line} in #{@file}"
          exit 1
        end
      end
    end

  end
end
