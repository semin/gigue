$:.unshift(File.dirname(__FILE__)) unless
  $:.include?(File.dirname(__FILE__)) || $:.include?(File.expand_path(File.dirname(__FILE__)))

require 'narray'
require File.join(File.dirname(__FILE__), 'enumerable_extensions')
require File.join(File.dirname(__FILE__), 'gigue', 'profile')
require File.join(File.dirname(__FILE__), 'gigue', 'profile_position')
require File.join(File.dirname(__FILE__), 'gigue', 'profile_sequence_aligner')
require File.join(File.dirname(__FILE__), 'gigue', 'profile_sequence_alignment')

module Gigue
  VERSION = '0.0.1'
  NONE, UP, LEFT, DIAG = 0, 1, 2, 3
end
