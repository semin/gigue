$:.unshift(File.dirname(__FILE__)) unless
  $:.include?(File.dirname(__FILE__)) || $:.include?(File.expand_path(File.dirname(__FILE__)))

require 'bio'
require 'narray'
require 'logger'
require 'inline'
require 'ostruct'
require 'pathname'
require 'stringio'

require File.join(File.dirname(__FILE__), 'array_extensions')
require File.join(File.dirname(__FILE__), 'enumerable_extensions')
require File.join(File.dirname(__FILE__), 'gigue', 'esst')
require File.join(File.dirname(__FILE__), 'gigue', 'essts')
require File.join(File.dirname(__FILE__), 'gigue', 'joy_tem')
require File.join(File.dirname(__FILE__), 'gigue', 'sequence')
require File.join(File.dirname(__FILE__), 'gigue', 'substitution_table')
require File.join(File.dirname(__FILE__), 'gigue', 'sequence_profile')
require File.join(File.dirname(__FILE__), 'gigue', 'sequence_profile_position')
require File.join(File.dirname(__FILE__), 'gigue', 'fugue_profile')
require File.join(File.dirname(__FILE__), 'gigue', 'structural_profile')
require File.join(File.dirname(__FILE__), 'gigue', 'structural_profile_position')
require File.join(File.dirname(__FILE__), 'gigue', 'sequence_sequence_aligner')
require File.join(File.dirname(__FILE__), 'gigue', 'sequence_sequence_alignment')
require File.join(File.dirname(__FILE__), 'gigue', 'profile_sequence_aligner')
require File.join(File.dirname(__FILE__), 'gigue', 'profile_sequence_alignment')
require File.join(File.dirname(__FILE__), 'gigue', 'profile_profile_aligner')
require File.join(File.dirname(__FILE__), 'gigue', 'profile_profile_alignment')
require File.join(File.dirname(__FILE__), 'gigue', 'multiple_sequence_alignment')
require File.join(File.dirname(__FILE__), 'gigue', 'multiple_sequence_alignment_column')
require File.join(File.dirname(__FILE__), 'gigue', 'environment_type_class_definition')

$logger = Logger.new(STDOUT)
$logger.level = Logger::DEBUG

module Gigue
  VERSION = '0.0.1'
  NONE, UP, LEFT, DIAG = 0, 1, 2, 3
  #AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWYJ'
end
