$:.unshift(File.dirname(__FILE__)) unless $:.include?(File.dirname(__FILE__)) || $:.include?(File.expand_path(File.dirname(__FILE__)))

require 'bio'
require 'logger'
require 'narray'
require 'inline'
require 'ostruct'
require 'pathname'
require 'stringio'
require 'optparse'
require 'parallel'

require_relative 'array_extensions'
require_relative 'gigue/esst'
require_relative 'gigue/essts'
require_relative 'gigue/joy_tem'
require_relative 'gigue/sequence'
require_relative 'gigue/substitution_table'
require_relative 'gigue/sequence_profile'
require_relative 'gigue/sequence_profile_position'
require_relative 'gigue/fugue_profile'
require_relative 'gigue/structural_profile'
require_relative 'gigue/structural_profile_position'
require_relative 'gigue/sequence_sequence_aligner'
require_relative 'gigue/sequence_sequence_alignment'
require_relative 'gigue/profile_sequence_aligner'
require_relative 'gigue/profile_sequence_alignment'
require_relative 'gigue/profile_profile_aligner'
require_relative 'gigue/profile_profile_alignment'
require_relative 'gigue/multiple_sequence_alignment'
require_relative 'gigue/multiple_sequence_alignment_column'
require_relative 'gigue/environment_type_class_definition'

module Gigue
  # logger setting
  $logger = Logger.new(STDERR)
  $logger.level = Logger::ERROR

  # global constants
  VERSION = '0.0.1'
  NONE, UP, LEFT, DIAG = 0, 1, 2, 3
end
