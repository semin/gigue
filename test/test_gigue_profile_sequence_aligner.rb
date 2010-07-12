require File.join(File.dirname(__FILE__), "test_helper.rb")
require 'gigue/profile_sequence_aligner'

class TestProfileSequenceAligner < Test::Unit::TestCase
  include Gigue

  def setup
    file      = File.join(File.dirname(__FILE__), "d.240.1.1.fug")
    @profile  = FugueProfile.new(file)
    @seq      = 'VRKSIGRIVTMKRNSRNLEEIKPYLFRAIEESYYKLDKRIPKAIHVVAVTEDLDIVSRGRTFPHGISKETAYSESVKLLQKILEEDERKIRRIGVRFSKFI'
  end

  def test_global_align
    ps_aligner = ProfileSequenceAligner.new(@profile, @seq)
    ps_alignment = ps_aligner.global_alignment
  end
end
