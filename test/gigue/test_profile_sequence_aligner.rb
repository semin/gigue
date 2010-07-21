require File.join(File.dirname(__FILE__), '..', 'test_helper.rb')

class TestProfileSequenceAligner < Test::Unit::TestCase
  include Gigue

  def setup
    $logger.level = Logger::WARN
    file          = File.join(File.dirname(__FILE__), '..', 'd.240.1.1.fug')
    @profile      = FugueProfile.new(file)
    @sequence     = Sequence.new('test', 'VRKSIGRIVTMKRNSRNLEEIKPYLFRAIEESYYKLDKRIPKAIHVVAVTEDLDIVSRGRTFPHGISKETAYSESVKLLQKILEEDERKIRRIGVRFSKFI')
  end

  def test_global_alignment_with_affine_gap_penalty
    psa = ProfileSequenceAligner.new(@profile, @sequence)
    gal = psa.global_alignment_with_affine_gap_penalty
    assert_instance_of(ProfileSequenceAlignmentAffineGap, gal)
  end

  def test_global_alignment_with_linear_gap_penalty
    psa = ProfileSequenceAligner.new(@profile, @sequence)
    gal = psa.global_alignment_with_linear_gap_penalty
    assert_instance_of(ProfileSequenceAlignmentLinearGap, gal)
  end

  def test_local_alignment_with_linear_gap_penalty
    psa = ProfileSequenceAligner.new(@profile, @sequence)
    gal = psa.local_alignment_with_linear_gap_penalty
    assert_instance_of(ProfileSequenceAlignmentLinearGap, gal)
  end

end
