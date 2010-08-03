require File.join(File.dirname(__FILE__), '..', 'test_helper.rb')

class TestSequenceProfile < Test::Unit::TestCase

  include Gigue

  def setup
    @blt = File.join(File.dirname(__FILE__), '..', 'test.blast.out6')
    @msa = MultipleSequenceAlignment.create_from_psiblast_output_style6(@blt)
    @sqp = @msa.to_sequence_profile
  end

  def test_instance_type
    assert_instance_of(SequenceProfile, @sqp)
  end

  def test_length
    assert_equal(@sqp.msa.length, @sqp.length)
    assert_equal(@sqp.msa.length, @sqp.positions.length)
  end

  def test_first_position_frequency
    # VI--II---IL---II-I-I------II-V---V--I-III----IV---III--V-VIV----I-------I---VI----------I-----I--VII-----V-V------I---V---I----------I---I--------V-I------------II--V---------I----------I----V-I------V--------I---------I---I---VI-I--II-II---V-I-V-----
    assert_equal(19,  @sqp.positions[0].frequency_of('V'))
    assert_equal(46,  @sqp.positions[0].frequency_of('I'))
    assert_equal(1,   @sqp.positions[0].frequency_of('L'))
    assert_equal(185, @sqp.positions[0].frequency_of('-'))
  end

  def test_first_position_probability
    seqcnt = @sqp.msa.sequences.size.to_f
    assert_equal(19/seqcnt,   @sqp.positions[0].probability_of('V'))
    assert_equal(46/seqcnt,   @sqp.positions[0].probability_of('I'))
    assert_equal(1/seqcnt,    @sqp.positions[0].probability_of('L'))
    assert_equal(185/seqcnt,  @sqp.positions[0].probability_of('-'))
  end

  def test_last_position_frequency
    # I-V-V-V-----V-I------VV-------L------V---V-----------V--L----V---I-----V---------------------V--------V-V-------VV------V--IVIV----V--V-------V-----I-------VVV-I-------VM--V-------V------V--V----V-----VI---V-I-VV-----V--------V---------------V--------
    assert_equal(9,   @sqp.positions[-1].frequency_of('I'))
    assert_equal(38,  @sqp.positions[-1].frequency_of('V'))
    assert_equal(2,   @sqp.positions[-1].frequency_of('L'))
    assert_equal(1,   @sqp.positions[-1].frequency_of('M'))
    assert_equal(201, @sqp.positions[-1].frequency_of('-'))
  end

  def test_last_position_probability
    seqcnt = @sqp.msa.sequences.size.to_f
    assert_equal(9/seqcnt,    @sqp.positions[-1].probability_of('I'))
    assert_equal(38/seqcnt,   @sqp.positions[-1].probability_of('V'))
    assert_equal(2/seqcnt,    @sqp.positions[-1].probability_of('L'))
    assert_equal(1/seqcnt,    @sqp.positions[-1].probability_of('M'))
    assert_equal(201/seqcnt,  @sqp.positions[-1].probability_of('-'))
  end

end
