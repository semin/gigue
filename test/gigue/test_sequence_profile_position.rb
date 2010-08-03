require File.join(File.dirname(__FILE__), '..', 'test_helper.rb')

class TestSequenceProfilePosition < Test::Unit::TestCase

  include Gigue

  def setup
    @spp = SequenceProfilePosition.new('--MMIIV-MVV-')
  end

  def test_frequency_of
    assert_equal(3, @spp.frequency_of('M'))
    assert_equal(2, @spp.frequency_of('I'))
    assert_equal(3, @spp.frequency_of('V'))
    assert_equal(4, @spp.frequency_of('-'))
  end

  def test_probability_of
    size = @spp.probe.size.to_f
    assert_equal(3/size, @spp.probability_of('M'))
    assert_equal(2/size, @spp.probability_of('I'))
    assert_equal(3/size, @spp.probability_of('V'))
    assert_equal(4/size, @spp.probability_of('-'))
  end

end
