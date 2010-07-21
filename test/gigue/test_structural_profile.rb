require File.join(File.dirname(__FILE__), '..', 'test_helper.rb')

class TestStructuralProfile < Test::Unit::TestCase

  include Gigue

  def setup
    @joy  = File.join(File.dirname(__FILE__), '..', 'd.240.1.1.tem')
    @esst = File.join(File.dirname(__FILE__), '..', 'ulla-mlr-logo-toccata-maskA.mat')
    @prof = StructuralProfile.new(@joy, @esst)
  end

  def test_length
    assert_equal(157, @prof.length)
  end

  def test_blosum_weighting
    i = 0
    @prof.weights.each do |k, v|
      assert_equal(0.2, v)
      assert_equal(@prof.sequences[i].code, k)
      i += 1
    end
  end

end
