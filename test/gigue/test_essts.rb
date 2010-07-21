require File.join(File.dirname(__FILE__), '..', 'test_helper.rb')

class TestEssts < Test::Unit::TestCase

  include Gigue

  def setup
    file    = File.join(File.dirname(__FILE__), '..', 'subst-logo-smooth-toccata-maskA.mat')
    @essts  = Essts.new(file)
  end

  def test_type
    assert(:logo, @essts.type)
  end

end

