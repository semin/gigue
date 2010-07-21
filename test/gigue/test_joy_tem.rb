require File.join(File.dirname(__FILE__), '..', 'test_helper.rb')

class TestJoyTem < Test::Unit::TestCase

  include Gigue

  def setup
    @file   = File.join(File.dirname(__FILE__), '..', 'd.240.1.1.tem')
    @joytem = JoyTem.new(@file)
  end

  def test_file
    assert_equal(@file, @joytem.file)
  end

  def test_first_entry
    code  = '106363'
    desc  = 'sequence'
    data  = '-Q----SFSEED-SFK-KC-SSEV--EAKNKIEELLASLLNRV-CQDG----R-K--PHTVRLIIRR--YSSEKHYGRESRQ--C--PIPSHVIQKLGTGNY--DVMTPM-VDILMKLF-RNMVNV-KMPF--HLTLLSVCFCNLK-----------'
    assert(@joytem.entries.has_key?(code))
    assert(@joytem.entries[code].has_key?(desc))
    assert_equal(data, @joytem.entries[code][desc])
  end

  def test_last_entry
    code  = '106701'
    desc  = 'Ooi number'
    data  = '----012343433-2-2-233--3322443453344333-2--221----1333----33445553322------22233333--2---2221----10----4--5434543454343-2--21111-24--34345555432-------------'
    assert(@joytem.entries.has_key?(code))
    assert(@joytem.entries[code].has_key?(desc))
    assert_equal(data, @joytem.entries[code][desc])
  end

end
