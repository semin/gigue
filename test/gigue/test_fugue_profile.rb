require File.join(File.dirname(__FILE__), '..', 'test_helper.rb')

class TestFugueProfile < Test::Unit::TestCase

  include Gigue

  def setup
    file      = File.join(File.dirname(__FILE__), '..', 'd.240.1.1.fug')
    @profile  = FugueProfile.new(file)
  end

  def test_command
    assert_equal("melody -t d.240.1.1.tem -c ../classdef_canonical_env5.dat -s ../toccata/MaskA.dat", @profile.command)
  end

  def test_length
    assert_equal(157, @profile.length)
  end

  def test_number_of_sequences
    assert_equal(5, @profile.sequences.size)
  end

  def test_weighting
    assert_equal(1, @profile.weighting)
  end

  def test_entry_names
    assert_equal('106363',  @profile.sequences[0].code)
    assert_equal('90374',   @profile.sequences[1].code)
    assert_equal('99680',   @profile.sequences[2].code)
    assert_equal('90378',   @profile.sequences[3].code)
    assert_equal('106701',  @profile.sequences[4].code)
  end

  def test_entry_weights
    assert_equal(0.2,  @profile.weights[@profile.sequences[0].code])
    assert_equal(0.2,  @profile.weights[@profile.sequences[1].code])
    assert_equal(0.2,  @profile.weights[@profile.sequences[2].code])
    assert_equal(0.2,  @profile.weights[@profile.sequences[3].code])
    assert_equal(0.2,  @profile.weights[@profile.sequences[4].code])
  end

  def test_multiple_factor
    assert_equal(10.0, @profile.multiple_factor)
  end

  def test_format
    assert_equal('0-FUGUE', @profile.format)
  end

  def test_row_symbols
    assert_equal('ACDEFGHIKLMNPQRSTVWYJU'.split(''), @profile.row_symbols)
  end

  def test_col_symbols
    assert_equal('ACDEFGHIKLMNPQRSTVWYJ'.split(''), @profile.col_symbols)
  end

  def test_env_symbols
    assert_equal('HEPCAaSsOoNn'.split(''), @profile.env_symbols)
  end

  def test_gap_open_ins_term
    assert_equal(100, @profile.gap_open_ins_term)
  end

  def test_gap_open_del_term
    assert_equal(100, @profile.gap_open_del_term)
  end

  def test_gap_ext_ins_term
    assert_equal(100, @profile.gap_ext_ins_term)
  end

  def test_gap_ext_del_term
    assert_equal(100, @profile.gap_ext_del_term)
  end

  def test_aa_colnames
    assert_equal('ACDEFGHIKLMNPQRSTVWYJU'.split(''), @profile.aa_colnames)
  end

  def test_gap_colnames
    assert_equal('InsO InsE DelO DelE COIL HNcp HCcp HInt SNcp SCcp SInt NRes  Ooi  Acc'.split, @profile.gap_colnames)
  end

  def test_env_colnames
    assert_equal('H   E   P   C   A   a   S   s   O   o   N   n'.split, @profile.env_colnames)
  end

  def test_positions
    ps = @profile.positions
    assert_equal(157,     ps.size)
    assert_equal('-P---', ps[0].probe)
    assert_equal(-20,     ps[0].mat_score('A'))
    assert_equal(-50,     ps[0].mat_score('U'))
    assert_equal(100,     ps[0].gap_score('InsO'))
    assert_equal(11,      ps[0].gap_score('Acc'))
    assert_equal(0,       ps[0].env_score('H'))
    assert_equal(100,     ps[0].env_score('n'))
    assert_equal('--L--', ps[-1].probe)
    assert_equal(-20,     ps[-1].mat_score('A'))
    assert_equal(-20,     ps[-1].mat_score('U'))
    assert_equal(100,     ps[-1].gap_score('InsO'))
    assert_equal(11,      ps[-1].gap_score('Acc'))
    assert_equal(0,       ps[-1].env_score('H'))
    assert_equal(100,     ps[-1].env_score('n'))
  end
end
