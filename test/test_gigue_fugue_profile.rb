require File.join(File.dirname(__FILE__), "test_helper.rb")
require 'gigue/fugue_profile'

class TestFugueProfile < Test::Unit::TestCase
  def setup
    file = File.join(File.dirname(__FILE__), "d.240.1.1.fug")
    @profile = FugueProfile.new(file)
  end

  def test_command
    assert_equal("melody -t d.240.1.1.tem -c ../classdef_canonical_env5.dat -s ../toccata/MaskA.dat", @profile.command)
  end

  def test_length
    assert_equal(157, @profile.length)
  end

  def test_seq_cnt
    assert_equal(5, @profile.seq_cnt)
  end

  def test_weighting
    assert_equal(1, @profile.weighting)
  end

  def test_entry_names
    assert_equal('106363',  @profile.entry_names[0])
    assert_equal('90374',   @profile.entry_names[1])
    assert_equal('99680',   @profile.entry_names[2])
    assert_equal('90378',   @profile.entry_names[3])
    assert_equal('106701',  @profile.entry_names[4])
  end

  def test_entry_weights
    assert_equal(0.2,  @profile.entry_weights[0])
    assert_equal(0.2,  @profile.entry_weights[1])
    assert_equal(0.2,  @profile.entry_weights[2])
    assert_equal(0.2,  @profile.entry_weights[3])
    assert_equal(0.2,  @profile.entry_weights[4])
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

  def test_probes
    assert_equal('-P---', @profile.probes[0])
    assert_equal('--L--', @profile.probes[-1])
    assert_equal(157, @profile.probes.size)
  end
end
