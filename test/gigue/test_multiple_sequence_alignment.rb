require File.join(File.dirname(__FILE__), '..', 'test_helper.rb')

class TestMultipleSequenceAlignment < Test::Unit::TestCase

  include Gigue

  def setup
    @file = File.join(File.dirname(__FILE__), '..', 'test.blast.out6')
    @msa  = MultipleSequenceAlignment.create_from_psiblast_output_style6(@file)
  end

  def test_instance_type
    assert_instance_of(MultipleSequenceAlignment, @msa)
  end

  def test_first_sequence_code
    assert_equal('test|1-101', @msa.sequences[0].code)
  end

  def test_first_sequence_data
    assert_equal('VRKSIGRIVTMKRNSR---NLEEIKPYLFRAIEESYY--K--L--D--K--R-------I--P--K--A--IHVVAVTEDLDI-VSRGRT----F-P-HG-I-S-K---E-T-A-Y-S-E-SVKLLQKI-L--E------EDE-RKIRRIGVRFSKFI', @msa.sequences[0].data)
  end

  def test_last_sequence_code
    assert_equal('UniRef90_Q74H50|249-347', @msa.sequences[-1].code)
  end

  def test_last_sequence_data
    assert_equal('--KSVGHSMTLDRDLT---ARRDILKYLLQLSEMVGR--R--A--RRYG--V-------A--G--K--T--VHLTIRYADFTT-VGKQQT----R-N-QA-TNS-T---E-E-I-Y-A-E-AVKILDTF-----------ELL-QPVRLLGVRITNL-', @msa.sequences[-1].data)
  end

  def test_to_sequence_profile
    sqp = @msa.to_sequence_profile
    assert_instance_of(SequenceProfile, sqp)
  end
end
