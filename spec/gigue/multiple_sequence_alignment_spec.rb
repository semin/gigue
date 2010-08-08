require File.join(File.dirname(__FILE__), '..', 'spec_helper.rb')
include Gigue

describe MultipleSequenceAlignment do

  before(:all) do
    @file = File.join(File.dirname(__FILE__), '..', 'test.blast.out6')
    @msa  = MultipleSequenceAlignment.create_from_psiblast_output_style6(@file)
  end

  it "is an instance of MultipleSequenceAlignment class" do
    @msa.should be_an_instance_of(MultipleSequenceAlignment)
  end

  it "has correct information for the first sequence" do
    @msa.sequences[0].code.should == 'test|1-101'
    @msa.sequences[0].data.should == 'VRKSIGRIVTMKRNSR---NLEEIKPYLFRAIEESYY--K--L--D--K--R-------I--P--K--A--IHVVAVTEDLDI-VSRGRT----F-P-HG-I-S-K---E-T-A-Y-S-E-SVKLLQKI-L--E------EDE-RKIRRIGVRFSKFI'
  end

  it "has correct information for the last sequence" do
    @msa.sequences[-1].code.should == 'UniRef90_Q74H50|249-347'
    @msa.sequences[-1].data.should == '--KSVGHSMTLDRDLT---ARRDILKYLLQLSEMVGR--R--A--RRYG--V-------A--G--K--T--VHLTIRYADFTT-VGKQQT----R-N-QA-TNS-T---E-E-I-Y-A-E-AVKILDTF-----------ELL-QPVRLLGVRITNL-'
  end

  it "#length returns a length of multiple sequence alignment" do
    @msa.length.should == '--KSVGHSMTLDRDLT---ARRDILKYLLQLSEMVGR--R--A--RRYG--V-------A--G--K--T--VHLTIRYADFTT-VGKQQT----R-N-QA-TNS-T---E-E-I-Y-A-E-AVKILDTF-----------ELL-QPVRLLGVRITNL-'.length
  end

  it "#columns returns a collection of MultipleSequenceAlignmentColumn" do
    @msa.columns[0].should be_an_instance_of(MultipleSequenceAlignmentColumn)
    @msa.columns[0].probe.should == 'VI--II---IL---II-I-I------II-V---V--I-III----IV---III--V-VIV----I-------I---VI----------I-----I--VII-----V-V------I---V---I----------I---I--------V-I------------II--V---------I----------I----V-I------V--------I---------I---I---VI-I--II-II---V-I-V-----'
  end

  it "#to_sequence_profile returns an instnace of SequenceProfile" do
    @msa.to_sequence_profile.should be_an_instance_of(SequenceProfile)
  end

end
