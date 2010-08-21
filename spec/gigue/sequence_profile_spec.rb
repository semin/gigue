require_relative '../spec_helper.rb'
include Gigue

describe SequenceProfile do

  before(:all) do
    @blt = File.join(File.dirname(__FILE__), '..', 'test.blast.out6')
    @msa = MultipleSequenceAlignment.create_from_psiblast_output_style6(@blt)
    @sqp = @msa.to_sequence_profile
  end

  it "is an instance of SequenceProfile" do
    @sqp.should be_an_instance_of(SequenceProfile)
  end

  it "#length returns the length of sequence profile" do
    @sqp.length.should == @msa.length
  end

  it "has correct frequency information for the first position" do
    # VI--II---IL---II-I-I------II-V---V--I-III----IV---III--V-VIV----I-------I---VI----------I-----I--VII-----V-V------I---V---I----------I---I--------V-I------------II--V---------I----------I----V-I------V--------I---------I---I---VI-I--II-II---V-I-V-----
    @sqp.positions[0].raw_frequency_of('V').should == 19
    @sqp.positions[0].raw_frequency_of('I').should == 46
    @sqp.positions[0].raw_frequency_of('L').should == 1
    @sqp.positions[0].raw_frequency_of('-').should == 185
  end

  it "has correct probability information for the first position" do
    seqcnt = @sqp.depth.to_f
    @sqp.positions[0].raw_probability_of('V').should == 19/seqcnt
    @sqp.positions[0].raw_probability_of('I').should == 46/seqcnt
    @sqp.positions[0].raw_probability_of('L').should == 1/seqcnt
    @sqp.positions[0].raw_probability_of('-').should == 185/seqcnt
  end

  it "has correct frequency information for the last position" do
    # I-V-V-V-----V-I------VV-------L------V---V-----------V--L----V---I-----V---------------------V--------V-V-------VV------V--IVIV----V--V-------V-----I-------VVV-I-------VM--V-------V------V--V----V-----VI---V-I-VV-----V--------V---------------V--------
    @sqp.positions[-1].raw_frequency_of('I').should == 9
    @sqp.positions[-1].raw_frequency_of('V').should == 38
    @sqp.positions[-1].raw_frequency_of('L').should == 2
    @sqp.positions[-1].raw_frequency_of('M').should == 1
    @sqp.positions[-1].raw_frequency_of('-').should == 201
  end

  it "has correct probability information for the last position" do
    seqcnt = @sqp.depth.to_f
    @sqp.positions[-1].raw_probability_of('I').should == 9/seqcnt
    @sqp.positions[-1].raw_probability_of('V').should == 38/seqcnt
    @sqp.positions[-1].raw_probability_of('L').should == 2/seqcnt
    @sqp.positions[-1].raw_probability_of('M').should == 1/seqcnt
    @sqp.positions[-1].raw_probability_of('-').should == 201/seqcnt
  end

end

