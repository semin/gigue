require File.join(File.dirname(__FILE__), '..', 'spec_helper.rb')
include Gigue

describe ProfileSequenceAligner do

  before(:all) do
    $logger.level = Logger::WARN

    file  = File.join(File.dirname(__FILE__), '..', 'd.240.1.1.fug')
    @prf  = FugueProfile.new(file)
    @seq  = Sequence.new('VRKSIGRIVTMKRNSRNLEEIKPYLFRAIEESYYKLDKRIPKAIHVVAVTEDLDIVSRGRTFPHGISKETAYSESVKLLQKILEEDERKIRRIGVRFSKFI', 'test')
    @psa  = ProfileSequenceAligner.new(@prf, @seq)
  end

  it "has profile" do
    @psa.profile.should == @prf
  end

  it "has sequence" do
    @psa.sequence.should == @seq
  end

  it "#global_alignment_with_linear_gap_penalty returns an instance of ProfileSequenceAlignmentLinearGap class" do
    ali = @psa.global_alignment_with_linear_gap_penalty
    ali.should be_an_instance_of(ProfileSequenceAlignmentLinearGap)
  end

  it "#global_alignment_with_affine_gap_penalty returns an instance of ProfileSequenceAlignmentAffineGap class" do
    ali = @psa.global_alignment_with_affine_gap_penalty
    ali.should be_an_instance_of(ProfileSequenceAlignmentAffineGap)
  end

  it "#local_alignment_with_linear_gap_penalty returns an instance of ProfileSequenceAlignmentLinearGap class" do
    ali = @psa.local_alignment_with_linear_gap_penalty
    ali.should be_an_instance_of(ProfileSequenceAlignmentLinearGap)
  end

  it "#local_alignment_with_affine_gap_penalty returns an instance of ProfileSequenceAlignmentAffineGap class" do
    ali = @psa.local_alignment_with_affine_gap_penalty
    ali.should be_an_instance_of(ProfileSequenceAlignmentAffineGap)
  end

end
