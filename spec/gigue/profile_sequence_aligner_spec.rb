require_relative '../spec_helper.rb'
include Gigue

describe ProfileSequenceAligner do

  before(:all) do
    $logger.level = Logger::WARN

    file  = File.join(File.dirname(__FILE__), '..', 'd.240.1.1.fug')
    @prf  = FugueProfile.new(file)
    @seq  = Sequence.new('VRKSIGRIVTMKRNSRNLEEIKPYLFRAIEESYYKLDKRIPKAIHVVAVTEDLDIVSRGRTFPHGISKETAYSESVKLLQKILEEDERKIRRIGVRFSKFI', 'test')
    @psa  = ProfileSequenceAligner.new(@prf, @seq)
  end

  it "has a profile" do
    @psa.structural_profile.should == @prf
  end

  it "has a sequence" do
    @psa.sequence.should == @seq
  end

  it "#local_alignment_linear_gap returns an instance of ProfileSequenceLocalAlignmentLinearGap class" do
    ali = @psa.local_alignment_linear_gap
    ali.should be_an_instance_of(ProfileSequenceLocalAlignmentLinearGap)
  end

  it "#local_alignment_linear_gap_rb returns an instance of ProfileSequenceLocalAlignmentLinearGap class" do
    ali = @psa.local_alignment_linear_gap_rb
    ali.should be_an_instance_of(ProfileSequenceLocalAlignmentLinearGap)
  end

  it "#local_alignment_linear_gap_cpp returns an instance of ProfileSequenceLocalAlignmentLinearGap class" do
    ali = @psa.local_alignment_linear_gap_cpp
    ali.should be_an_instance_of(ProfileSequenceLocalAlignmentLinearGap)
  end

  it "#local_alignment_affine_gap returns an instance of ProfileSequenceLocalAlignmentAffineGap class" do
    ali = @psa.local_alignment_affine_gap
    ali.should be_an_instance_of(ProfileSequenceLocalAlignmentAffineGap)
  end

  it "#local_alignment_affine_gap_rb returns an instance of ProfileSequenceLocalAlignmentAffineGap class" do
    ali = @psa.local_alignment_affine_gap_rb
    ali.should be_an_instance_of(ProfileSequenceLocalAlignmentAffineGap)
  end

  it "#local_alignment_affine_gap_cpp returns an instance of ProfileSequenceLocalAlignmentAffineGap class" do
    ali = @psa.local_alignment_affine_gap_cpp
    ali.should be_an_instance_of(ProfileSequenceLocalAlignmentAffineGap)
  end

  it "#global_alignment_linear_gap_cpp returns an instance of ProfileSequenceGlobalAlignmentLinearGap class" do
    ali = @psa.global_alignment_linear_gap
    ali.should be_an_instance_of(ProfileSequenceGlobalAlignmentLinearGap)
  end

  it "#global_alignment_linear_gap_rb returns an instance of ProfileSequenceGlobalAlignmentLinearGap class" do
    ali = @psa.global_alignment_linear_gap_rb
    ali.should be_an_instance_of(ProfileSequenceGlobalAlignmentLinearGap)
  end

  it "#global_alignment_linear_gap_cpp returns an instance of ProfileSequenceGlobalAlignmentLinearGap class" do
    ali = @psa.global_alignment_linear_gap_cpp
    ali.should be_an_instance_of(ProfileSequenceGlobalAlignmentLinearGap)
  end

  it "#global_alignment_affine_gap_cpp returns an instance of ProfileSequenceGlobalAlignmentAffineGap class" do
    ali = @psa.global_alignment_affine_gap
    ali.should be_an_instance_of(ProfileSequenceGlobalAlignmentAffineGap)
  end

  it "#global_alignment_affine_gap_rb returns an instance of ProfileSequenceGlobalAlignmentAffineGap class" do
    ali = @psa.global_alignment_affine_gap_rb
    ali.should be_an_instance_of(ProfileSequenceGlobalAlignmentAffineGap)
  end

  it "#global_alignment_affine_gap_cpp returns an instance of ProfileSequenceGlobalAlignmentAffineGap class" do
    ali = @psa.global_alignment_affine_gap_cpp
    ali.should be_an_instance_of(ProfileSequenceGlobalAlignmentAffineGap)
  end

end
