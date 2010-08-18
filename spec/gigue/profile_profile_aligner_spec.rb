require File.join(File.dirname(__FILE__), '..', 'spec_helper.rb')
include Gigue

describe ProfileProfileAligner do

  before(:all) do
    $logger.level = Logger::WARN

    # build structural profile
    tem = File.join(File.dirname(__FILE__), '..', 'd.240.1.1.tem')
    sub = File.join(File.dirname(__FILE__), '..', 'ulla-logo-smooth-toccata-maskA.mat')
    @stp = StructuralProfile.new(tem, sub)

    # build sequence profile
    blt = File.join(File.dirname(__FILE__), '..', 'test.blast.out6')
    msa = MultipleSequenceAlignment.create_from_psiblast_output_style6(blt)
    @sqp = msa.to_sequence_profile

    @ppa = ProfileProfileAligner.new(@stp, @sqp)
  end

  it "has a structural profile object" do
    @ppa.structural_profile.should == @stp
  end

  it "has a sequence profile object" do
    @ppa.sequence_profile.should == @sqp
  end

  it "#local_alignment_with_linear_gap_penalty returns an instance of ProfileProfileAlignmentLinearGap class" do
    ali = @ppa.local_alignment_with_linear_gap_penalty
    ali.should be_an_instance_of(ProfileProfileLocalAlignmentLinearGap)
  end

  it "#local_alignment_with_affine_gap_penalty returns an instance of ProfileProfileAlignmentAffineGap class" do
    ali = @ppa.local_alignment_with_affine_gap_penalty
    ali.should be_an_instance_of(ProfileProfileLocalAlignmentAffineGap)
  end

  it "#global_alignment_with_linear_gap_penalty returns an instance of ProfileProfileAlignmentLinearGap class" do
    ali = @ppa.global_alignment_with_linear_gap_penalty
    ali.should be_an_instance_of(ProfileProfileGlobalAlignmentLinearGap)
  end

  it "#global_alignment_with_affine_gap_penalty returns an instance of ProfileProfileAlignmentAffineGap class" do
    ali = @ppa.global_alignment_with_affine_gap_penalty
    ali.should be_an_instance_of(ProfileProfileGlobalAlignmentAffineGap)
  end

end
