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

  it "#local_alignment_linear_gap returns an instance of ProfileProfileAlignmentLinearGap class" do
    ali = @ppa.local_alignment_linear_gap
    ali.should be_an_instance_of(ProfileProfileLocalAlignmentLinearGap)
  end

  it "#local_alignment_linear_gap_rb returns an instance of ProfileProfileAlignmentLinearGap class" do
    ali = @ppa.local_alignment_linear_gap_rb
    ali.should be_an_instance_of(ProfileProfileLocalAlignmentLinearGap)
  end

  it "#local_alignment_linear_gap_cpp returns an instance of ProfileProfileAlignmentLinearGap class" do
    ali = @ppa.local_alignment_linear_gap_cpp
    ali.should be_an_instance_of(ProfileProfileLocalAlignmentLinearGap)
  end

  it "both #local_alignment_linear_gap_rb and #local_alignment_linear_gap_cpp return the same raw_score" do
    ali_rb = @ppa.local_alignment_linear_gap_rb
    ali_cpp = @ppa.local_alignment_linear_gap_cpp
    ali_cpp.raw_score.should == ali_rb.raw_score
  end

  it "#local_alignment_affine_gap returns an instance of ProfileProfileAlignmentAffineGap class" do
    ali = @ppa.local_alignment_affine_gap
    ali.should be_an_instance_of(ProfileProfileLocalAlignmentAffineGap)
  end

  it "#local_alignment_affine_gap_rb returns an instance of ProfileProfileAlignmentLinearGap class" do
    ali = @ppa.local_alignment_affine_gap_rb
    ali.should be_an_instance_of(ProfileProfileLocalAlignmentAffineGap)
  end

  it "#local_alignment_affine_gap_cpp returns an instance of ProfileProfileAlignmentAffineGap class" do
    ali = @ppa.local_alignment_affine_gap_cpp
    ali.should be_an_instance_of(ProfileProfileLocalAlignmentAffineGap)
  end

  it "both #local_alignment_affine_gap_rb and #local_alignment_affine_gap_cpp return the same raw_score" do
    ali_rb = @ppa.local_alignment_affine_gap_rb
    ali_cpp = @ppa.local_alignment_affine_gap_cpp
    ali_cpp.raw_score.should == ali_rb.raw_score
  end

  it "#global_alignment_linear_gap returns an instance of ProfileProfileAlignmentLinearGap class" do
    ali = @ppa.global_alignment_linear_gap
    ali.should be_an_instance_of(ProfileProfileGlobalAlignmentLinearGap)
  end

  it "#global_alignment_linear_gap_rb returns an instance of ProfileProfileAlignmentLinearGap class" do
    ali = @ppa.global_alignment_linear_gap_rb
    ali.should be_an_instance_of(ProfileProfileGlobalAlignmentLinearGap)
  end

  it "#global_alignment_linear_gap_cpp returns an instance of ProfileProfileAlignmentLinearGap class" do
    ali = @ppa.global_alignment_linear_gap_cpp
    ali.should be_an_instance_of(ProfileProfileGlobalAlignmentLinearGap)
  end

  it "both #global_alignment_linear_gap_rb and #global_alignment_linear_gap_cpp return the same raw_score" do
    ali_rb = @ppa.global_alignment_linear_gap_rb
    ali_cpp = @ppa.global_alignment_linear_gap_cpp
    ali_cpp.raw_score.should == ali_rb.raw_score
  end

  it "#global_alignment_affine_gap returns an instance of ProfileProfileAlignmentAffineGap class" do
    ali = @ppa.global_alignment_affine_gap
    ali.should be_an_instance_of(ProfileProfileGlobalAlignmentAffineGap)
  end

  it "#global_alignment_affine_gap_rb returns an instance of ProfileProfileAlignmentAffineGap class" do
    ali = @ppa.global_alignment_affine_gap_rb
    ali.should be_an_instance_of(ProfileProfileGlobalAlignmentAffineGap)
  end

  it "#global_alignment_affine_gap_cpp returns an instance of ProfileProfileAlignmentAffineGap class" do
    ali = @ppa.global_alignment_affine_gap_cpp
    ali.should be_an_instance_of(ProfileProfileGlobalAlignmentAffineGap)
  end

  it "both #global_alignment_affine_gap_rb and #global_alignment_affine_gap_cpp return the same raw_score" do
    ali_rb = @ppa.global_alignment_affine_gap_rb
    ali_cpp = @ppa.global_alignment_affine_gap_cpp
    ali_cpp.raw_score.should == ali_rb.raw_score
  end

end
