require_relative '../spec_helper.rb'
include Gigue

shared_examples_for "profile-profile alignment" do

  it "has a structural profile object" do
    @ali.structural_profile.should == @stp
    @ali.structural_profile.should be_an_instance_of(StructuralProfile)
  end

  it "has a sequence profile object" do
    @ali.sequence_profile.should == @sqp
    @ali.sequence_profile.should be_an_instance_of(SequenceProfile)
  end

  it "#raw_score returns raw alignment score as an integer" do
    @ali.raw_score.should be_a_kind_of(Integer)
  end

  #it "#calculate_z_score returns Z-score as a float" do
    #@ali.calculate_z_score(10).should be_a_kind_of(Float)
  #end

end

describe ProfileProfileAlignment do

  before(:all) do
    $logger.level = Logger::WARN

    tem   = File.join(File.dirname(__FILE__), '..', 'd.240.1.1.tem')
    sub   = File.join(File.dirname(__FILE__), '..', 'ulla-logo-smooth-toccata-maskA.mat')
    @stp  = StructuralProfile.create_from_joy_tem_and_essts(tem, sub)

    blt   = File.join(File.dirname(__FILE__), '..', 'test.blast.out6')
    msa   = MultipleSequenceAlignment.create_from_psiblast_output_style6(blt)
    @sqp  = msa.to_sequence_profile

    @ppa = ProfileProfileAligner.new(@stp, @sqp)
  end

  describe ProfileProfileLocalAlignmentLinearGap do
    before(:all) do
      @ali = @ppa.local_alignment_linear_gap
    end

    it_should_behave_like "profile-profile alignment"
  end

  describe ProfileProfileLocalAlignmentAffineGap do
    before(:all) do
      @ali = @ppa.local_alignment_affine_gap
    end

    it_should_behave_like "profile-profile alignment"
  end

  describe ProfileProfileGlobalAlignmentLinearGap do
    before(:all) do
      @ali = @ppa.global_alignment_linear_gap
    end

    it_should_behave_like "profile-profile alignment"
  end

  describe ProfileProfileGlobalAlignmentAffineGap do
    before(:all) do
      @ali = @ppa.global_alignment_affine_gap
    end

    it_should_behave_like "profile-profile alignment"
  end

end
