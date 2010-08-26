require_relative '../spec_helper.rb'
include Gigue

shared_examples_for "profile-sequence alignment" do

  it "has a structural profile object" do
    @ali.structural_profile.should == @prf
    @ali.structural_profile.should be_an_instance_of(StructuralProfile)
  end

  it "has a sequence" do
    @ali.sequence.should == @seq
    @ali.sequence.should be_an_instance_of(Sequence)
  end

  it "#raw_score returns raw alignment score as an integer" do
    @ali.raw_score.should be_a_kind_of(Integer)
  end

  it "#calculate_z_score returns Z-score as a float" do
    @ali.calculate_z_score(5).should be_a_kind_of(Float)
  end

end

describe ProfileSequenceAlignment do

  before(:all) do
    $logger.level = Logger::WARN

    @tem = File.join(File.dirname(__FILE__), "..", "d.240.1.1.tem")
    @sub = File.join(File.dirname(__FILE__), "..", "ulla-logo-smooth-toccata-maskA.mat")
    @prf = StructuralProfile.create_from_joy_tem_and_essts(@tem, @sub)
    @seq = Sequence.new('VRKSIGRIVTMKRNSRNLEEIKPYLFRAIEESYYKLDKRIPKAIHVVAVTEDLDIVSRGRTFPHGISKETAYSESVKLLQKILEEDERKIRRIGVRFSKFI', 'test')
    @psa = ProfileSequenceAligner.new(@prf, @seq)
  end

  describe ProfileSequenceLocalAlignmentLinearGap do
    before(:all) do
      @ali = @psa.local_alignment_linear_gap
    end

    it_should_behave_like "profile-sequence alignment"

  end

  describe ProfileSequenceLocalAlignmentAffineGap do
    before(:all) do
      @ali = @psa.local_alignment_affine_gap
    end

    it_should_behave_like "profile-sequence alignment"
  end

  describe ProfileSequenceGlobalAlignmentLinearGap do
    before(:all) do
      @ali = @psa.global_alignment_linear_gap
    end

    it_should_behave_like "profile-sequence alignment"
  end

  describe ProfileSequenceGlobalAlignmentAffineGap do
    before(:all) do
      @ali = @psa.global_alignment_affine_gap
    end

    it_should_behave_like "profile-sequence alignment"
  end

end
