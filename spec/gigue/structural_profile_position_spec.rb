require File.join(File.dirname(__FILE__), '..', 'spec_helper.rb')
include Gigue

describe StructuralProfilePosition do

  before(:all) do
    @probe  = '-TTA-TTTT-TRQ'
    @aas    = 'ACDEFGHIKLMNPQRSTVWYJU'.split('')
    @mat_s  = @aas.to_hash([3, -79, -2, 1, -22, -19, -7, -15, 3, -20, -11, -2, -10, 3, 8, 16, 37, -12, -29, -12, -3, -12])
    @gaps   = %w[InsO InsE DelO DelE]
    @gap_s  = @gaps.to_hash([100, 90, 80, 70])
    @pos    = StructuralProfilePosition.new(@probe, @mat_s, @gap_s)
  end

  it "#probe returns it probe string" do
    @pos.probe.should == @probe
  end

  it "#mat_score('C') returns a match score with 'C' for the position" do
    @pos.mat_score('C').should == -79
  end

  it "#gap_score(GAPID) returns a gap score of GAPID for the position" do
    @pos.gap_score('InsO').should == 100
  end

  it "#gap_ins_open returns a insertion gap opening score for the position" do
    @pos.gap_ins_open.should == 100
  end

  it "#gap_ins_ext returns a insertion gap extenstion score for the position" do
    @pos.gap_ins_ext.should == 90
  end

  it "#gap_del_open returns a deletion gap opening score for the position" do
    @pos.gap_del_open.should == 80
  end

  it "#gap_del_ext returns a deletion gap extenstion score for the position" do
    @pos.gap_del_ext.should == 70
  end

end
