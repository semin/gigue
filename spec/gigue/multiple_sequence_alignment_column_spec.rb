require File.join(File.dirname(__FILE__), '..', 'spec_helper.rb')
include Gigue

describe MultipleSequenceAlignmentColumn do

  before(:all) do
    @probe  = '--AAABB-AB'
    @msacol = MultipleSequenceAlignmentColumn.new(@probe)
  end

  it "#probe returns its probe string" do
    @msacol.probe.should == @probe
  end

end
