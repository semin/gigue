require File.join(File.dirname(__FILE__), '..', 'spec_helper.rb')
include Gigue

describe SequenceProfilePosition do

  before(:all) do
    @spp = SequenceProfilePosition.new('--MMIIV-MVV-')
  end

  it "#frequency_of(A) returns observation count for amino acid A" do
    @spp.frequency_of('M').should == 3
    @spp.frequency_of('I').should == 2
    @spp.frequency_of('V').should == 3
    @spp.frequency_of('-').should == 4
  end

  it "#probability_of(A) returns observation probability of amino acid A" do
    size = @spp.probe.size.to_f
    @spp.probability_of('M').should == 3/size
    @spp.probability_of('I').should == 2/size
    @spp.probability_of('V').should == 3/size
    @spp.probability_of('-').should == 4/size
  end

end
