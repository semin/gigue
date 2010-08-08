require File.join(File.dirname(__FILE__), '..', 'spec_helper.rb')
include Gigue

describe StructuralProfile do

  before(:all) do
    @joy  = File.join(File.dirname(__FILE__), '..', 'd.240.1.1.tem')
    @esst = File.join(File.dirname(__FILE__), '..', 'ulla-mlr-logo-toccata-maskA.mat')
    @prof = StructuralProfile.new(@joy, @esst)
  end

  it "#length returns the length of structural profile" do
    @prof.length.should == 157
  end

  it "weights sequences using BLOSUM like weighting scheme" do
    i = 0
    @prof.weights.each do |k, v|
      k.should == @prof.sequences[i].code
      v.should == 0.2
      i += 1
    end
  end

end

