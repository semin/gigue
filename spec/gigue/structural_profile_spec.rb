require_relative '../spec_helper.rb'
include Gigue

describe StructuralProfile do

  before(:all) do
    @joy  = File.join(File.dirname(__FILE__), '..', 'd.240.1.1.tem')
    @esst = File.join(File.dirname(__FILE__), '..', 'ulla-mlr-logo-toccata-maskA.mat')
    @stp  = StructuralProfile.new(@joy, @esst)
  end

  it "#length returns the length of structural profile" do
    @stp.length.should == 157
  end

  it "#depth returns the number of structures in a profile" do
    @stp.depth.should == 5
  end

  it "can weights sequences using a VA (Vingron and Argos) weighting scheme" do
    va_stp = StructuralProfile.new(@joy, @esst, :weighting => :va)
    va_stp.sequences[0].weight.should == 0.20124031688419095
    va_stp.sequences[1].weight.should == 0.19963716993679606
    va_stp.sequences[2].weight.should == 0.2041590477774935
    va_stp.sequences[3].weight.should == 0.20117639883404523
    va_stp.sequences[4].weight.should == 0.19378706656747433
  end

  it "can weights sequences using a BLOSUM-like weighting scheme" do
    blosum_stp = StructuralProfile.new(@joy, @esst, :weighting => :blosum)
    blosum_stp.sequences[0].weight.should == 0.2
    blosum_stp.sequences[1].weight.should == 0.2
    blosum_stp.sequences[2].weight.should == 0.2
    blosum_stp.sequences[3].weight.should == 0.2
    blosum_stp.sequences[4].weight.should == 0.2
  end

end

