require File.join(File.dirname(__FILE__), '..', 'spec_helper.rb')
include Gigue

describe Essts do

  before(:all) do
    @file   = File.join(File.dirname(__FILE__), '..', 'subst-logo-smooth-toccata-maskA.mat')
    @essts  = Essts.new(@file)
  end

  it "#file returns its source file location" do
    @essts.file.should == @file
  end

  it "#type return its type" do
    @essts.type.should == :logo
  end

  it "has environments" do
    @essts.should have(6).environments
  end

  it "has environment-specific substition tables" do
    @essts.should have(64).essts
  end

  it "#colnames returns column names of ESSTs" do
    @essts.colnames.should == @essts[0].colnames
  end

  it "#rownames returns row names of ESSTs" do
    @essts.rownames.should == @essts[0].rownames
  end

  it "#[i] returns ith ESST if i is integer" do
    @essts[0].label.should == 'HASON'
  end

  it "#[l] returns ESST having label, l as String" do
    @essts['HASON'].no.should     == 0
    @essts['HASON'].label.should  == 'HASON'
  end

end
