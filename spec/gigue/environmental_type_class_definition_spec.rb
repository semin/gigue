require File.join(File.dirname(__FILE__), '..', 'spec_helper.rb')
include Gigue

describe EnvironmentTypeClassDefinition do

  before(:all) do
    @file   = File.join(File.dirname(__FILE__), '..', 'classdef_std64.dat')
    @envdef = EnvironmentTypeClassDefinition.new(@file)
    @envs   = @envdef.environments
  end

  it "#file returns its source file location" do
    @envdef.file.should == @file
  end

  it "has environments" do
    @envdef.should have(5).environments
  end

  it "returns correct description of the first environment" do
    @envs[0].name.should == 'secondary structure and phi angle'
  end

  it "returns correct class values of the first environment" do
    @envs[0].values.should == %w[H E P C]
  end

  it "returns correct class labels of the first environment" do
    @envs[0].labels.should == %w[H E P C]
  end

  it "returns correct constraint status of the first environment" do
    @envs[0].constrained.should be_false
  end

  it "returns correct silency status of the first environment" do
    @envs[0].silent.should be_false
  end

  it "returns correct description of the last environment" do
    @envs[-1].name.should == 'hydrogen bond to mainchain NH'
  end

  it "returns correct class values of the last environment" do
    @envs[-1].values.should == %w[T F]
  end

  it "returns correct class labels of the last environment" do
    @envs[-1].labels.should == %w[N n]
  end

  it "returns correct constraint status of the last environment" do
    @envs[-1].constrained.should be_false
  end

  it "returns correct silency status of the last environment" do
    @envs[-1].silent.should be_false
  end

end
