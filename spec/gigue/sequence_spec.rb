require File.join(File.dirname(__FILE__), '..', 'spec_helper.rb')
include Gigue

describe Sequence do

  before(:all) do
    @code = 'test'
    @desc = 'sequence'
    @data = 'VRKSIGRIVTMKRNSRNLEEIKPYLFRAIEESYYKLDKRIPKAIHVVAVTEDLDIVSRGRTFPHGISKETAYSESVKLLQKILEEDERKIRRIGVRFSKFI'
    @seq = Sequence.new(@data, @code, @desc)
  end

  it "#data returns sequence's data" do
    @seq.data.should == @data
  end

  it "#code returns sequence's code" do
    @seq.code.should == @code
  end

  it "#description returns sequence's description" do
    @seq.description.should == @desc
  end

  it "#amino_acids returns an array of amino acid characters constituting a sequence" do
    @seq.amino_acids.should == @data.split('')
  end

  it "#length returns the length of a sequence" do
    @seq.length.should == @data.length
  end

  it "#[i] returns an amino acid character at a position, i in a sequence" do
    @seq[0].should == @data[0]
  end

  it "#gapless returns a new Sequence object with gapless data" do
    data  = '---ABCDEFG--AABB--'
    code  = 'seqgap'
    desc  = 'sequence'
    seqg  = Sequence.new(data, code, desc)
    seqgl = seqg.gapless
    seqgl.should be_an_instance_of(Sequence)
    seqgl.data.should == data.gsub('-', '')
    seqgl.code.should == code
    seqgl.description.should == desc
  end

  it "#gapless! returns a Sequence object itself with gapless data" do
    data  = '---ABCDEFG--AABB--'
    code  = 'seqgap'
    desc  = 'sequence'
    seqg  = Sequence.new(data, code, desc)
    seqgl = seqg.gapless!
    seqgl.should be_an_instance_of(Sequence)
    seqgl.data.should == data.gsub('-', '')
    seqgl.code.should == code
    seqgl.description.should == desc
    seqgl.object_id.should == seqg.object_id
  end

  it "#shuffle_c returns a Sequnece object with shuffled data" do
    shuffled_seq = @seq.shuffle_c
    shuffled_seq.should be_an_instance_of(Sequence)
    shuffled_seq.data.should_not == @data
    shuffled_seq.amino_acids.should_not == @data.split('')
  end

end
