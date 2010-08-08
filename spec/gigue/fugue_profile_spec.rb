require File.join(File.dirname(__FILE__), '..', 'spec_helper.rb')
include Gigue

describe FugueProfile do

  before(:all) do
    file      = File.join(File.dirname(__FILE__), '..', 'd.240.1.1.fug')
    @profile  = FugueProfile.new(file)
  end

  it "has sequences" do
    @profile.should have(5).sequences
  end

  it "stores sequences with correct codes" do
    @profile.sequences[0].code == '106363'
    @profile.sequences[1].code == '90374'
    @profile.sequences[2].code == '99680'
    @profile.sequences[3].code == '90378'
    @profile.sequences[4].code == '106701'
  end

  it "stores weights for sequences" do
    @profile.weights[@profile.sequences[0].code] == 0.2
    @profile.weights[@profile.sequences[1].code] == 0.2
    @profile.weights[@profile.sequences[2].code] == 0.2
    @profile.weights[@profile.sequences[3].code] == 0.2
    @profile.weights[@profile.sequences[4].code] == 0.2
  end

  it "has positions" do
    @profile.should have(157).positions
  end

  it "has correct information for the first position" do
    @profile.positions[0].probe.should              == '-P---'
    @profile.positions[0].mat_score('A').should     == -20
    @profile.positions[0].mat_score('U').should     == -50
    @profile.positions[0].gap_score('InsO').should  == 100
    @profile.positions[0].gap_score('Acc').should   == 11
    @profile.positions[0].env_score('H').should     == 0
    @profile.positions[0].env_score('n').should     == 100
  end

  it "has correct information for the last position" do
    @profile.positions[-1].probe.should              == '--L--'
    @profile.positions[-1].mat_score('A').should     == -20
    @profile.positions[-1].mat_score('U').should     == -20
    @profile.positions[-1].gap_score('InsO').should  == 100
    @profile.positions[-1].gap_score('Acc').should   == 11
    @profile.positions[-1].env_score('H').should     == 0
    @profile.positions[-1].env_score('n').should     == 100
  end

  it "#command returns program and its arguments used to generate fugue profile" do
    @profile.command.should == "melody -t d.240.1.1.tem -c ../classdef_canonical_env5.dat -s ../toccata/MaskA.dat"
  end

  it "#length returns profile length" do
    @profile.length.should == 157
  end

  it "#weighting returns weighting scheme" do
    @profile.weighting.should == 1
  end

  it "#multiple_factor returns multiplication factor" do
    @profile.multiple_factor.should == 10.0
  end

  it "#format returns profile type" do
    @profile.format.should == '0-FUGUE'
  end

  it "#rowsymbols returns symbols used for profile rows" do
    @profile.rowsymbols.should == 'ACDEFGHIKLMNPQRSTVWYJU'.split('')
  end

  it "#colsymbols returns symbols used for profile columns" do
    @profile.colsymbols.should == 'ACDEFGHIKLMNPQRSTVWYJ'.split('')
  end

  it "#envsymbols returns symbols used for environments class labels" do
    @profile.envsymbols.should == 'HEPCAaSsOoNn'.split('')
  end

  it "#gap_open_ins_term returns terminal insertion gap opening penalty" do
    @profile.gap_open_ins_term == 100
  end

  it "#gap_open_del_term returns terminal deletion gap opening penalty" do
    @profile.gap_open_del_term == 100
  end

  it "#gap_ext_ins_term returns terminal insertion gap extension penalty" do
    @profile.gap_ext_ins_term == 100
  end

  it "#gap_ext_del_term returns terminal deletion gap extension penalty" do
    @profile.gap_ext_del_term == 100
  end

  it "#aa_colnames returns column names of amino acids" do
    @profile.aa_colnames.should == 'ACDEFGHIKLMNPQRSTVWYJU'.split('')
  end

  it "#gap_colnames returns column names of gaps" do
    @profile.gap_colnames.should == 'InsO InsE DelO DelE COIL HNcp HCcp HInt SNcp SCcp SInt NRes  Ooi  Acc'.split
  end

  it "#env_colnamse returns column names of environments" do
    @profile.env_colnames.should == 'H   E   P   C   A   a   S   s   O   o   N   n'.split
  end

end
