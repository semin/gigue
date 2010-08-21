require_relative '../spec_helper.rb'
include Gigue

describe JoyTem do

  before(:all) do
    @file   = File.join(File.dirname(__FILE__), '..', 'd.240.1.1.tem')
    @joytem = JoyTem.new(@file)
  end

  it "#file returns source file location" do
    @joytem.file.should == @file
  end

  it "has correct information for the fist entry (structure)" do
    code  = '106363'
    desc  = 'sequence'
    data  = '-Q----SFSEED-SFK-KC-SSEV--EAKNKIEELLASLLNRV-CQDG----R-K--PHTVRLIIRR--YSSEKHYGRESRQ--C--PIPSHVIQKLGTGNY--DVMTPM-VDILMKLF-RNMVNV-KMPF--HLTLLSVCFCNLK-----------'
    @joytem.entries.has_key?(code).should be_true
    @joytem.entries[code].has_key?(desc).should be_true
    @joytem.entries[code][desc].should == data
  end

  it "has correct information for the last entry (structure)" do
    code  = '106701'
    desc  = 'Ooi number'
    data  = '----012343433-2-2-233--3322443453344333-2--221----1333----33445553322------22233333--2---2221----10----4--5434543454343-2--21111-24--34345555432-------------'
    @joytem.entries.has_key?(code).should be_true
    @joytem.entries[code].has_key?(desc).should be_true
    @joytem.entries[code][desc].should == data
  end

end
