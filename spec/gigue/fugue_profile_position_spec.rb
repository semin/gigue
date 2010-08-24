require_relative '../spec_helper.rb'
include Gigue

describe FugueProfilePosition do

  before(:all) do
    # Seq /            A    C    D    E    F    G    H    I    K    L    M    N    P    Q    R    S    T    V    W    Y    J    U / InsO InsE DelO DelE COIL HNcp HCcp HInt SNcp SCcp SInt NRes  Ooi  Acc /  H   E   P   C   A   a   S   s   O   o   N   n
    # -TTA-TTTT-TRQ    3  -79   -2    1  -22  -19   -7  -15    3  -20  -11   -2  -10    3    8   16   37  -12  -29  -12   -3  -12   100  100  100  100    0    0    0    0    0    0    0   87    1    7     0   0   0   0  87   0  27  60   0  87   0  87
    @probe  = '-TTA-TTTT-TRQ'
    @aas    = 'ACDEFGHIKLMNPQRSTVWYJU'.split('')
    @mat_s  = @aas.to_hash('3  -79   -2    1  -22  -19   -7  -15    3  -20  -11   -2  -10    3    8   16   37  -12  -29  -12   -3  -12'.split.map(&:to_i))
    @gaps   = 'InsO InsE DelO DelE COIL HNcp HCcp HInt SNcp SCcp SInt NRes  Ooi  Acc'.split(/\s+/)
    @gap_s  = @gaps.to_hash('100  100  100  100    0    0    0    0    0    0    0   87    1    7'.split.map(&:to_i))
    @envs   = 'H   E   P   C   A   a   S   s   O   o   N   n'.split(/\s+/)
    @env_s  = @envs.to_hash('0   0   0   0  87   0  27  60   0  87   0  87'.split.map(&:to_i))
    @pos    = FugueProfilePosition.new(@probe, @mat_s, @gap_s, @env_s)
  end

  it "#probe returns it probe string" do
    @pos.probe.should == @probe
  end

  it "#mat_score('A') returns a match score with 'A' for the position" do
    @pos.mat_score('A').should == 3
  end

  it "#gap_score(GAPID) returns a gap score of GAPID for the position" do
    @pos.gap_score('InsO').should == 100
  end

  it "#env_score(ENVID) returns a environment score of ENVID for the position" do
    @pos.env_score('n').should == 87
  end

  it "#gap_ins_open returns a insertion gap opening score for the position" do
    @pos.gap_ins_open.should == 100
  end

  it "#gap_ins_ext returns a insertion gap extenstion score for the position" do
    @pos.gap_ins_ext.should == 100
  end

  it "#gap_del_open returns a deletion gap opening score for the position" do
    @pos.gap_del_open.should == 100
  end

  it "#gap_del_ext returns a deletion gap extenstion score for the position" do
    @pos.gap_del_ext.should == 100
  end

end
