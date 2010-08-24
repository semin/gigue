require_relative '../spec_helper.rb'
include Gigue

describe Esst do

  before(:all) do
    matrix = NMatrix[
      [4, -4, -1, -1, -1,  0, -2, 0, -1,  0,  0, -1,  0,  0, -1,  1, -1,  1, -2, -2,  3],
      [9, 23, -7, -7, -6,-10,-14,-8,-10, -8, -8, -5,-10, -9, -9,-15, -3, -6, -8, -1, -2],
      [1, -8,  9,  3, -5,  0, -2,-5,  0, -5, -4,  1,  0,  1, -1,  1,  2, -4, -5, -7, -3],
      [1, -7,  3,  8, -3, -1,  1,-2,  2, -2, -2,  1,  1,  3,  1,  1,  0, -2, -4, -4, -5],
      [3, -4, -8, -5,  7, -5, -2, 0, -4,  1,  1, -5, -5, -2, -3, -5,  1,  0,  3,  5, -1],
      [1, -7, -5, -4, -5,  8, -2,-6, -3, -5, -4, -3, -2, -1, -3, -1, -2, -4, -6, -8,  1],
      [1, -7, -1,  0,  0, -2,  9,-2,  0, -2, -1,  1, -2,  1,  1, -1,  0, -1, -1,  5, -2],
      [2, -7, -8, -3,  0, -6, -6, 4, -4,  2,  1, -6, -4, -4, -3, -2, -4,  2, -2, -2,  0],
      [1, -6,  0, -1, -3, -1,  1,-2,  6, -2, -2,  1,  0, -2,  3, -1,  0, -1, -4, -3, -3],
      [1, -6, -7, -4,  2, -5, -3, 3, -3,  4,  2, -4, -4, -2, -3, -4, -3,  1,  0, -2,  0],
      [1, -6,-10, -4,  2, -4,  1, 2, -2,  3,  7, -2, -4, -2, -1,  0, -3,  1,  0, -3, -2],
      [1, -5,  2,  0, -4,  0,  3,-4,  0, -3, -3,  9, -1, -2,  0,  0,  1, -3, -4, -6, -1],
      [1, -7, -6, -5, -5, -2, -7,-4, -2, -4, -4, -3,  9, -6, -2, -5, -5, -3, -5, -5, -4],
      [1, -6,  4,  2, -2, -1,  2,-2,  2, -1, -1,  1,  0,  9,  2,  2,  2, -1, -3,  0, -1],
      [0, -6, -1,  0, -2, -2,  1,-2,  3, -1, -1,  1, -1,  0,  7, -2, -1, -1, -3,  0, -1],
      [1, -4, -2, -2, -3,  0,  1,-3,  0, -3, -2,  0,  0, -1,  0,  5,  3, -2, -4, -7,  0],
      [1, -6, -4, -3, -2, -2, -1,-1, -1, -2, -1, -1, -1, -4, -1,  2,  5,  0, -4, -6, -1],
      [2, -5, -7, -6,  0, -5, -3, 2, -3,  0,  1, -4, -4, -4, -3,  0, -3,  3, -1, -4,  0],
      [3, -5, -6, -4,  3, -6, -1,-1, -3,  0,  0, -6, -4, -1, -2, -3, -5, -1, 14,  7, -5],
      [2, -3, -6, -3,  4, -4,  4, 0, -2,  0,  0, -2, -4, -3, -1, -5,  0, -1,  4, 10, -3],
      [0,  2, -6, -5,  0, -3,  0, 0, -3,  0,  1, -3, -5, -5, -3, -4,  2,  1, -1, -2, 10],
      [1, 17, -6, -5, -1, -4, -1,-1, -4, -1,  0, -3, -6, -6, -4, -6,  1,  0, -2, -2,  8]
    ]
    colnames  = 'ACDEFGHIKLMNPQRSTVWYJ'.split('')
    rownames  = 'ACDEFGHIKLMNPQRSTVWYJU'.split('')
    @esst     = Esst.new(:logo, 'HASON', 0, colnames, rownames, matrix)
  end

  it "#type returns its type" do
    @esst.type.should == :logo
  end

  it "#label returns label of environment class combination" do
    @esst.label.should == 'HASON'
  end

  it "#no returns number assigned to it" do
    @esst.no.should == 0
  end

  it "#colnames returns column names" do
    @esst.colnames.should == 'ACDEFGHIKLMNPQRSTVWYJ'.split('')
  end

  it "#rownames returns row names" do
    @esst.rownames.should == 'ACDEFGHIKLMNPQRSTVWYJU'.split('')
  end

  it "#matrix returns substitution matrix as an instance of NMatrix class" do
    @esst.matrix.should be_an_instance_of(NMatrix)
  end

  it "#scores_from('A') returns substitution scores from amino acid 'A'" do
    @esst.scores_from('A').should == NMatrix[[4], [9], [1], [1], [3], [1], [1], [2], [1], [1], [1], [1], [1], [1], [0], [1], [1], [2], [3], [2], [0], [1]]
  end

  it "#scores_to('A') returns substitution scores to amino acid 'A'" do
    @esst.scores_to('A').should == NMatrix[[4, -4, -1, -1, -1, 0, -2, 0, -1, 0, 0, -1, 0, 0, -1, 1, -1, 1, -2, -2, 3]]
  end

  it "#score('A', 'A') returns substitution score between amino acid 'A' and 'A'" do
    @esst.score('A', 'A').should == 4
  end

end

