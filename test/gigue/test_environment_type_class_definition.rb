require File.join(File.dirname(__FILE__), '..', 'test_helper.rb')

class TestEnvironmentTypeClassDefinition < Test::Unit::TestCase

  include Gigue

  def setup
    file  = File.join(File.dirname(__FILE__), '..', 'classdef_canonical_env5.dat')
    @def  = EnvironmentTypeClassDefinition.new(file)
    @envs = @def.environments
  end

  def test_environment_number
    assert_equal(5, @envs.size)
  end

  def test_first_environment
    assert_equal('secondary structure and phi angle', @envs[0].name)
    assert_equal(%w[H E P C], @envs[0].values)
    assert_equal(%w[H E P C], @envs[0].labels)
    assert(!@envs[0].constrained)
    assert(!@envs[0].silent)
  end

  def test_last_environment
    assert_equal('hydrogen bond to mainchain NH', @envs[-1].name)
    assert_equal(%w[T F], @envs[-1].values)
    assert_equal(%w[N n], @envs[-1].labels)
    assert(!@envs[-1].constrained)
    assert(!@envs[-1].silent)
  end

end
