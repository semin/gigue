require 'rubygems'
require "rspec/core/rake_task"
gem 'hoe', '>= 2.1.0'
require 'hoe'
require 'fileutils'
require './lib/gigue'

Hoe.plugin :newgem
Hoe.plugin :website
# Hoe.plugin :cucumberfeatures

# Generate all the Rake tasks
# Run 'rake -T' to see list of generated tasks (from gem root directory)
$hoe = Hoe.spec 'gigue' do
  self.developer 'Semin Lee', 'seminlee@gmail.com'
  self.post_install_message = 'PostInstall.txt' # TODO remove if post-install message not required
  self.rubyforge_name       = self.name # TODO this is default value
  self.extra_deps           = [ ['narray','>= 0.5.9.7'],
                                ['RubyInline', '>=3.8.5'],
                                ['bio', '>=1.4.0'],
                                ['parallel', '>=0.4.6'] ]
end

require 'newgem/tasks'
#Dir['tasks/**/*.rake'].each { |t| load t }
#Dir['tasks/*.rake'].each { |t| load t }

# TODO - want other tests/tasks run by default? Add them to the list
#remove_task :default
#task :default => [:spec, :features]
desc "Run all specs"
RSpec::Core::RakeTask.new(:spec) do |t|
  t.pattern = 'spec/**/*_spec.rb'
end
remove_task :default
task :default => [:spec]
