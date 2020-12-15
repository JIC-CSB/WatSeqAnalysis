require 'csv'
require 'zlib'
require 'fileutils'
require 'optparse'

class Region

	attr_accessor :chromosome, :start, :end, :orientation

	def initialize(chromosome , start, finish, orientation: nil )
		@chromosome = chromosome
		@start = start
		@end = finish
		@orientation = orientation
	end

	def length
		self.end - self.start
	end

	def to_r
		"#{@chromosome}:#{@start}-#{@end}"
	end

	def overlap(other)
		ret = other.chromosome == self.chromosome
		ret &= (other.start.between?(self.start, self.end) or other.end.between?(self.start, self.end) )
		ret
	end

	def <=>(other)
		return self.chromosome <=> other.chromosome if self.chromosome != other.chromosome 
		return self.start <=> other.start if self.start != other.start 
		return self.end <=> other.end
	end

	def to_bed
		[
			self.chromosome, 
			self.start, 
			self.end, 
			self.orientation
		].compact.join("\t")
	end

	def contains(other)
		return false if other.chromosome != self.chromosome
		return other.start >= self.start && other.end <= self.end
	end

	def contains_all(others)
		others.each {|e| return false unless self.contains(e) }
		return true
	end

	def in_range(start, finish)
		left      = self.start <= finish && self.end >= finish 
		right     = self.start <= finish && self.end >= finish 
		contained = self.start >= start  && self.end <= finish
		return  left || right || contained; 
	end

	def self.bed(line) 
		bed_arr = line.split("\t")
		Region.new(bed_arr[0],bed_arr[1],bed_arr[2])
	end

	def self.mat(line)
		#chr1A_part1:68238:91244
		line_arr = line.chomp.split("\t")
		line_arr = line_arr[0].split(":")
		Region.new(line_arr[0],line_arr[1],line_arr[2])
	end
end


options = {
	:folder => "."
	:input => "allMergedCoverages.tab from the scripts to generate the coverage from the bed files"
}

opts = OptionParser.new do |o|
	o.banner = "Usage: #{File.basename($0)} [options]"
	o.on("-i", '--input [FILE]', " 'allMergedCoverages.tab'") do |arg|
		options[:input] = arg
	end
	o.on('-b', '--bed [FILE]', "BED file with the regions to keep") do |arg|
		options[:bed] = arg
	end
end

opts.parse!

cov_path="#{options[:input]}"
bed_path="#{options[:bed]}"

cov_io = Zlib::GzipReader.open(cov_path)
cov_each = cov_io.each 

bed_io = File.open(bed_path)
bed_each = bed_io.each
bed = Region.bed(bed_each.next)

out = $stdout
mat_line = mat_each.next
out.puts mat_line


mat_line = mat_each.next
mat = Region.mat(mat_line)
i = 0
while true
	begin
		while mat < bed #Advance mat table
			mat_line = mat_each.next
			mat = Region.mat(mat_line)
		end

		while bed.contains(mat)
			out.puts mat_line
			mat_line = mat_each.next
			mat = Region.mat(mat_line)
			i += 1
		end
		bed = Region.bed(bed_each.next)
		break if i > 10
	rescue StopIteration => e
		$stderr.puts "Iteration done"
	end
	
end

cov_io.close
bed_io.close
