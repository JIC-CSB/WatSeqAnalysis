require 'csv'
require 'zlib'
require 'fileutils'
require 'optparse'

options = {
	:folder => "."
}

opts = OptionParser.new do |o|
  o.banner = "Usage: #{File.basename($0)} [options]"
  o.on("-f", '--input_folder [DIR]', "Folder containging mat.csv.gz and df.csv.gz as produced by  the R package bio.tilling") do |arg|
  	options[:folder] = arg
  end
end

opts.parse!

mat_path ="#{ options[:folder] }/mat.csv.gz"
df_path ="#{ options[:folder] }/df.csv.gz"
cov_path="#{options[:folder]}/allMergedCoverages.tab.gz"

def compare_names(left, right) 
	
	left.each_with_index do |l, i|
		if l != right[i]
			$stderr.puts "#{l} != #{right[i]}"
			return false
		end
	end
	return true
end

def each_region(mat_path, df_path, cov_path)

	puts mat_path
	puts df_path
	puts cov_path

	df_io = Zlib::GzipReader.open(df_path)
	mat_io = Zlib::GzipReader.open(mat_path)
	cov_io = Zlib::GzipReader.open(cov_path)

	df_each = df_io.each
	mat_each = mat_io.each
	cov_each = cov_io.each

	names = mat_each.next.gsub("\"","").gsub("merge\.rmdup\.", "").split(",").drop(1)
	df_names = df_each.next.gsub("\"","").split(",")
	cov_names = cov_each.next.gsub("merge\.rmdup\.", "").gsub("-",".").split("\t").drop(1)


	raise Exception.new "Columns not in the same order" unless compare_names(names, cov_names) 
	i = 1
	while true
		begin
			df_line = df_each.next.gsub("\"","").split(",")
			mat_line = mat_each.next.gsub("\"","").split(",").drop(1)
			cov_line = cov_each.next.split("\t").drop(1)
			yield names, mat_line, df_line, cov_line, i
			i += 1
			#break if i > 10
		rescue StopIteration
	  		$stderr.puts "Done iteration"
	  		break
	  	end
	end
	
  	puts "EXITING REGION"
  	cov_io.close
	mat_io.close
	df_io.close
end


out_path = "#{options[:folder]}/long_tables/"
FileUtils.mkdir_p out_path
last_names = []

df_out = File.open("#{out_path}/windows.bed","w")
df_out.puts ["chrom", "chromStart", "chromEnd", "reg_id", "sd"].join("\t")

norm_covs_out = File.open("#{out_path}/norm_covs.tsv", "w")
norm_covs_out.puts ["reg_id", "line_id", "norm_cov"].join("\t")

covs_out = File.open("#{out_path}/covs.tsv", "w")
covs_out.puts ["reg_id", "line_id", "norm_cov"].join("\t")

each_region(mat_path, df_path, cov_path) do |names, mat, df, cov, i|
	df_out.puts [df[1], df[2], df[3], i, df[6]].join("\t") 

	cov.each_with_index do |val, j|
		covs_out.puts [i, j+1, val].join("\t")
	end

	mat.each_with_index do |val, j|
		norm_covs_out.puts [i, j+1, val].join("\t")
	end

	#puts names.first(5).join("\t")
	#puts mat.first(5).join("\t")
	last_names = names
	#break if i > 20
end

covs_out.close
df_out.close
norm_covs_out.close

lines_out = File.open("#{out_path}/lines.tsv", "w")
lines_out.puts ["line_id", "line"].join("\t")
last_names.each_with_index do |l,i|
	lines_out.puts "#{i+1}\t#{l}"
end
lines_out.close


