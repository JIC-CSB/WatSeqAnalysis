#!/bin/bash
#SBATCH -p RG-Cristobal-Uauy,jic-medium,nbi-medium
#SBATCH -c 1
#SBATCH -N 1
#SBATCH --mem 16000
#SBATCH -o log/merge_cov_\%A_\%a.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ricardo.ramirez-gonzalez@jic.ac.uk
#SBATCH --array=100,150,200
#SBATCH --time=2-00:00:00 
#SBATCH -J merge_cov

function log_line() {
	echo $(date) "$1" >&2
}

function merge_bed_cov_gz(){
	
	folder="$1"
	cd $folder

	for fullfile in `ls *.bedCov.gz`; do
	    filename=$(basename "$fullfile")
	    extension="${filename##*.}"
	    filename="${filename%.*.gz}"
	    log_line "Extracting $filename"
	    echo $filename > cov.${filename}.txt
	    gunzip -c $fullfile | awk '{print $4}'  >> cov.${filename}.txt
	    break
	done


	echo "" > exonNames.txt
	gunzip -c $fullfile | awk '{print $1":"$2":"$3}' >> exonNames.txt
	#head exonNames.txt
	#This line may be needed if your user can't have many file opens. This is because to paste all the lines we have all the bed files open at the same time. 
	ulimit -S -n 2048
	paste exonNames.txt cov*.txt | gzip -c >  allMergedCoverages.tab.gz
}

bp=$SLURM_ARRAY_TASK_ID
#bp=100
log_line "Starting for ${bp}bp"
merge_bed_cov_gz ./bedCov/window_${bp}bp

log_line "DONE"