#!/bin/bash -e
#SBATCH -p nbi-largemem
#SBATCH --mem=700G
#SBATCH --mail-type=END,FAIL 
#SBATCH -N 1 
#SBATCH -J dels
#SBATCH --array=1-3
#SBATCH -o log/norm_cov_v3_\%A_\%a.out
#SBATCH --time=2-00:00:00 

function log_line() {
	echo $(date) "$1" >&2
}

i=$SLURM_ARRAY_TASK_ID

declare -a sizes=("100bp" "150bp" "200bp")


size=$(($i%3))

size=${sizes[$size]}

log_line "Starting $i: $size"

source bio.tilling-20180922
dir="/jic/scratch/projects/watseq/window_coverage"

Rscript scripts/detect_deletions_simple.R -f "$dir/bedCov/window_$size/allMergedCoverages.tab.gz"  -z  -o "$dir/bedCov/deletions_20200702/$size/"
#Rscript scripts/detect_deletions_simple_v2.R -f "$dir/bedCov/window_$size/allMergedCoverages.tab.gz"  -z  -o "$dir/bedCov/deletions_20200623/$size/"


log_line "DONE $size"
