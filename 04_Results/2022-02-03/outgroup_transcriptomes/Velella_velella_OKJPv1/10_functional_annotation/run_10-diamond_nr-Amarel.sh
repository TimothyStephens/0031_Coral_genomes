#!/bin/bash
#SBATCH --partition=main       	       	       	# Partition (job queue)
#SBATCH --requeue   	       	       	       	# Re-run job if preempted
#SBATCH --export=NONE  	       	       	       	# Export current environment variables to the launched application
#SBATCH --nodes=1      	       	       	       	# Number of nodes
#SBATCH --ntasks=1     	       	       	       	# Number of tasks (usually = cores) on each node
#SBATCH --cpus-per-task=52     	       	       	# Cores per task (>1 if multithread tasks)
#SBATCH --output=%x.slurm_out.%j-%2t.%a        	# STDOUT output file (will also contain STDERR if --error is not specified)
#SBATCH --mem=180G     	       	       	       	# Real memory (RAM) per node required (MB)
#SBATCH --time=72:00:00        	       	       	# Total run time limit (HH:MM:SS)
#SBATCH --job-name=run_10-diamond_nr-Amarel.sh	# Replace with your jobname


#### Pre-run setup
set -e -o pipefail
source ~/slurm_config/slurm_config
prerun_info # Print info
ulimit -n 10000 # Need to set high to avoid 'Too many open files' error
cd $SLURM_SUBMIT_DIR

export PATH="$PATH:/home/ts942/PROGRAMS/DIAMOND_v2.0.15/bin"
DB="/scratch/ts942/databases/ncbi/2022_07/nr.dmnd"
NCPUS=${SLURM_CPUS_PER_TASK}

PREFIX="Velella_velella_OKJPv1"
QUERY="${PREFIX}.transcripts.pep.faa"
OUT="${QUERY}.diamond_blastp_nr.outfmt6_short"

if [ ! -f "${OUT}.gz" ];
then
	run_cmd "md5sum ${QUERY} ../00_databases/${PREFIX}.transcripts.pep.faa | tee ${OUT}.job_md5sum_list.txt"
	run_cmd "diamond blastp --ultra-sensitive --max-target-seqs 0 --evalue 1e-5 --query ${QUERY} --db ${DB} --out ${OUT}.tmp --outfmt 6 qseqid sseqid evalue bitscore staxids --compress 1 --threads ${NCPUS} && mv ${OUT}.tmp.gz ${OUT}.gz"
else
	log "Already finished running this part!"
fi


