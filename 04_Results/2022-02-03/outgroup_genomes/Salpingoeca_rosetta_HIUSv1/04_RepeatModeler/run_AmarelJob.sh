#!/bin/bash
#SBATCH --partition=main 			# Partition (job queue) 
#SBATCH --no-requeue				# Do not re-run job if preempted
#SBATCH --export=ALL 				# Export current environment variables to the launched application
#SBATCH --nodes=1 				# Number of nodes 
#SBATCH --ntasks=1 				# Number of tasks (usually = cores) on each node
#SBATCH --cpus-per-task=24 			# Cores per task (>1 if multithread tasks)
#SBATCH --output=%x.slurm_out.%j-%2t.%a 	# STDOUT output file (will also contain STDERR if --error is not specified)
#SBATCH --mem=120G 				# Real memory (RAM) per node required (MB) 
#SBATCH --time=72:00:00 			# Total run time limit (HH:MM:SS) 
#SBATCH --job-name=AmarelJob 			# Replace with your jobname


#### Pre-run setup
source ~/slurm_config/slurm_config_v0.4.sh
set +eu; conda activate py27; set -eu


#### Start Script
RUNSCRIPT=$(ls -1 run_Repeat*.sh)
CHECKPOINT="${RUNSCRIPT}.checkpoint.done"
bash "${RUNSCRIPT}" \
  && (touch "${CHECKPOINT}" && log "  - Done!") \
  || log "ERROR:${RUNSCRIPT} failed to finish correctly!"  


