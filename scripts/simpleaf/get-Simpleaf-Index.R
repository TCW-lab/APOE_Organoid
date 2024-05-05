#!/bin/bash
#PBS -l walltime=2:00:00   # Adjust the time limit as needed
#PBS -P tcwlab    	#specified project name
#PBS -l nodes=1:ppn=12      # Specify the number of nodes and CPU cores
#PBS -N "Simpleaf_Alevin_Fry"        # job name
#PBS -o job_output.log  # output log file
#PBS -e job_error.log   # error log file

# Load Conda environment
source /usr3/graduate/akg/miniconda3/bin/activate /usr3/graduate/akg/miniconda3/envs/af

# Set environment variable
export ALEVIN_FRY_HOME="$PWD"

# Set ulimit
ulimit -n 4096

# Run commands
simpleaf set-paths
simpleaf index --output /projectnb/tcwlab/LabMember/akg/scRNASeq/index_dir \
  --fasta /projectnb/tcwlab/RawData/Genome/hg38/Homo_sapiens_assembly38.fasta \
  --gtf /projectnb/tcwlab/RawData/Genome/hg38/gencode.v26.annotation.gtf