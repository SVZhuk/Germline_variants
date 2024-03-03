#!/bin/bash

###### Request resources ######

#SBATCH --exclude=it03
#SBATCH --job-name=sarek-test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBARCH --cpus-per-task=2
#SBATCH --mem=6000
#SBATCH --partition=long
#SBATCH --time=160:00:00
#SBATCH --output=sarek-test.out
#SBATCH --mail-type=NONE
#SBATCH --mail-user=szhuk@ku.edu.tr

###############################

set -ex -o pipefail
ulimit -s unlimited
ulimit -l unlimited
ulimit -a

module load singularity/3.10.4
module load java/17.0.1

source activate nf-core

echo "=============================== SAREK PIPELINE =========================="

PipelineDir="/kuttam_fg/refdata/zhuk/pipelines/3_3_2/"
nf_config="/scratch/users/szhuk/hpc_run/workspace/scr/kuacc.config"
param_file="/kuacc/users/szhuk/hpc_run/workspace/Germline_variants/nf-params-kuacc.json"

nextflow run ${PipelineDir} \
    -profile singularity \
    -c ${nf_config} \
    -params-file ${param_file} \
    -resume
