#!/usr/bin/env bash

###### Request resources ######

#SBATCH --exclude=it03
#SBATCH --job-name=sarek-test
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBARCH --cpus-per-task=2
#SBATCH --mem=6000
#SBATCH --partition=long
#SBATCH --time=160:00:00
#SBATCH --output=repeats.out
#SBATCH --mail-type=NONE
#SBATCH --mail-user=szhuk@ku.edu.tr

###############################

set -ex -o pipefail
ulimit -s unlimited
ulimit -l unlimited
ulimit -a

module load singularity/3.10.4
module load java/17.0.1


while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -i|--input_csv)
      INPUT_CSV="$2"
      shift
      shift
      ;;
    -o|--output_dir)
      OUTPUT_DIR="$2"
      shift
      shift
      ;;
    -h|--help)
      echo "Usage: run_exomes_sarek.sh -i|--input_csv <input_csv> -o|--output_dir <output_folder_name>"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

echo "Input spread sheet used : $INPUT_CSV"
echo "Name of output folder: $OUTPUT_DIR"

source activate nf-core

echo "=============================== SAREK PIPELINE =========================="

PipelineDir="/kuttam_fg/refdata/zhuk/pipelines/3_3_2/"
nf_config="/kuacc/users/szhuk/hpc_run/workspace/Germline_variants/configs/sarek_koc.config"
param_file="/kuacc/users/szhuk/hpc_run/workspace/Germline_variants/configs/nf-params-kuacc.json"

nextflow run ${PipelineDir} \
    -profile singularity \
    -c ${nf_config} \
    -params-file ${param_file} \
    --input ${INPUT_CSV} \
    --outdir results-${OUTPUT_DIR} \
    -resume
