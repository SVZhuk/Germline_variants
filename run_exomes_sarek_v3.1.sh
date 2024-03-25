#!/bin/bash

set -ex -o pipefail

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
      echo "Usage: run_exomes_sarek_v3.1.sh -i|--input_csv <input_csv> -o|--output_dir <output_folder_name>"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

PIPELINE_DIR="/home/students/nextflow_dir/pipelines/nf-core-sarek-3.1/workflow/"
PARAM_FILE="/mnt/hdd/nextflow_dir/workspace/scripts/Germline_variants/configs/wes_almazov_params.yaml"
CONFIG_FILE="/mnt/hdd/nextflow_dir/workspace/scripts/Germline_variants/configs/wes_almazov.config"

echo "Activating nf-core conda environment."
#source /mnt/hdd/nextflow_dir/workspace/anaconda3/etc/profile.d/conda.sh
conda activate nf-core

docker system prune --all --force
docker image prune --all --force

nextflow run ${PIPELINE_DIR}  \
    -profile docker \
    -c $CONFIG_FILE
    -params-file $PARAM_FILE \
    --input ${INPUT_CSV} \
    --outdir results-${OUTPUT_DIR} \
