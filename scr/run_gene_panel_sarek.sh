#!/usr/bin/env bash

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
    -p| --gene_panel)
      GENE_PANEL="$2"
      shift
      shift
      ;;
    -h|--help)
      echo "Usage: run_gene_panel_sarek.sh -i|--input_csv <input_csv> -o|--output_dir <output_folder_name> -p|--gene_panel <gene_panel_name> -h|--help"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

echo "Input file used : $INPUT_CSV"
echo "Name of output folder: $OUTPUT_DIR"
echo "Gene panel used: $GENE_PANEL"

if [[ $GENE_PANEL == "hyperchol" ]]; then
    PARAM_FILE="/mnt/hdd/nextflow_dir/workspace/scripts/Germline_variants/configs/hyperchol_almazov_params.json"
elif [[ $GENE_PANEL == "ccp" ]]; then
    PARAM_FILE="/mnt/hdd/nextflow_dir/workspace/scripts/Germline_variants/configs/ccp_almazov_params.json"
elif [[ $GENE_PANEL == "hcmp" ]]; then
    PARAM_FILE="/mnt/hdd/nextflow_dir/workspace/scripts/Germline_variants/configs/hcmp_almazov_params.json"
elif [[ $GENE_PANEL == "lqt" ]]; then
    PARAM_FILE="/mnt/hdd/nextflow_dir/workspace/scripts/Germline_variants/configs/lqt_almazov_params.json"
else
    echo "Unknown gene panel: $GENE_PANEL"
    echo "Please use one of the following gene panels: hyperchol, ccp, hcmp, lqt."
    exit 1
fi

CONFIG_FILE="/mnt/hdd/nextflow_dir/workspace/scripts/Germline_variants/configs/wes_almazov.config"

echo "Activating nf-core conda environment."
source /mnt/hdd/nextflow_dir/workspace/anaconda3/etc/profile.d/conda.sh
conda activate nf-core-new

docker system prune --all --force
docker image prune --all --force

nextflow run nf-core/sarek -r 3.4.0 \
    -profile docker \
    -c $CONFIG_FILE \
    -params-file $PARAM_FILE \
    --input $INPUT_CSV \
    --outdir results-$OUTPUT_DIR \
    -resume
