#!/usr/bin/env bash

set -o pipefail
set -o errexit

SELF_PATH=$(dirname "${BASH_SOURCE[0]}")

source "${SELF_PATH}/common/job.sh"
source "${SELF_PATH}/preparation/vcf_preparation.sh"

function usage() {
  echo "Usage: $0 -i|--input_csv <input_csv> -o|--output_dir <output_folder_name> -p|--gene_panel <gene_panel_name> -h|--help"
  exit 1
}

function source_env() {
  if [[ -f pipeline.env ]]; then
    source pipeline.env
  else
    echo "ERROR: Missing pipeline.env file"
    exit 1
  fi

  required_vars=("EXOME_INTERVALS" "HYPERCHOL_INTERVALS" "CCP_INTERVALS" "HCMP_INTERVALS" "LQT_INTERVALS" "CONDA_ENV_PATH" "CONDA_ENV_NAME" "CONFIG_FILE" "PARAM_FILE")
  for var in "${required_vars[@]}"; do
    if [[ -z "${!var:-}" ]]; then
      echo "ERROR: Environment variable $var is not set"
      exit 1
    fi
  done
}

function parse_args() {
  while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
    -i | --input_csv)
      INPUT_CSV="$2"
      shift
      shift
      ;;
    -o | --output_dir)
      OUTPUT_DIR="$2"
      shift
      shift
      ;;
    -p | --gene_panel)
      GENE_PANEL="$2"
      shift
      shift
      ;;
    -h | --help)
      usage
      ;;
    *)
      echo "Unknown option: $1"
      usage
      ;;
    esac
  done

  if [[ -z "${INPUT_CSV:-}" || -z "${OUTPUT_DIR:-}" || -z "${GENE_PANEL:-}" ]]; then
    usage
  fi

  if [[ $GENE_PANEL == "exome" ]]; then
    INTERVALS_BED="$EXOME_INTERVALS"
  elif [[ $GENE_PANEL == "hyperchol" ]]; then
    INTERVALS_BED="$HYPERCHOL_INTERVALS"
  elif [[ $GENE_PANEL == "ccp" ]]; then
    INTERVALS_BED="$CCP_INTERVALS"
  elif [[ $GENE_PANEL == "hcmp" ]]; then
    INTERVALS_BED="$HCMP_INTERVALS"
  elif [[ $GENE_PANEL == "lqt" ]]; then
    INTERVALS_BED="$LQT_INTERVALS"
  else
    echo "ERROR: Unknown gene panel: $GENE_PANEL"
    echo "Use one of the following gene panels: exome, hyperchol, ccp, hcmp, lqt."
    exit 1
  fi

  TIMESTAMP=$(date +'%Y%m%d')
  LOG_FILE="${TIMESTAMP}_${OUTPUT_DIR}_sarek_pipeline.log"

}

function run_sarek() {
  log "Input sample sheet: $INPUT_CSV"
  log "Pipeline output folder: $OUTPUT_DIR"
  log "Using intervals: $INTERVALS_BED"

  if [[ ! -f "$INPUT_CSV" ]]; then
    echo "ERROR: Input CSV file not found: $INPUT_CSV"
    exit 1
  fi

  if [[ ! -f "$INTERVALS_BED" ]]; then
    echo "ERROR: Intervals BED file not found: $INTERVALS_BED"
    exit 1
  fi

  local cmd="nextflow run nf-core/sarek -r 3.4.0 \
    -profile docker \
    -c $CONFIG_FILE \
    -params-file $PARAM_FILE \
    --intervals INTERVALS_BED \
    --input $INPUT_CSV \
    --outdir results-$OUTPUT_DIR \
    -resume"

  log_and_run "$cmd"
}

function run_vcf_preparation() {
  log "Running Nexflow VCF preparation script."
  input_table="results-${OUTPUT_DIR}/csv/variantcalled.csv"
  output_dir="results-${OUTPUT_DIR}/variant_calling/vcf_preparation/"
  local cmd="nextflow run ${SELF_PATH}/preparation/vcf_preparation.nf \
    --csv_input $input_table \
    --output_dir $output_dir \
    --reference $REFERENCE_GENOME \
    -c ${SELF_PATH}/preparation/nextflow.config \
    -resume \
    -dump-channels"

  log_and_run "$cmd"
}

function main() {
  source_env
  parse_args "$@"
  activate_conda
  # cleanup_docker
  # run_sarek
  run_vcf_preparation
}

main "$@"
