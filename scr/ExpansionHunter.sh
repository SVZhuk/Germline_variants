#!/usr/bin/env bash

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -i|--input_cram)
      INPUT_CRAM="$2"
      shift
      shift
      ;;
    -o|--output_index)
      OUTPUT_IDX="$2"
      shift
      shift
      ;;
    -h|--help)
      echo "Usage: ExpansionHunter.sh -i|--input_cram <input_cram/bam> -o|--output_index <output_file_index>"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

ExpansionHunter="/scratch/users/szhuk/hpc_run/workspace/ExpansionHunter/build/install/bin/ExpansionHunter"
VARIANT_CATALOG="/scratch/users/szhuk/hpc_run/workspace/ExpansionHunter/variant_catalog/grch38/variant_catalog.json"
REFERENCE_FASTA="/kuttam_fg/refdata/zhuk/iGenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"

ExpansionHunter --reads $INPUT_CRAM \
                --reference $REFERENCE_FASTA \
                --variant-catalog $VARIANT_CATALOG \
                --output-prefix $OUTPUT_IDX

module load samtools/1.14