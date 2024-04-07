#!/usr/bin/env bash

###### Request resources ######

#SBATCH --job-name=expansionHunter
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBARCH --cpus-per-task=40
#SBATCH --partition=mid
#SBATCH --time=24:00:00
#SBATCH --output=atxn-repeats.out
#SBATCH --mail-type=NONE
#SBATCH --mail-user=szhuk@ku.edu.tr

###############################

set -euo pipefail

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -i|--input_dir)
      INPUT_DIR="$2"
      shift
      shift
      ;;
    -h|--help)
      echo "Usage: ExpansionHunter.sh -i|--input_dir <input_dir with cram/bam>"
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      exit 1
      ;;
  esac
done

module load samtools/1.14

VARIANT_CATALOG="/scratch/users/szhuk/hpc_run/workspace/ExpansionHunter/variant_catalog/grch38/variant_catalog.json"
REFERENCE_FASTA="/kuttam_fg/refdata/zhuk/iGenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"

echo -e "VariantId\tGenotype\tAlleleDepth" > metrics_file.tsv

for file in $(find ${INPUT_DIR} -type f -name "*.cram"); do
    if [[ ! -f "${file}" ]] || [[ ! -r "${file}" ]]; then
        echo "File not found or not readable: ${file}"
        exit 1
    fi

    ID=$(basename "${file}" | cut -d'_' -f1)

    ExpansionHunter --reads $file \
                    --reference $REFERENCE_FASTA \
                    --variant-catalog $VARIANT_CATALOG \
                    --output-prefix $ID \
                    --threads 30 || { echo "ExpansionHunter failed for ID: $ID"; continue; }

    BAM_FILE="${ID}_realigned.bam"

    samtools sort $BAM_FILE -o ${ID}_sorted.bam
    samtools index ${ID}_sorted.bam

    REViewer-v0.2.7-linux_x86_64 --reads ${ID}_sorted.bam \
        --vcf ${ID}.vcf \
        --reference $REFERENCE_FASTA \
        --catalog $VARIANT_CATALOG \
        --locus ATXN1,ATXN2,ATXN3 --output-prefix ATXN_${ID}

    sed -n '2,4p' ATXN_${ID}.metrics.tsv >> metrics_file.tsv

    rm ${ID}_sorted.bam*

done
