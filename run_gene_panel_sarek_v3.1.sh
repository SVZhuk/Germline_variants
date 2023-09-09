#!/bin/bash

set -ex -o pipefail

# Script for variant calling from SureSelect capture experiments
# Uses nf-core conda env
# Uses local sarek pipeline

PipelineDir="/home/students/nextflow_dir/pipelines/nf-core-sarek-3.1/workflow/"

InputCsv=$1
Probes=$2
Name=$3

echo "Input file used is: $InputCsv."
echo "Probes used for capture are: $Probes."
echo "Name of output folder for results is: $Name."

if [[ $# -eq 0 ]]; then
  echo "No arguments provided! Exiting..."
  exit 1
fi

if [[ -z "$InputCsv" ]]; then
  echo "Input first positional argument: TSV with samples! Exiting..."
  exit 1
fi

if [[ -z "$Probes" ]]; then
  echo "Input second positional argument: bed file with sorted probe intervals! Exiting..."
  exit 1
fi

if [[ -z "$Name" ]]; then
  echo "Input third positional argument: name for results output folder! Exiting..."
  exit 1
fi

echo "Activating nf-core conda environment."
source /mnt/hdd/nextflow_dir/workspace/anaconda3/etc/profile.d/conda.sh
conda activate nf-core

echo "Currently running shell from $(pwd)"

docker system prune --all --force
docker image prune --all --force

iGenomes_base="/mnt/hdd/nextflow_dir/workspace/REF/iGenomes/"

nextflow run ${PipelineDir} \
  -profile docker \
  --input ${InputCsv} \
  --outdir results-${Name} \
  --genome GATK.GRCh38 \
  --igenomes_base ${iGenomes_base} \
  --wes true \
  --intervals ${Probes} \
  --tools haplotypecaller,deepvariant \
  --save_output_as_bam \
  --max_cpus 80 \
  --max_memory 180.GB \
  -resume
