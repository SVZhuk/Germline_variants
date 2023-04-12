#!/bin/bash

# Script for variant calling from SureSelect capture experiments
# Uses nf-core conda env
# Uses local sarek pipeline

PipelineDir="/home/students/nextflow_dir/pipelines/nf-core-sarek-3.1/workflow/"

InputCsv=$1
Probes="/home/students/nextflow_dir/workspace/REF/agilent/SureSelectHumanAllExonV6r2/SureSelectAllExonV6r2_100bp_hg38.sorted.bed"
Name=$2

echo "Input file used is: $InputCsv."
echo "Probes used for capture are: $Probes."
echo "Name of output folder for results is: $Name."

if [[ -z "$InputCsv" ]]; then
  echo "Input first positional argument: TSV with samples! Exiting..."
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

iGenomes_base="/mnt/hdd/nextflow_dir/workspace/REF/iGenomes/"

nextflow run ${PipelineDir}  \
    -profile docker \
    --input ${InputCsv} \
    --joint_germline false \
    --outdir results-${Name} \
    --genome GATK.GRCh38 \
    --igenomes_base ${iGenomes_base} \
    --wes true \
    --intervals ${Probes} \
    --tools haplotypecaller,strelka,freebayes,deepvariant \
    --save_output_as_bam \
    --max_cpus 80 \
    --max_memory 180.GB \
    -resume

