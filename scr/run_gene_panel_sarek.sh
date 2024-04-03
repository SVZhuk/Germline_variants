#!/bin/bash

# Script for variant calling from SureSelect capture experiments
# Uses nf-core conda env
# Uses local sarek pipeline

Ref="/mnt/hdd/nextflow_dir/workspace/REF/iGenomes/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
RefInd="/mnt/hdd/nextflow_dir/workspace/REF/iGenomes/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta.fai"
RefDict="/mnt/hdd/nextflow_dir/workspace/REF/iGenomes/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict"
BWAind="/mnt/hdd/nextflow_dir/workspace/REF/iGenomes/GATK/GRCh38/Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64.{alt,amb,ann,bwt,pac,sa}"
dbsnp="/mnt/hdd/nextflow_dir/workspace/REF/iGenomes/GATK/GRCh38/Annotation/GATKBundle/dbsnp_138.hg38.vcf.gz"
dbsnpIndex="/mnt/hdd/nextflow_dir/workspace/REF/iGenomes/GATK/GRCh38/Annotation/GATKBundle/dbsnp_138.hg38.vcf.gz.tbi"
indels="/mnt/hdd/nextflow_dir/workspace/REF/iGenomes/GATK/GRCh38/Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
indelsIndex="/mnt/hdd/nextflow_dir/workspace/REF/iGenomes/GATK/GRCh38/Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
intervals="/mnt/hdd/nextflow_dir/workspace/REF/iGenomes/GATK/GRCh38/Annotation/intervals/wgs_calling_regions.hg38.bed"

docker system prune --all --force
docker image prune --all --force

InputTsv=$1
Probes=$2
Name=$3

echo "Input file used is: $InputTsv."
echo "Probes used for capture are: $Probes."
echo "Name of output folder for results is: $Name."

if [[ -z "$InputTsv" ]]; then
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

nextflow run nf-core/sarek -r 2.7.1 -profile docker \
  --input $InputTsv \
  --step mapping \
  --generate_gvcf \
  --tools HaplotypeCaller \
  --outdir ./results-$Name \
  --genome GRCh38 \
  --target_bed $Probes \
  --intervals $intervals \
  --dict $RefDict \
  --fasta $Ref \
  --fasta_fai $RefInd \
  --bwa $BWAind \
  --dbsnp $dbsnp \
  --dbsnp_index $dbsnpIndex \
  --known_indels $indels \
  --known_indels_index $indelsIndex \
  --max_cpus 80 \
  --single_cpu_mem '8.GB' \
  --max_memory '150.GB'

exit
