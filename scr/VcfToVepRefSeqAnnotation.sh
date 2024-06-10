#!/bin/bash

set -ex -o pipefail

# shellcheck source=/mnt/hdd/nextflow_dir/workspace/anaconda3/etc/profile.d/conda.sh
source /mnt/hdd/nextflow_dir/workspace/anaconda3/etc/profile.d/conda.sh
conda activate gatk

InputDir=$1
OutputDir=$2

# Make separate folders for VCF files and TSV tables
mkdir -p "$OutputDir"/{TSVs,VCFs}

echo "Input VCF files are located at: ${InputDir}"
echo "Output directory for per patient VCF files and TSV tables is: ${OutputDir}"

for file in $(find ${InputDir} \( -name "*.haplotypecaller.filtered.vcf.gz" -o -name "*.deepvariant.vcf.gz" \)); do

  echo "Splitting multiallelic sites."
  ID=$(basename ${file%.vcf.gz})
  Reference="/home/students/nextflow_dir/workspace/REF/iGenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"

  bcftools norm -m-both -o ${OutputDir}/VCFs/${ID}.step1.vcf ${file}

  echo "Left-normalization of variants."

  bcftools norm -f ${Reference} -o ${OutputDir}/VCFs/${ID}.lnorm.vcf ${OutputDir}/VCFs/${ID}.step1.vcf

  rm ${OutputDir}/VCFs/${ID}.step1.vcf

  echo "Indexing $sample VCF file for gatk tool."

  gatk IndexFeatureFile \
    -I ${OutputDir}/VCFs/${ID}.lnorm.vcf

  echo "Writing genotype information from ${file} VCF file to Table."

  gatk VariantsToTable \
    -V ${OutputDir}/VCFs/${ID}.lnorm.vcf \
    -F CHROM -F POS -F TYPE -F REF -F ALT -F FILTER -GF AD -GF DP -GF GQ \
    -O ${OutputDir}/TSVs/${ID}.table \
    --show-filtered

done

echo "Deactivating gatk environment."
conda deactivate

InputDir=${OutputDir}/VCFs/

# Make directory for VEP outputs
mkdir -p ${OutputDir}/VEP

VepCache="/mnt/hdd/nextflow_dir/workspace/REF/VEP_108_cache/"
Plugins="/mnt/hdd/nextflow_dir/workspace/REF/VEP_108_cache/Plugins/"
CaddSnps="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/CADD_whole_genome_SNVs.tsv.gz"
CaddIndels="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/CADD_gnomad.genomes.r3.0.indel.tsv.gz"
LofteePath="/mnt/hdd/nextflow_dir/workspace/REF/VEP_108_cache/Plugins/loftee/"
human_ancestor_fa="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/human_ancestor.fa.gz"
conservation_file="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/phylocsf_gerp.sql"
ClinVar="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/clinvar_20240528.vcf.gz"
TopMed="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/TOPMED_GRCh38_20180418.vcf.gz"
dbNSFPfile="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/dbNSFP4.2a_grch38.gz"
replacement_logic="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/dbNSFP_replacement_logic.txt"
dbscSNVfile="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/dbscSNV1.1_GRCh38.txt.gz"
spliceai_snv="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/spliceai_scores.raw.snv.hg38.vcf.gz"
spliceai_indel="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/spliceai_scores.raw.indel.hg38.vcf.gz"

for file in ${InputDir}/*.lnorm.vcf; do

  ID=$(basename ${file%.lnorm.vcf})
  perl /home/students/nextflow_dir/workspace/ensembl-vep/vep \
    --offline \
    --cache \
    --refseq \
    --dir_cache ${VepCache} \
    --dir_plugins ${Plugins} \
    --no_stats \
    --hgvs \
    --hgvsg \
    --fork 8 \
    --buffer_size 80000 \
    --assembly GRCh38 \
    --symbol \
    --show_ref_allele \
    --af_gnomad \
    --af \
    --af_1kg \
    --max_af \
    --numbers \
    --terms SO \
    --variant_class \
    --canonical \
    --mane \
    --exclude_predicted \
    --pick \
    --individual all \
    --sift b \
    --polyphen b \
    --pubmed \
    --plugin CADD,$CaddSnps,$CaddIndels \
    --plugin LoF,loftee_path:${Plugins},human_ancestor_fa:$human_ancestor_fa,conservation_file:$conservation_file \
    --custom $ClinVar,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
    --custom $TopMed,,vcf,exact,0,TOPMED \
    --plugin dbNSFP,$dbNSFPfile,$replacement_logic,gnomAD_genomes_AF,gnomAD_genomes_NFE_AF,gnomAD_genomes_POPMAX_AF,ClinPred_pred,clinvar_trait,clinvar_OMIM_id \
    --plugin dbscSNV,$dbscSNVfile \
    --plugin SpliceAI,snv=$spliceai_snv,indel=$spliceai_indel,cutoff=0.5 \
    --format vcf \
    --tab \
    --input_file ${file} \
    --output_file ${OutputDir}/VEP/${ID}.VEP.RefSeq.tsv \
    --force_overwrite

  echo "Done with annotation of $file file!"

done

exit
