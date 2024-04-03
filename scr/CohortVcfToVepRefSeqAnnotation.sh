#!/bin/bash

# shellcheck source=/mnt/hdd/nextflow_dir/workspace/anaconda3/etc/profile.d/conda.sh
source /mnt/hdd/nextflow_dir/workspace/anaconda3/etc/profile.d/conda.sh
conda activate gatk

InputVCF=$1
OutputDir=$2
FileName=$(basename ${InputVCF%.vcf.gz})

if [[ $# -eq 0 ]]; then
  echo "No arguments provided!
        Requires positional arguments:
        1: Path to cohort VCF after VQSR recalibration
        2: Path to output directory
        Exiting..."
fi

# Make separate folders for VCF files and TSV tables
mkdir -p "$OutputDir"/{TSVs,VCFs}

echo "Input Cohort VCF file is located at: $InputVCF."
echo "Output directory for per patient VCF files and TSV tables is: $OutputDir"

Reference="/home/students/nextflow_dir/workspace/REF/iGenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"

bcftools view -i 'FILTER="PASS"' ${InputVCF} >${OutputDir}/${FileName}.bcftools.PASS.vcf
InputVCF=${OutputDir}/${FileName}.bcftools.PASS.vcf

for sample in $(bcftools query -l ${InputVCF}); do

  echo "Writing $sample patients's data to a separate VCF file in $OutputDir."

  # -c1 option drops empty genotypes
  bcftools view -c1 -s "$sample" -Ov -o "$OutputDir"/VCFs/"$sample".vcf "$InputVCF" \
    2>&1 | tee -a "$OutputDir"/SplittingVcfLog.txt

  echo "Splitting multiallelic sites."

  bcftools norm -m-both -o "$OutputDir"/VCFs/"$sample".step1.vcf "$OutputDir"/VCFs/"$sample".vcf \
    2>&1 | tee -a "$OutputDir"/SplittingVcfLog.txt

  echo "Left-normalization of variants."

  bcftools norm -f "$Reference" -o "$OutputDir"/VCFs/"$sample".lnorm.vcf "$OutputDir"/VCFs/"$sample".step1.vcf \
    2>&1 | tee -a "$OutputDir"/SplittingVcfLog.txt

  rm "$OutputDir"/VCFs/"$sample".vcf
  rm "$OutputDir"/VCFs/"$sample".step1.vcf

  echo "Indexing $sample VCF file for gatk tool."

  gatk IndexFeatureFile \
    -I "$OutputDir"/VCFs/"$sample".lnorm.vcf

  echo "Writing genotype information from $sample VCF file to Table."

  gatk VariantsToTable \
    -V "$OutputDir"/VCFs/"$sample".lnorm.vcf \
    -F CHROM -F POS -F TYPE -F REF -F ALT -F FILTER -GF AD -GF DP -GF GQ \
    -O "$OutputDir"/TSVs/"$sample".table \
    2>&1 | tee -a "$OutputDir"/SplittingVcfLog.txt

done

echo "Deactivating gatk environment."
conda deactivate

InputDir=$OutputDir/VCFs/

# Make directory for VEP outputs
mkdir -p "$OutputDir"/VEP

VepCache="/mnt/hdd/nextflow_dir/workspace/REF/VEP_108_cache/"
Plugins="/mnt/hdd/nextflow_dir/workspace/REF/VEP_108_cache/Plugins/"
CaddSnps="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/CADD_whole_genome_SNVs.tsv.gz"
CaddIndels="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/CADD_gnomad.genomes.r3.0.indel.tsv.gz"
LofteePath="/mnt/hdd/nextflow_dir/workspace/REF/VEP_108_cache/Plugins/loftee/"
human_ancestor_fa="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/human_ancestor.fa.gz"
conservation_file="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/phylocsf_gerp.sql"
ClinVar="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/clinvar_20221211.vcf.gz"
TopMed="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/TOPMED_GRCh38_20180418.vcf.gz"
dbNSFPfile="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/dbNSFP4.2a_grch38.gz"
replacement_logic="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/dbNSFP_replacement_logic.txt"
dbscSNVfile="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/dbscSNV1.1_GRCh38.txt.gz"
spliceai_snv="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/spliceai_scores.raw.snv.hg38.vcf.gz"
spliceai_indel="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/spliceai_scores.raw.indel.hg38.vcf.gz"

for file in "$InputDir"/*.lnorm.vcf; do

  ID=$(basename "${file%.lnorm.vcf}")
  perl /home/students/nextflow_dir/workspace/ensembl-vep/vep \
    --offline \
    --cache \
    --refseq \
    --dir_cache ${VepCache} \
    --dir_plugins ${Plugins} \
    --no_stats \
    --hgvs \
    --hgvsg \
    --fork 20 \
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
    --check_existing \
    --check_frequency \
    --freq_pop 1KG_EUR \
    --freq_freq 0.01 \
    --freq_gt_lt gt \
    --freq_filter exclude \
    --plugin CADD,$CaddSnps,$CaddIndels \
    --plugin LoF,loftee_path:$LofteePath,human_ancestor_fa:$human_ancestor_fa,conservation_file:$conservation_file \
    --custom $ClinVar,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
    --custom $TopMed,,vcf,exact,0,TOPMED \
    --plugin dbNSFP,$dbNSFPfile,$replacement_logic,gnomAD_genomes_AF,gnomAD_genomes_NFE_AF,gnomAD_genomes_POPMAX_AF,ClinPred_pred,clinvar_trait,clinvar_OMIM_id \
    --plugin dbscSNV,$dbscSNVfile \
    --plugin SpliceAI,snv=$spliceai_snv,indel=$spliceai_indel,cutoff=0.5 \
    --format vcf \
    --tab \
    --input_file "$file" \
    --output_file "$OutputDir"/VEP/"$ID".VEP.RefSeq.tsv \
    --warning_file "$OutputDir"/VEP/"$ID".VEP.RefSeq.error.log.txt \
    --force_overwrite

  echo "Done with annotation of $file file!"

done

exit
