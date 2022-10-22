#!/bin/bash

# Script takes 4 parameters as input arguments:
# 1 - path to sorted bed file with capture probe information
# 2 - path to folder with GVCF files (uses find to make a list of files)
# 3 - creates a folder with Cohort VCF
# 4 - specifies ID for output cohort VCF

Probes=$1
GvcfDir=$2
CohortVcfDir=$3
CohortVcfId=$4

echo "Probe Set is located at: $Probes."
echo "GVCFs directory is located at: $GvcfDir."
echo "Cohort VCF file output directory is located at: $CohortVcfDir."
echo "Final file ID is: $CohortVcfId."

if [[ -z "$Probes" ]]; then
  echo "Provide bed file with probe information as first positional argument! Exiting..."
  exit 1
fi

if [[ -z "$GvcfDir" ]]; then
  echo "Provide path to directory with GVCFs for aggregation as second positional argument! Exiting..."
  exit 1
fi

if [[ -z "$CohortVcfDir" ]]; then
  echo "Needed output folder as third positional argument!"
  exit 1
fi

if [[ -z "$CohortVcfId" ]]; then
  echo "Provide ID for final file as forth positional argument! Exiting..."
  exit 1
fi

mkdir -p "$CohortVcfDir"

# shellcheck source=/mnt/hdd/nextflow_dir/workspace/anaconda3/etc/profile.d/conda.sh
source /mnt/hdd/nextflow_dir/workspace/anaconda3/etc/profile.d/conda.sh
conda activate gatk

# find and list all gVCF files of interst
find "$GvcfDir" -iname "*.g.vcf.gz" > "$CohortVcfDir"/gvcf.list

# Reference files
Reference="/mnt/hdd/nextflow_dir/workspace/REF/iGenomes/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta"
dbSNP="/mnt/hdd/nextflow_dir/workspace/REF/iGenomes/GATK/GRCh38/Annotation/GATKBundle/dbsnp_138.hg38.vcf.gz"
RefGenomeDict="/mnt/hdd/nextflow_dir/workspace/REF/iGenomes/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict"

# Convert gVCF files into cohort VCF file:
gatk --java-options "-Xms24G -Xmx64G -XX:ParallelGCThreads=2" GenomicsDBImport \
  -R "$Reference" \
  --intervals "$Probes" \
  --merge-input-intervals \
  --genomicsdb-workspace-path "$CohortVcfDir/$CohortVcfId"_database \
  --batch-size 100 \
  --max-num-intervals-to-import-in-parallel 10 \
  -V "$CohortVcfDir"/gvcf.list \
  2>&1 | tee -a "$CohortVcfId"_Annotation_Log.txt

gatk --java-options "-Xmx64G -XX:ParallelGCThreads=2" GenotypeGVCFs \
  -R "$Reference" \
  -V gendb://"$CohortVcfDir/$CohortVcfId"_database \
  -O "$CohortVcfDir/$CohortVcfId".joint.vcf \
  2>&1 | tee -a "$CohortVcfId"_Annotation_Log.txt

# Hard filtering of variants
# Select SNPs and INDELs separately
mkdir -p "$CohortVcfDir"/HardFiltering
HfDir="$CohortVcfDir"/HardFiltering

gatk --java-options "-Xmx60G" SelectVariants \
  -R "$Reference" \
  --intervals "$Probes" \
  -V "$CohortVcfDir/$CohortVcfId".joint.vcf \
  -O "$HfDir/$CohortVcfId".snp.unfiltered.vcf \
  -select-type SNP \
  2>&1 | tee -a "$CohortVcfId"_Annotation_Log.txt

gatk --java-options "-Xmx60G" SelectVariants \
  -R "$Reference" \
  --intervals "$Probes" \
  -V "$CohortVcfDir/$CohortVcfId".joint.vcf \
  -O "$HfDir/$CohortVcfId".indel.unfiltered.vcf \
  -select-type INDEL \
  2>&1 | tee -a "$CohortVcfId"_Annotation_Log.txt

# Filter variants
gatk --java-options "-Xmx64G" VariantFiltration \
  -R "$Reference" \
  -V "$HfDir/$CohortVcfId".snp.unfiltered.vcf \
  --filter-expression "QD < 2.0" --filter-name "QD_lt_2" \
  --filter-expression "SOR > 3.0" --filter-name "SOR_gt_3" \
  --filter-expression "QUAL < 30.0" --filter-name "QUAL_lt_30" \
  --filter-expression "FS > 60.0" --filter-name "FS_gt_60" \
  --filter-expression "MQ < 40.0" --filter-name "MQ_lt_40" \
  --filter-expression "MQRankSum < -12.5" --filter-name "MQRS_lt_n12.5" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "RPRS_lt_n8" \
  -O "$HfDir/$CohortVcfId".snp.filtered.vcf \
  2>&1 | tee -a "$CohortVcfId"_Annotation_Log.txt

gatk --java-options "-Xmx64G" VariantFiltration \
  -R "$Reference" \
  -V "$HfDir/$CohortVcfId".indel.unfiltered.vcf \
  --filter-expression "QD < 2.0" --filter-name "QD_lt_2" \
  --filter-expression "QUAL < 30.0" --filter-name "QUAL_lt_30" \
  --filter-expression "FS > 200.0" --filter-name "FS_gt_200" \
  --filter-expression "ReadPosRankSum < -20.0" --filter-name "RPRSum_lt_n20" \
  -O "$HfDir/$CohortVcfId".indel.filtered.vcf \
  2>&1 | tee -a "$CohortVcfId"_Annotation_Log.txt

gatk --java-options "-Xmx64G" MergeVcfs \
  -I "$HfDir/$CohortVcfId".indel.filtered.vcf \
  -I "$HfDir/$CohortVcfId".snp.filtered.vcf \
  -O "$HfDir/$CohortVcfId".all.filtered.vcf \
  2>&1 | tee -a "$CohortVcfId"_Annotation_Log.txt

echo "Collecting statistics on Hard Filtered Files..."

gatk CollectVariantCallingMetrics \
  -I "$HfDir/$CohortVcfId".all.filtered.vcf \
  --DBSNP "$dbSNP" \
  -SD "$RefGenomeDict" \
  -O "$HfDir/$CohortVcfId".all.filtered.metrics \
  2>&1 | tee -a "$CohortVcfId"_Annotation_Log.txt

rm "$HfDir"/*.snp.filtered.vcf
rm "$HfDir"/*.snp.filtered.vcf.idx
rm "$HfDir"/*.indel.filtered.vcf
rm "$HfDir"/*.indel.filtered.vcf.idx
rm "$HfDir"/*.unfiltered.vcf
rm "$HfDir"/*.unfiltered.vcf.idx

InputVCF="$HfDir/$CohortVcfId".all.filtered.vcf
OutputDir="$HfDir"

# Make separate folders for VCF files and TSV tables
mkdir -p "$OutputDir"/{TSVs,VCFs}

for sample in $(bcftools query -l "$InputVCF"); do

  echo "Writing $sample patients's data to a separate VCF file in $OutputDir."

  # -c1 option drops empty genotypes
  bcftools view -c1 -s "$sample" -Ov -o "$OutputDir/VCFs/$sample".vcf "$InputVCF" \
    2>&1 | tee -a "$CohortVcfId"_Annotation_Log.txt

  echo "Splitting multiallelic sites."

  bcftools norm -m-both -o "$OutputDir/VCFs/$sample".step1.vcf "$OutputDir/VCFs/$sample".vcf \
    2>&1 | tee -a "$CohortVcfId"_Annotation_Log.txt

  echo "Left-normalization of variants."

  bcftools norm -f "$Reference" -o "$OutputDir/VCFs/$sample".lnorm.vcf "$OutputDir/VCFs/$sample".step1.vcf \
    2>&1 | tee -a "$CohortVcfId"_Annotation_Log.txt

  rm "$OutputDir/VCFs/$sample".vcf
  rm "$OutputDir/VCFs/$sample".step1.vcf

  echo "Indexing $sample VCF file for gatk tool."

  gatk IndexFeatureFile \
    -I "$OutputDir/VCFs/$sample".lnorm.vcf

  echo "Writing genotype information from $sample VCF file to Table."

  gatk VariantsToTable \
    -V "$OutputDir/VCFs/$sample".lnorm.vcf \
    -F CHROM -F POS -F TYPE -F REF -F ALT -F FILTER -GF AD -GF DP -GF GQ \
    -O "$OutputDir/TSVs/$sample".table \
    --show-filtered \
    2>&1 | tee -a "$CohortVcfId"_Annotation_Log.txt

done

echo "Deactivating gatk environment."
conda deactivate

InputDir=$OutputDir/VCFs/

conda activate vep

# Make directory for VEP outputs
mkdir -p "$OutputDir"/VEP

VepCache="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/"
Plugins="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/"
Ref="/mnt/hdd/nextflow_dir/workspace/REF/Ref_GRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
Bam="/mnt/hdd/nextflow_dir/workspace/REF/Ref_GRCh38/GCF_000001405.39_GRCh38.p13_knownrefseq_alns.bam"
CaddSnps="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/CADD_whole_genome_SNVs.tsv.gz"
CaddIndels="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/CADD_gnomad.genomes.r3.0.indel.tsv.gz"
LofteePath="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/loftee/"
human_ancestor_fa="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/human_ancestor.fa.gz"
conservation_file="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/phylocsf_gerp.sql"
ClinVar="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/clinvar_20220709.vcf.gz"
TopMed="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/TOPMED_GRCh38_20180418.vcf.gz"
dbNSFPfile="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/dbNSFP4.2a_grch38.gz"
replacement_logic="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/dbNSFP_replacement_logic.txt"
dbscSNVfile="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/dbscSNV1.1_GRCh38.txt.gz"
spliceai_snv="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/spliceai_scores.raw.snv.hg38.vcf.gz"
spliceai_indel="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/spliceai_scores.raw.indel.hg38.vcf.gz"

for file in "$InputDir"/*.lnorm.vcf; do

  ID=$(basename "${file%.lnorm.vcf}")
  vep --offline \
    --cache \
    --refseq \
    --dir_cache "$VepCache" \
    --dir_plugins "$Plugins" \
    --no_stats \
    --fasta "$Ref" \
    --bam "$Bam" \
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

echo "Deactivating vep environment."
conda deactivate

exit
