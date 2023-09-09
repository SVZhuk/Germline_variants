#!/bin/bash

FILE=$1

VepCache="/mnt/hdd/nextflow_dir/workspace/REF/VEP_108_cache/"
Plugins="/mnt/hdd/nextflow_dir/workspace/REF/VEP_108_cache/Plugins/"
CaddSnps="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/CADD_whole_genome_SNVs.tsv.gz"
CaddIndels="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/CADD_gnomad.genomes.r3.0.indel.tsv.gz"
LofteePath="/mnt/hdd/nextflow_dir/workspace/REF/VEP_108_cache/Plugins/loftee/"
human_ancestor_fa="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/human_ancestor.fa.gz"
conservation_file="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/phylocsf_gerp.sql"
ClinVar="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/clinvar_20230520.vcf.gz"
TopMed="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/TOPMED_GRCh38_20180418.vcf.gz"
dbNSFPfile="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/dbNSFP4.2a_grch38.gz"
replacement_logic="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/dbNSFP_replacement_logic.txt"
dbscSNVfile="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/dbscSNV1.1_GRCh38.txt.gz"
spliceai_snv="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/spliceai_scores.raw.snv.hg38.vcf.gz"
spliceai_indel="/mnt/hdd/nextflow_dir/workspace/REF/VEP_cache/Plugins/spliceai_scores.raw.indel.hg38.vcf.gz"

ID=$(basename ${FILE%.*})
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
    --format vcf \
    --plugin CADD,$CaddSnps,$CaddIndels \
    --plugin LoF,loftee_path:${Plugins},human_ancestor_fa:$human_ancestor_fa,conservation_file:$conservation_file \
    --custom $ClinVar,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
    --custom $TopMed,,vcf,exact,0,TOPMED \
    --plugin dbNSFP,$dbNSFPfile,$replacement_logic,gnomAD_genomes_AF,gnomAD_genomes_NFE_AF,gnomAD_genomes_POPMAX_AF,ClinPred_pred,clinvar_trait,clinvar_OMIM_id \
    --plugin dbscSNV,$dbscSNVfile \
    --plugin SpliceAI,snv=$spliceai_snv,indel=$spliceai_indel,cutoff=0.5 \
    --tab \
    --input_file $FILE \
    --output_file ${ID}.VEP.RefSeq.tsv \
    --force_overwrite
