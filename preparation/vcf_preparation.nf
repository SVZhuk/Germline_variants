#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.csv_input = params.csv_input
params.output_dir = params.output_dir
params.reference = params.reference
params.metadata = params.metadata
params.expression_table = params.expression_table

// Define the main workflow
workflow {
    // Create input channel from CSV file
    input_ch = Channel
        .fromPath(params.csv_input)
        .splitCsv(header: true)
        .map { row -> tuple(row.patient, row.sample, row.variantcaller, file(row.vcf)) }

    adjust_script = file("report_gen/adjust_genotype_table.py")
    report_script = file("report_gen/prepare_report_table.py")

    // Split multiallelic sites
    SPLIT_MULTIALLELIC(input_ch)

    // Left-normalize variants
    LEFT_NORMALIZE(SPLIT_MULTIALLELIC.out.split_vcf, params.reference)

    // Index VCF files
    INDEX_VCF(LEFT_NORMALIZE.out.normalized_vcf)

    // Convert variants to table
    VARIANTS_TO_TABLE(INDEX_VCF.out.indexed_vcf)

    // Annotate with VEP
    VEP_ANNOTATION(LEFT_NORMALIZE.out.normalized_vcf)

    // Adjust genotype table
    ADJUST_GENOTYPE_TABLE(VARIANTS_TO_TABLE.out.table_out, adjust_script)

    // Prepare report table
    vep_and_genotype = VEP_ANNOTATION.out.vep_out
        .join(ADJUST_GENOTYPE_TABLE.out.adjusted_table, by: [0,1,2])
        .map { patient, sample, variantcaller, vep_table, adjusted_table ->
            println "Debug: $patient, $sample, $variantcaller, $vep_table, $adjusted_table"
        tuple(patient, sample, variantcaller, vep_table, adjusted_table)
    }

    PREPARE_REPORT_TABLE(
        vep_and_genotype,
        report_script,
        file(params.metadata),
        file(params.expression_table)
    )
}

// Define processes
process SPLIT_MULTIALLELIC {
    container 'community.wave.seqera.io/library/bcftools:812c399a55984a80'
    
    input:
    tuple val(patient), val(sample), val(variantcaller), path(vcf)

    output:
    tuple val(patient), val(sample), val(variantcaller), path("${sample}.step1.vcf"), emit: split_vcf

    script:
    """
    bcftools norm -m-both -o ${sample}.step1.vcf ${vcf}
    """
}

process LEFT_NORMALIZE {
    container 'community.wave.seqera.io/library/bcftools:812c399a55984a80'

    publishDir "${params.output_dir}/VCFs", mode: 'copy'

    input:
    tuple val(patient), val(sample), val(variantcaller), path(vcf)
    path reference

    output:
    tuple val(patient), val(sample), val(variantcaller), path("${sample}.lnorm.vcf"), emit: normalized_vcf

    script:
    """
    bcftools norm -f ${reference} -o ${sample}.lnorm.vcf ${vcf}
    """
}

process INDEX_VCF {
    container 'community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867'

    publishDir "${params.output_dir}/VCFs", mode: 'copy'

    input:
    tuple val(patient), val(sample), val(variantcaller), path(vcf)

    output:
    tuple val(patient), val(sample), val(variantcaller), path("${vcf}"), path("${vcf}.idx"), emit: indexed_vcf

    script:
    """
    gatk IndexFeatureFile -I ${vcf}
    """
}

process VARIANTS_TO_TABLE {
    container 'community.wave.seqera.io/library/gatk4:4.5.0.0--730ee8817e436867'

    publishDir "${params.output_dir}/TSVs", mode: 'copy'

    input:
    tuple val(patient), val(sample), val(variantcaller), path(vcf), path(idx)

    output:
    tuple val(patient), val(sample), val(variantcaller), path("${sample}.table"), emit: table_out

    script:
    """
    gatk VariantsToTable \
        -V ${vcf} \
        -F CHROM -F POS -F TYPE -F REF -F ALT -F FILTER -GF AD -GF DP -GF GQ \
        -O ${sample}.table \
        --show-filtered
    """
}

process ADJUST_GENOTYPE_TABLE {
    publishDir "${params.output_dir}/TSVs", mode: 'copy'

    input:
    tuple val(patient), val(sample), val(variantcaller), path(table)
    path(script)

    output:
    tuple val(patient), val(sample), val(variantcaller), 
    path("${sample}.mod.tsv"), emit: adjusted_table

    script:
    """
    cp ${table} ./input_table.tsv
    python ${script} --file_path input_table.tsv
    mv input_table.mod.tsv ${sample}.mod.tsv
    """
}

process VEP_ANNOTATION {
    container 'ensemblorg/ensembl-vep:release_108.0'

    errorStrategy 'terminate'

    publishDir "${params.output_dir}/VEP", mode: 'copy'

    input:
    tuple val(patient), val(sample), val(variantcaller), path(vcf)

    output:
    tuple val(patient), val(sample), val(variantcaller), 
    path("${sample}.VEP.RefSeq.tsv"), emit: vep_out

    script:
    """
    mkdir -p \${PWD}/vep_output
    chmod 777 \${PWD}/vep_output

    vep \
    --offline \
    --cache \
    --refseq \
    --dir_cache ${params.vep_cache} \
    --dir_plugins ${params.vep_plugins} \
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
    --plugin CADD,${params.cadd_snps},${params.cadd_indels} \
    --plugin LoF,loftee_path:${params.vep_plugins},human_ancestor_fa:${params.human_ancestor_fa},conservation_file:${params.conservation_file} \
    --custom ${params.clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
    --custom ${params.topmed},,vcf,exact,0,TOPMED \
    --plugin dbNSFP,${params.dbnsfp_file},${params.replacement_logic},gnomAD_genomes_AF,gnomAD_genomes_NFE_AF,gnomAD_genomes_POPMAX_AF,ClinPred_pred,clinvar_trait,clinvar_OMIM_id \
    --plugin dbscSNV,${params.dbscsnv_file} \
    --plugin SpliceAI,snv=${params.spliceai_snv},indel=${params.spliceai_indel},cutoff=0.5 \
    --format vcf \
    --tab \
    --input_file ${vcf} \
    --output_file \${PWD}/vep_output/${sample}.VEP.RefSeq.tsv \
    --warning_file \${PWD}/vep_output/${sample}.VEP.RefSeq.tsv_warnings.txt \
    --force_overwrite

    mv \${PWD}/vep_output/${sample}.VEP.RefSeq.tsv .
    """
}

process PREPARE_REPORT_TABLE {
    publishDir "${params.output_dir}/annotation", mode: 'copy'

    input:
    tuple val(patient), val(sample), val(variantcaller), path(vep_table), path(adjusted_table)
    path(script)
    path(metadata)
    path(expression_table)

    output:
    path "${sample}.report.tsv", emit: report

    script:
    """
    python ${script} --vep_table ${vep_table} \
        --genotype_table ${adjusted_table} \
        --metadata ${metadata} \
        --expression_table ${expression_table} \
        --output ${sample}.report.tsv
    """
}
