# Load required libraries
library(readr)
library(dplyr)
library(stringr)
library(janitor)
library(openxlsx)

# set working directory to folder with VEP and GATK reports
# final reports are written here:
setwd("~/Desktop/AAK_lab/Variant_Annotation/data/2023.04.06_CCP_panel")
# Main directory for particular panel
# It needs to contain TSVs/ and VEP/ folders with respective outputs for proper parsing
MainDir <- "./"
# this will turn off scientific notation , need this for proper conversion of text to numbers 
options(scipen = 999) 
# Specify name of excel workbook:
file_ID <- "2023.04.06_CCP-panel"

################################################################################
# Prep files with genotype information for merging with VEP output
# Process all files in a loop.
# Point program to the folder where all genotype programs are stored:
deepvariant <- '.deepvariant'
haplotypecaller <- '.haplotypecaller.filtered'

caller <- deepvariant

files <- list.files(path = paste0(MainDir, "TSVs/"), full.names = TRUE, pattern = paste0(caller,'.table'))
for (file in files){
  file_path <- file
  file_name <- basename(file_path) %>% str_replace_all(paste0(caller, '.table'), paste0(caller, '.mod.tsv'))
  print(paste('Processing', file))
  df <- read.table(file = file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  df[, 'uploaded_variation'] <- NA
  df_indel <- df %>% filter(TYPE == 'INDEL')
  df$uploaded_variation <- with(df, paste0(CHROM, '_', POS, '_', REF, '/', ALT))
  df <- within(df, uploaded_variation[TYPE == 'INDEL'] <- ifelse(str_length(df_indel$REF) < str_length(df_indel$ALT), 
                                                                 paste0(df_indel$CHROM, '_', df_indel$POS, '_', '-', '/', sub('^.','', df_indel$ALT)), 
                                                                 paste0(df_indel$CHROM, '_', df_indel$POS + 1, '_', sub('^.', '', df_indel$REF), '/', '-')))
  genot_df_trimmed <- select(df, TYPE, 6:10)
  write.table(genot_df_trimmed, file = paste0(MainDir, "TSVs/", file_name), sep = '\t', row.names = FALSE)
  print(paste('Done processing', file_name))
}

################################################################################ 
# Merge vep reports with genotype info tables 
dir.create(file.path(MainDir, "Final/"))  
gen_files <- list.files(path = paste0(MainDir, "TSVs/"), full.names = TRUE, pattern = paste0(caller, ".mod.tsv"))
vep_files <- list.files(path = paste0(MainDir, "VEP/"), full.names = TRUE, pattern = paste0(caller, '.VEP.RefSeq.tsv')) %>%
    str_sort(numeric = TRUE)
patient_ID <- basename(vep_files) %>% str_replace_all(paste0(caller, '.VEP.RefSeq.tsv'), '')

# Load reference dbNSFP4.1 table with more meta data about genes:
xref_table <- read.table("~/Desktop/AAK_lab/Variant_Annotation/data/Annotation_DBs/dbNSFP4.1_short_gene_meta.txt",
                         sep = "\t", quote = "", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)

################################################################################
# Create excel workbook with short list of variants:
wb <- createWorkbook(
  creator = "Sergei_Zhuk",
  title = paste("SureSelect", file_ID),
  subject = file_ID,
  category = "Report"
)
callername <- strsplit(caller, '\\.')[[1]][2]
saveWorkbook(wb, file = paste0(file_ID, "_short_", callername, ".xlsx"), overwrite = TRUE)

# Write report tables into a workbook
for (file in vep_files){
  ID <- basename(file) %>% str_replace_all(paste0(caller, '.VEP.RefSeq.tsv'), '')
  print(paste('Processing', ID))
  print(paste('Loading VEP annotation from', file))
  vep_df <- read_tsv(file = file, col_names = TRUE, skip = 108)
  print(paste('Loading genotype info by using ID of', ID))
  gen_table <- read_tsv(file = paste0(MainDir, "TSVs/", ID, caller, '.mod.tsv')) %>%
    clean_names()
  vep_df <- vep_df %>%
    rename(uploaded_variation = `#Uploaded_variation`)
  vep_df$Location <- vep_df$Location %>% gsub("-.*","", .)
  vep_df$uploaded_variation <- paste0(vep_df$Location, '_', vep_df$REF_ALLELE, '/', vep_df$Allele) %>% str_replace(':', '_')
  comb_df <- merge(x = vep_df, y = gen_table, by = 'uploaded_variation', all.x = TRUE)
  comb_df <- merge(x =comb_df, y = xref_table, by.x = "SYMBOL", by.y = "RefGene", all.x = TRUE)
  comb_df$HGVS_desc <- paste(comb_df$SYMBOL, comb_df$HGVSg, comb_df$HGVSc, comb_df$HGVSp, sep = ";")
  comb_df$Consequence_AA <- str_extract(comb_df$HGVS_desc, regex("(?<=p\\.).*"))
  comb_df_reorder <- relocate(comb_df,
                              Chr = Location,
                              Ref = REF_ALLELE,
                              Alt = Allele,
                              RefGene = SYMBOL,
                              Exon = EXON,
                              HGVS_description = HGVS_desc,
                              Zyg = ZYG,
                              rsID = Existing_variation,
                              GATK_FILTER = filter,
                              CoverageDepth = ends_with('_dp'),
                              Consequence,
                              Consequence_AA,
                              gnomAD_exome_NFE = gnomADe_NFE_AF,
                              gnomAD_exome_Comb = gnomADe_AF,
                              gnomAD_genome_NFE = gnomAD_genomes_NFE_AF,
                              gnomAD_genome_Comb = gnomAD_genomes_AF, 
                              Thousand_EUR = EUR_AF,
                              TopMed_AF = TOPMED_GRCh38_20180418.vcf.gz_TOPMED,
                              CADD = CADD_PHRED,
                              VEP_Impact = IMPACT,
                              ClinVar_CLNSIG,
                              SIFT = SIFT,
                              PolyPhen = PolyPhen,
                              ClinPred = ClinPred_pred,
                              LOFTEE_LoF = LoF,
                              scSNV_ada_score_splicing = ada_score,
                              scSNV_rf_score_splicing = rf_score,
                              SpliceAI_cutoff = SpliceAI_cutoff,
                              ExAC_pLI,
                              gnomAD_pLI,
                              Full_Name = Gene_full_name,
                              Gene_Function = Function_description,
                              ClinVar_publications = ClinVar,
                              ClinVar_CLNREVSTAT,
                              ClinVar_CLNDN,
                              PubMed_records = PUBMED,
                              Orphanet_disorder,
                              Disease_description,
                              MIM_disease,
                              GO_biological_process,
                              GO_molecular_function,
                              HGVSg,
                              HGVSc,
                              HGVSp,
                              AlleleDepth = ends_with('_ad'),
                              GenotypeQual = ends_with('_gq'))
  cols.num <- c("gnomAD_exome_NFE", "CADD", "CoverageDepth", "GenotypeQual",
                "gnomAD_exome_Comb", "gnomAD_genome_NFE", "gnomAD_genome_Comb", "Thousand_EUR", "TopMed_AF")
  comb_df_reorder[cols.num] <- sapply(comb_df_reorder[cols.num], as.numeric)
  subset <- comb_df_reorder %>%
    filter(Consequence != "intron_variant" & (ClinVar_CLNSIG != "Pathogenic" | ClinVar_CLNSIG != "Pathogenic/Likely_pathogenic")) %>%
    filter(Consequence != "non_coding_transcript_exon_variant" & (ClinVar_CLNSIG != "Pathogenic" | ClinVar_CLNSIG != "Pathogenic/Likely_pathogenic")) %>%
    filter(Consequence != "intergenic_variant" & (ClinVar_CLNSIG != "Pathogenic" | ClinVar_CLNSIG != "Pathogenic/Likely_pathogenic")) %>%
    filter(Consequence != "intron_variant,NMD_transcript_variant" & (ClinVar_CLNSIG != "Pathogenic" | ClinVar_CLNSIG != "Pathogenic/Likely_pathogenic")) %>%
    filter(Consequence != "intron_variant,non_coding_transcript_variant" & (ClinVar_CLNSIG != "Pathogenic" | ClinVar_CLNSIG != "Pathogenic/Likely_pathogenic")) %>%
    filter(Consequence != "3_prime_UTR_variant" & (ClinVar_CLNSIG != "Pathogenic" | ClinVar_CLNSIG != "Pathogenic/Likely_pathogenic")) %>%
    filter(Consequence != "5_prime_UTR_variant" & (ClinVar_CLNSIG != "Pathogenic" | ClinVar_CLNSIG != "Pathogenic/Likely_pathogenic")) %>%
    filter(Consequence != "upstream_gene_variant" & (ClinVar_CLNSIG != "Pathogenic" | ClinVar_CLNSIG != "Pathogenic/Likely_pathogenic")) %>%
    filter(Consequence != "downstream_gene_variant" & (ClinVar_CLNSIG != "Pathogenic" | ClinVar_CLNSIG != "Pathogenic/Likely_pathogenic")) %>%
    filter(Consequence != 'synonymous_variant') %>%
    filter(Consequence != "synonymous_variant,NMD_transcript_variant") %>% 
    filter(GenotypeQual >= 20) %>%
    filter(gnomAD_exome_NFE < 0.01 | is.na(gnomAD_exome_NFE)) %>%
    filter(case_when(type == 'SNP' ~ CoverageDepth >= 5,
                     type == 'INDEL' ~ CoverageDepth >= 10)) %>%
    #filter(CANONICAL == 'YES') %>%
    arrange(!is.na(gnomAD_exome_NFE), gnomAD_exome_NFE)
  #subset[is.na(subset)] <- "."
  #subset[,cols.num] <- apply(subset[,cols.num], 2,
  #                           function(x) as.numeric(as.character(x)))
  print("Completed formatting table")
  name_sheet <- ID %>% str_extract("[:graph:]{1,31}")
  workbook <- loadWorkbook(paste0(file_ID, "_short_", callername, ".xlsx"))
  addWorksheet(workbook, name_sheet)
  writeData(workbook, name_sheet, subset, startCol = 1, startRow = 2, rowNames = FALSE,
            keepNA = TRUE, na.string = ".")
  saveWorkbook(workbook, paste0(file_ID, "_short_", callername, ".xlsx"), overwrite = TRUE)
  # write extended unfiltered dataframe in a separate files
  write.table(comb_df_reorder, file = paste0(MainDir, "Final/", ID, "_", callername, '.final.tsv'), sep = '\t', row.names = FALSE)
  print(paste('Finished processing of', ID))
}

################################################################################
# Make report with all variants 

################################################################################
# Create excel workbook with list of all found variants:
wb <- createWorkbook(
  creator = "Sergei_Zhuk",
  title = paste("SureSelect", file_ID),
  subject = file_ID,
  category = "Report"
)
saveWorkbook(wb, file = paste0(file_ID, "_all_", callername, ".xlsx"), overwrite = TRUE)

# Write report tables into a workbook
for (file in vep_files){
  ID <- basename(file) %>% str_replace_all(paste0(caller, '.VEP.RefSeq.tsv'), '')
  print(paste('Processing', ID))
  print(paste('Loading VEP annotation from', file))
  vep_df <- read_tsv(file = file, col_names = TRUE, skip = 108)
  print(paste('Loading genotype info by using ID of', ID))
  gen_table <- read_tsv(file = paste0(MainDir, "TSVs/", ID, caller, '.mod.tsv')) %>%
    clean_names()
  vep_df <- vep_df %>%
    rename(uploaded_variation = `#Uploaded_variation`)
  vep_df$Location <- vep_df$Location %>% gsub("-.*","", .)
  vep_df$uploaded_variation <- paste0(vep_df$Location, '_', vep_df$REF_ALLELE, '/', vep_df$Allele) %>% str_replace(':', '_')
  comb_df <- merge(x = vep_df, y = gen_table, by = 'uploaded_variation', all.x = TRUE)
  comb_df <- merge(x =comb_df, y = xref_table, by.x = "SYMBOL", by.y = "RefGene", all.x = TRUE)
  comb_df$HGVS_desc <- paste(comb_df$SYMBOL, comb_df$HGVSg, comb_df$HGVSc, comb_df$HGVSp, sep = ";")
  comb_df$Consequence_AA <- str_extract(comb_df$HGVS_desc, regex("(?<=p\\.).*"))
  comb_df_reorder <- relocate(comb_df,
                              Chr = Location,
                              Ref = REF_ALLELE,
                              Alt = Allele,
                              RefGene = SYMBOL,
                              Exon = EXON,
                              HGVS_description = HGVS_desc,
                              Zyg = ZYG,
                              rsID = Existing_variation,
                              GATK_FILTER = filter,
                              CoverageDepth = ends_with('_dp'),
                              Consequence,
                              Consequence_AA,
                              gnomAD_exome_NFE = gnomADe_NFE_AF,
                              gnomAD_exome_Comb = gnomADe_AF,
                              gnomAD_genome_NFE = gnomAD_genomes_NFE_AF,
                              gnomAD_genome_Comb = gnomAD_genomes_AF, 
                              Thousand_EUR = EUR_AF,
                              TopMed_AF = TOPMED_GRCh38_20180418.vcf.gz_TOPMED,
                              CADD = CADD_PHRED,
                              VEP_Impact = IMPACT,
                              ClinVar_CLNSIG,
                              SIFT = SIFT,
                              PolyPhen = PolyPhen,
                              ClinPred = ClinPred_pred,
                              LOFTEE_LoF = LoF,
                              scSNV_ada_score_splicing = ada_score,
                              scSNV_rf_score_splicing = rf_score,
                              SpliceAI_cutoff = SpliceAI_cutoff,
                              ExAC_pLI,
                              gnomAD_pLI,
                              Full_Name = Gene_full_name,
                              Gene_Function = Function_description,
                              ClinVar_publications = ClinVar,
                              ClinVar_CLNREVSTAT,
                              ClinVar_CLNDN,
                              PubMed_records = PUBMED,
                              Orphanet_disorder,
                              Disease_description,
                              MIM_disease,
                              GO_biological_process,
                              GO_molecular_function,
                              HGVSg,
                              HGVSc,
                              HGVSp,
                              AlleleDepth = ends_with('_ad'),
                              GenotypeQual = ends_with('_gq'))
  cols.num <- c("gnomAD_exome_NFE", "CADD", "CoverageDepth", "GenotypeQual",
                "gnomAD_exome_Comb", "gnomAD_genome_NFE", "gnomAD_genome_Comb", "Thousand_EUR", "TopMed_AF")
  comb_df_reorder[cols.num] <- sapply(comb_df_reorder[cols.num], as.numeric)
  subset <- comb_df_reorder %>%
    filter(GenotypeQual >= 20) %>%
    #filter(case_when(type == 'SNP' ~ CoverageDepth >= 5,
    #                 type == 'INDEL' ~ CoverageDepth >= 10)) %>%
    arrange(!is.na(gnomAD_exome_NFE), gnomAD_exome_NFE)
  #subset[is.na(subset)] <- "."
  #subset[,cols.num] <- apply(subset[,cols.num], 2,
  #                          function(x) as.numeric(as.character(x)))
  print("Completed formatting table")
  name_sheet <- ID %>% str_extract("[:graph:]{1,31}")
  workbook <- loadWorkbook(paste0(file_ID, "_all_", callername, ".xlsx"))
  addWorksheet(workbook, name_sheet)
  writeData(workbook, name_sheet, subset, startCol = 1, startRow = 2, rowNames = FALSE,
            keepNA = TRUE, na.string = ".")
  saveWorkbook(workbook, paste0(file_ID, "_all_", callername, ".xlsx"), overwrite = TRUE)
  # write extended unfiltered dataframe in a separate files
  write.table(comb_df_reorder, file = paste0(MainDir, "Final/", ID, "_", callername, '.final.tsv'), sep = '\t', row.names = FALSE)
  print(paste('Finished processing of', ID))
}
