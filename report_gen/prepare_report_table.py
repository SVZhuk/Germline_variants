import argparse
import re

import janitor  # noqa: F401
import pandas as pd


def ensure_string_type(df: pd.DataFrame, column: str) -> pd.DataFrame:
    if column in df.columns:
        df[column] = df[column].astype(str).str.strip()
    return df


def read_vep_table(vep_table: str) -> pd.DataFrame:
    df = pd.read_csv(vep_table, sep="\t", header=0, skiprows=104, index_col=False)
    df = df.rename(columns={"#Uploaded_variation": "uploaded_variation"})
    # VEP uses rsID instead of coordinates for some variants
    # for merge with genotype data, we need to extract the coordinates
    df["Location"] = df["Location"].str.replace(r"-.*", "", regex=True)
    df["uploaded_variation"] = (
        df["Location"]
        .str.cat(df["REF_ALLELE"], sep="_")
        .str.cat(df["Allele"], sep="/")
        .str.replace(":", "_")
    )
    df.columns = df.columns.str.lower()

    return df


def read_expression_table(expression_table: str) -> pd.DataFrame:
    df = pd.read_excel(expression_table, header=0, index_col=False)
    selected_df = df[["description", "rank_gse71613"]]
    selected_df.columns = selected_df.columns.str.lower()
    return selected_df


def read_genotype_table(genotype_table: str) -> pd.DataFrame:
    df = pd.read_csv(genotype_table, sep="\t", header=0, index_col=False)
    df.columns = df.columns.str.lower()

    return df


def read_gene_metadata(metadata: str) -> pd.DataFrame:
    df = pd.read_csv(
        metadata,
        sep="\t",
        quotechar='"',
        header=0,
        dtype=str,
        na_filter=False,
        index_col=False,
    )
    df.columns = df.columns.str.lower()

    return df


def merge_dataframes(
    vep_df: pd.DataFrame,
    genotype_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    expression_df: pd.DataFrame,
) -> pd.DataFrame:
    vep_df = ensure_string_type(vep_df, "uploaded_variation")
    genotype_df = ensure_string_type(genotype_df, "uploaded_variation")
    comb_df = pd.merge(vep_df, genotype_df, on="uploaded_variation", how="left")
    comb_df = pd.merge(
        comb_df, metadata_df, left_on="symbol", right_on="refgene", how="left"
    )
    comb_df = pd.merge(
        comb_df, expression_df, left_on="refgene", right_on="description", how="left"
    )

    return comb_df


def handle_duplicate_columns(df: pd.DataFrame) -> pd.DataFrame:
    duplicate_columns = df.columns[df.columns.duplicated()].unique()

    for col in duplicate_columns:
        cols = df.loc[:, df.columns == col]
        if cols.nunique(axis=1).eq(1).all():
            df = df.loc[:, ~df.columns.duplicated()]
        else:
            for i in range(1, len(cols.columns)):
                df.rename(
                    columns={cols.columns[i]: f"{cols.columns[i]}_{i}"}, inplace=True
                )

    return df


def rename_and_reorder_columns(comb_df: pd.DataFrame) -> pd.DataFrame:
    comb_df = comb_df.rename(
        columns={
            "location": "Chr",
            "ref_allele": "Ref",
            "allele": "Alt",
            "symbol": "RefGene",
            "exon": "Exon",
            "hgvs_desc": "HGVS_description",
            "zyg": "Zyg",
            "existing_variation": "rsID",
            "filter": "GATK_FILTER",
            "gnomade_nfe_af": "gnomAD_exome_NFE",
            "gnomade_af": "gnomAD_exome_Comb",
            "gnomad_genomes_nfe_af": "gnomAD_genome_NFE",
            "gnomad_genomes_af": "gnomAD_genome_Comb",
            "eur_af": "Thousand_EUR",
            "topmed_grch38_20180418.vcf.gz_topmed": "TopMed_AF",
            "cadd_phred": "CADD",
            "impact": "VEP_Impact",
            "clinpred_pred": "ClinPred",
            "lof": "LOFTEE_LoF",
            "ada_score": "scSNV_ada_score_splicing",
            "rf_score": "scSNV_rf_score_splicing",
            "spliceai_cutoff": "SpliceAI_cutoff",
            "gene_full_name": "Full_Name",
            "function_description": "Gene_Function",
            "clinvar": "ClinVar_publications",
            "pubmed": "PubMed_records",
        }
    )

    # Rename columns that end with '_ad', '_gq', and '_dp'
    comb_df = comb_df.rename(
        columns=lambda x: re.sub(r"_ad$", "AlleleDepth", x, flags=re.IGNORECASE)
    )
    comb_df = comb_df.rename(
        columns=lambda x: re.sub(r"_gq$", "GenotypeQual", x, flags=re.IGNORECASE)
    )
    comb_df = comb_df.rename(
        columns=lambda x: re.sub(r"_dp$", "CoverageDepth", x, flags=re.IGNORECASE)
    )

    comb_df = handle_duplicate_columns(comb_df)

    # Reorder columns
    columns_order = [
        "Chr",
        "Ref",
        "Alt",
        "RefGene",
        "Exon",
        "HGVS_description",
        "Zyg",
        "rsID",
        "GATK_FILTER",
        "Consequence",
        "Consequence_AA",
        "gnomAD_exome_NFE",
        "gnomAD_exome_Comb",
        "gnomAD_genome_NFE",
        "gnomAD_genome_Comb",
        "Thousand_EUR",
        "TopMed_AF",
        "CADD",
        "VEP_Impact",
        "ClinVar_CLNSIG",
        "SIFT",
        "PolyPhen",
        "ClinPred",
        "LOFTEE_LoF",
        "scSNV_ada_score_splicing",
        "scSNV_rf_score_splicing",
        "SpliceAI_cutoff",
        "ExAC_pLI",
        "gnomAD_pLI",
        "Full_Name",
        "Gene_Function",
        "ClinVar_publications",
        "ClinVar_CLNREVSTAT",
        "ClinVar_CLNDN",
        "PubMed_records",
        "Orphanet_disorder",
        "Disease_description",
        "MIM_disease",
        "GO_biological_process",
        "GO_molecular_function",
        "HGVSg",
        "HGVSc",
        "HGVSp",
        "CoverageDepth",
        "AlleleDepth",
        "GenotypeQual",
    ]

    # Ensure there are no duplicate columns in columns_order
    columns_order = list(dict.fromkeys(columns_order))
    print(f"Reordered columns: {columns_order}")

    # Reindex DataFrame
    comb_df_reorder = comb_df.reindex(columns=columns_order)

    # Convert specified columns to numeric
    cols_num = [
        "gnomAD_exome_NFE",
        "CADD",
        "CoverageDepth",
        "GenotypeQual",
        "gnomAD_exome_Comb",
        "gnomAD_genome_NFE",
        "gnomAD_genome_Comb",
        "Thousand_EUR",
        "TopMed_AF",
    ]
    comb_df_reorder[cols_num] = comb_df_reorder[cols_num].apply(
        pd.to_numeric, errors="coerce"
    )

    return comb_df_reorder


def format_report_table(comb_df: pd.DataFrame) -> pd.DataFrame:
    comb_df["hgvs_desc"] = (
        comb_df[["symbol", "hgvsg", "hgvsc", "hgvsp"]].astype(str).agg(";".join, axis=1)
    )
    comb_df["consequence_aa"] = comb_df["hgvs_desc"].str.extract(r"(?<=p\.)(.*)")

    return comb_df


def main(vep_table, genotype_table, metadata, expression_table, output):
    vep_df = read_vep_table(vep_table)
    genotype_df = read_genotype_table(genotype_table)
    metadata_df = read_gene_metadata(metadata)
    expression_df = read_expression_table(expression_table)

    total_df = merge_dataframes(vep_df, genotype_df, metadata_df, expression_df)
    total_df = format_report_table(total_df)
    total_df.to_csv(output, sep="\t", index=False)
    print(f"Report saved to {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Prepare a report table from VEP and genotype data"
    )
    parser.add_argument(
        "--vep_table", required=True, help="Path to the VEP annotation table"
    )
    parser.add_argument(
        "--genotype_table", required=True, help="Path to the adjusted genotype table"
    )
    parser.add_argument("--metadata", required=True, help="Path to the metadata file")
    parser.add_argument(
        "--expression_table", required=True, help="Path to the expression data table"
    )
    parser.add_argument(
        "--output", required=True, help="Path to save the output report"
    )

    args = parser.parse_args()

    main(
        args.vep_table,
        args.genotype_table,
        args.metadata,
        args.expression_table,
        args.output,
    )
