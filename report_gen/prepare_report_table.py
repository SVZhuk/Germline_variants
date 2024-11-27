import argparse

import janitor  # noqa: F401
import pandas as pd


def read_vep_table(vep_table: str) -> pd.DataFrame:
    df = pd.read_csv(vep_table, sep="\t", header=0, skiprows=104, index_col=False)
    df = df.rename(columns={"#Uploaded_variation": "uploaded_variation"})
    return df


def read_expression_table(expression_table: str) -> pd.DataFrame:
    df = pd.read_excel(expression_table, header=0, index_col=False)
    selected_df = df[["description", "rank_gse71613"]]
    return selected_df


def read_genotype_table(genotype_table: str) -> pd.DataFrame:
    df = pd.read_csv(genotype_table, sep="\t", header=0, index_col=False)
    df = df.clean_names()
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
    return df


def merge_dataframes(
    vep_df: pd.DataFrame,
    genotype_df: pd.DataFrame,
    metadata_df: pd.DataFrame,
    expression_df: pd.DataFrame,
) -> pd.DataFrame:
    # Merge VEP and genotype data
    comb_df = pd.merge(vep_df, genotype_df, on="uploaded_variation", how="left")
    comb_df = pd.merge(
        comb_df, metadata_df, left_on="SYMBOL", right_on="RefGene", how="left"
    )
    comb_df = pd.merge(
        comb_df, expression_df, left_on="RefGene", right_on="description", how="left"
    )

    return comb_df


def format_report_table(comb_df: pd.DataFrame) -> pd.DataFrame:
    comb_df["Location"] = comb_df["Location"].str.replace(r"-.*", "", regex=True)
    comb_df["uploaded_variation"] = (
        comb_df["Location"]
        .str.cat(comb_df["REF_ALLELE"], sep="_")
        .str.cat(comb_df["Allele"], sep="/")
        .str.replace(":", "_")
    )
    comb_df['HGVS_desc'] = comb_df[['SYMBOL', 'HGVSg', 'HGVSc', 'HGVSp']].astype(str).agg(';'.join, axis=1)
    comb_df['Consequence_AA'] = comb_df['HGVS_desc'].str.extract(r'(?<=p\.)(.*)')

    return comb_df


def main(vep_table, genotype_table, metadata, expression_table, output):
    # Read input files
    vep_df = read_vep_table(vep_table)
    genotype_df = read_genotype_table(genotype_table)
    metadata_df = read_gene_metadata(metadata)
    expression_df = read_expression_table(expression_table)

    total_df = merge_dataframes(vep_df, genotype_df, metadata_df, expression_df)
    total_df = format_report_table(total_df)

    # Add any additional columns or calculations here

    # Write the report
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
