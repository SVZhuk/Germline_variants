import argparse
import os

import pandas as pd


def process_tsv(file_path):
    # Read the TSV file
    df = pd.read_csv(file_path, sep="\t", header=0, dtype=str)

    # Add the 'uploaded_variation' column
    df["uploaded_variation"] = None

    # Process the 'uploaded_variation' column
    df_indel = df[df["TYPE"] == "INDEL"]
    df["uploaded_variation"] = df.apply(
        lambda row: f"{row['CHROM']}_{row['POS']}_{row['REF']}/{row['ALT']}", axis=1
    )
    df.loc[df["TYPE"] == "INDEL", "uploaded_variation"] = df_indel.apply(
        lambda row: f"{row['CHROM']}_{row['POS']}_{'-'}/{row['ALT'][1:]}"
        if len(row["REF"]) < len(row["ALT"])
        else f"{row['CHROM']}_{int(row['POS']) + 1}_{row['REF'][1:]}/-",
        axis=1,
    )

    # Select specific columns
    genot_df_trimmed = df[["TYPE"] + df.columns[5:10].tolist()]

    # Generate the output file name
    output_file_path = file_path.replace(".tsv", ".mod.tsv")
    if output_file_path == file_path:
        output_file_path = file_path + ".mod.tsv"

    # Write the modified DataFrame to a new TSV file
    genot_df_trimmed.to_csv(output_file_path, sep="\t", index=False)

    print(f"Done processing {os.path.basename(file_path)}")
    print(f"Output saved as {output_file_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Process a TSV file and produce a modified version."
    )
    parser.add_argument(
        "--file_path", type=str, required=True, help="Path to the input TSV file"
    )
    args = parser.parse_args()

    process_tsv(args.file_path)


if __name__ == "__main__":
    main()
