#!/usr/bin/env python

# uses argparse to read bed file and extract gene names:
# python genes_from_bed.py -i <input.bed> -o <output.txt>

import argparse
import os
from typing import List, Set, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract set of gene names from bed file"
    )
    parser.add_argument("-i", "--input", required=True, help="Path to input bed file")
    parser.add_argument("-o", "--output", required=True, help="Path to output file")
    return parser.parse_args()


def extract_genes(input_bed: str) -> Tuple[Set[str], int]:
    genes = set()
    exclude_prefixes: List[str] = [
        "ENST",
        "ens|",
        "miRNA|",
        "ccds|",
        "NR_",
        "NM_",
        "LOC",
        "LINC",
    ]
    with open(input_bed, "r", encoding="utf-8") as bed_file:
        for line in bed_file:
            fields: list[str] = line.strip().split("\t")
            if len(fields) > 3:
                gene_names: list[str] = fields[3].split(",")
                for gene_name in gene_names:
                    gene_name = gene_name.strip()
                    if gene_name.startswith("ref|"):
                        gene_name = gene_name.split("|")[1]
                    if ":" not in gene_name and not any(
                        gene_name.startswith(prefix) for prefix in exclude_prefixes
                    ):
                        genes.add(gene_name)
    gene_count = len(genes)
    return genes, gene_count


def write_output(output_txt: str, genes: Set[str], gene_count: int) -> None:
    output_file_name = os.path.splitext(os.path.basename(output_txt))[0]
    with open(output_txt, "w", encoding="utf-8") as output_file:
        output_file.write(f"Contents of {output_file_name} file:\n")
        output_file.write("List of genes:\n")
        output_file.write(", ".join(sorted(genes)) + "\n")
        output_file.write("\nTotal number of genes:\n")
        output_file.write(f"{gene_count}\n")

        print(f"Contents of {output_file_name} file:\n")
        print("List of genes:")
        print(", ".join(sorted(genes)))
        print("\nTotal number of genes:")
        print(f"{gene_count}")


def main() -> None:
    args = parse_args()
    genes, gene_count = extract_genes(args.input)
    write_output(args.output, genes, gene_count)


if __name__ == "__main__":
    main()
