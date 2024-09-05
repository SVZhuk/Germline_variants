#!/usr/bin/env python3

import os
import re
import argparse
import sys
import pathlib
import csv

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="List all files in a directory.")
    parser.add_argument("directory", help="Path to the directory.")
    return parser.parse_args()

def sanitize_sample(path: str, extension: str) -> dict:
    file_name = os.path.basename(path).replace(extension, "")
    id_pattern = r"^(.+?)_S\d+_(L\d{3})_(R[12])_\d+$"
    match = re.match(id_pattern, file_name)
    if match:
        sample_id, lane, read = match.groups()
        return {"ID": sample_id, 
                "lane": lane, 
                "read": read, 
                "file_path": path}
    else:
        print(f"Warning: Coult not extract infirmation from file name: {file_name}")
        return {}

def list_files(directory: str) -> list:
    files_list = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            if file.endswith(".fastq.gz"):
                files_list.append(file_path)
    sorted_files_list = sorted(files_list)
    for file in sorted_files_list:
        sample_info = sanitize_sample(file, ".fastq.gz")
        print(sample_info)

        
def make_sample_sheet(samples: list):
    csv_header = ["ID", "lane", "R1", "R2"]
    print(",".join(csv_header))

if __name__ == "__main__":
    args: argparse.Namespace = parse_args()
    list_samples = list_files(args.directory)
    make_sample_sheet(list_samples)
    print(f"Python version: {sys.version}")