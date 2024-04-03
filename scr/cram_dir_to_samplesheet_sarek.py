#!/usr/bin/env python3

import os
import csv
import argparse

def parse_args(args=None) -> argparse.Namespace:
    Description = "Generate nf-core/sarek samplesheet from a directory of CRAM files."
    Epilog = "Example usage: python cram_dir_to_samplesheet.py <CRAM_DIR> <SAMPLESHEET_FILE>"
    
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("CRAM_DIR", help="Folder containing raw CRAM files.")
    parser.add_argument("SAMPLESHEET_FILE", help="Output samplesheet file.")
    
    return parser.parse_args(args)
    
def cram_dir_to_samplesheet(cram_dir: str, samplesheet_file: str) -> None:
    files = []
    for root, dirs, filenames in os.walk(cram_dir):
        for filename in filenames:
            if filename.endswith(".cram"):
                patient = sample = filename.split(".")[0]
                cram_full_path = os.path.join(root, filename)
                crai_full_path = cram_full_path + ".crai"
                files.append({"patient": patient, "sample": sample, "cram": cram_full_path, "crai": crai_full_path})
                
    with open(samplesheet_file, "w") as fh:
        fieldnames = ["patient", "sample", "cram", "crai"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter=",")
        
        writer.writeheader()
        for file in files:
            cram = os.path.join(cram_dir, file["cram"])
            crai = os.path.join(cram_dir, file["cram"] + ".crai")
            writer.writerow({"patient": file["patient"], "sample": file["sample"], "cram": cram, "crai": crai})
            
if __name__ == "__main__":
    args = parse_args()
    cram_dir_to_samplesheet(args.CRAM_DIR, args.SAMPLESHEET_FILE)
