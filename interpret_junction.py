#!/usr/bin/env python3
import pysam
import os
import argparse
from collections import defaultdict, Counter

def load_annotated_barcodes(file_path):
    """
    Loads barcodes and their cell type annotations from a two-column CSV file.
    It automatically handles and strips common suffixes like '-1' or '_1' from the barcodes.
    Returns a dictionary mapping the barcode -> cell_type.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Error: Annotation file not found at '{file_path}'")

    annotations = {}
    with open(file_path, 'r') as f:
        # Verify the header structure to ensure the CSV is formatted correctly
        header = f.readline().strip().lower().split(',')
        if header != ['barcode', 'cell_type']:
            raise ValueError("Annotation file must be a CSV with a header line: 'barcode,cell_type'")
            
        for line in f:
            parts = line.strip().split(',')
            if len(parts) == 2:
                full_barcode, cell_type = parts
                # Clean the barcode to match the BAM format (e.g., removing '-1')
                simple_barcode = full_barcode.split('-')[0].split('_')[0]
                annotations[simple_barcode] = cell_type
                
    print(f"Loaded annotations for {len(annotations)} unique barcodes from {file_path}")
    return annotations

def analyze_annotated_cells(args):
    """
    Analyzes a specific list of high-quality, annotated cells to quantify junctions
    and determine zygosity, breaking down results by cell type.
    """
    print("\n--- Running Analysis on Annotated Cell List ---")
    
    # --- 1. Load and Parse Inputs ---
    # Load the cell metadata (type and barcode) and parse the genomic region string
    annotated_barcodes = load_annotated_barcodes(args.barcodes)
    chrom, start, end = args.region.split(':')[0], int(args.region.split(':')[1].split('-')[0]), int(args.region.split(':')[1].split('-')[1])

    # Convert the user-provided junction list into a dictionary for quick lookup
    junction_definitions = {}
    junction_map = {}
    for j in args.junctions:
        name, coords = j.split(':')
        j_start, j_end = map(int, coords.split('-'))
        junction_tuple = (j_start, j_end)
        junction_definitions[name] = junction_tuple
        junction_map[junction_tuple] = name
    
    print("\nAnalyzing the following junctions:")
    for name, coords in junction_definitions.items():
        print(f"  - {name}: {coords[0]}-{coords[1]}")

    # --- 2. Scan BAM and Collect Data ---
    print("\nScanning BAM file for matching reads...")
    cell_alleles = defaultdict(set) # Stores unique junctions found per cell
    with pysam.AlignmentFile(args.bam, "rb") as bamfile:
        # Efficiently fetch only reads in the gene region
        for read in bamfile.fetch(chrom, start, end):
            # Skip reads without a cell barcode (CB) tag
            if not read.has_tag('CB'):
                continue

            # Standardize the barcode from the BAM to match our annotation list
            full_barcode_from_bam = read.get_tag('CB')
            simple_barcode = full_barcode_from_bam.split('-')[0].split('_')[0]
            
            # Skip cells that were not included in the high-quality annotation list
            if simple_barcode not in annotated_barcodes:
                continue

            # Check if the read is 'spliced' (has multiple alignment blocks)
            blocks = read.get_blocks()
            if len(blocks) < 2:
                continue

            # Inspect the gaps between blocks to see if they match our target junctions
            for i in range(len(blocks) - 1):
                junction_tuple = (blocks[i][1], blocks[i+1][0])
                if junction_tuple in junction_map:
                    junction_name = junction_map[junction_tuple]
                    # Record the presence of this junction/allele for this cell
                    cell_alleles[simple_barcode].add(junction_name)

    # Exit if no relevant data was found
    if not cell_alleles:
        print("\nNo cells from your list were found to express any of the target junctions.")
        return
        
    # --- 3. Generate Final, Detailed Report ---
    report_lines = [] # Collects lines for printing/saving
    
    # Organize data: Genetic State -> Cell Type -> Frequency
    report = defaultdict(Counter)
    for cell, alleles_found in cell_alleles.items():
        cell_type = annotated_barcodes.get(cell, "Unknown")
        # Create a sorted key like "Deletion_B + Wild-Type" for heterozygotes
        state_key = " + ".join(sorted(list(alleles_found)))
        report[state_key][cell_type] += 1

    total_expressing_cells = len(cell_alleles)
    
    report_lines.append(f"--- Zygosity and Cell Type Report ---")
    report_lines.append(f"Found {total_expressing_cells} cells expressing any target allele within your list.")
    report_lines.append("--------------------------------------------------------------------------")
    
    # Sort the report categories by the total number of cells found in each state
    sorted_states = sorted(report.items(), key=lambda item: sum(item[1].values()), reverse=True)

    for state, cell_type_counts in sorted_states:
        total_in_state = sum(cell_type_counts.values())
        report_lines.append(f"\nState: {state} ({total_in_state} total cells)")
        report_lines.append("  |--Cell Type Breakdown:")
        
        # Breakdown the genetic state by biological cell type
        for cell_type, count in cell_type_counts.most_common():
            percentage = (count / total_in_state) * 100
            report_lines.append(f"    - {cell_type}: {count} cells ({percentage:.1f}%)")
            
    report_lines.append("--------------------------------------------------------------------------")

    # --- 4. Print to Terminal and Save to File ---
    final_report = "\n".join(report_lines)
    
    print("\n" + final_report) # Display the report in the console
    
    # Optionally save the result to the disk if an output file was specified
    if args.output_file:
        with open(args.output_file, 'w') as f:
            f.write(final_report)
        print(f"\nReport successfully saved to: {args.output_file}")


def main():
    # Setup the command-line argument parser
    parser = argparse.ArgumentParser(description="Analyze splice junctions for a specific list of annotated cells and break down results by cell type.")
    parser.add_argument("--bam", required=True, help="Path to the indexed input BAM file.")
    parser.add_argument("--region", required=True, help="Genomic region to analyze, e.g., 'chr17:29032000-29042500'.")
    parser.add_argument("--barcodes", required=True, help="Path to a two-column CSV file with a header: 'barcode,cell_type'.")
    parser.add_argument("--junctions", required=True, nargs='+', help="List of junctions to quantify, format: 'NAME:START-END'.")
    parser.add_argument("--output-file", "-o", help="(Optional) Path to save the final report text file.")
    
    args = parser.parse_args()
    analyze_annotated_cells(args)

if __name__ == "__main__":
    main()