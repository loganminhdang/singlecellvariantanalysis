#!/usr/bin/env python3
"""
Splicing and Zygosity Analyzer for scRNA-seq
-------------------------------------------
This script processes 10x Genomics BAM files to identify and quantify 
splice junctions. It helps distinguish between wild-type and variant alleles (e.g., deletions or insertions)
at the single-cell level.
"""

import pysam
import os
import argparse
from collections import defaultdict, Counter
import time

def load_barcodes_from_file(file_path):
    """
    Reads a list of cell barcodes.
    Standardizes them by removing 10x suffixes like '-1'.
    """
    if not os.path.exists(file_path):
        print(f"Error: Barcode file not found at '{file_path}'")
        return None
    with open(file_path, 'r') as f:
        # We use a set for 'O(1)' lookups, making it very fast to check if 
        # a read belongs to a high-quality cell.
        barcodes = {line.strip().split('-')[0] for line in f if line.strip()}
    print(f"Loaded {len(barcodes)} unique barcode stems from {file_path}")
    return barcodes

def parse_region(region_str):
    """
    Converts a string like 'chr17:29032000-29042500' into machine-readable integers.
    """
    try:
        chrom, coords = region_str.split(':')
        start, end = map(int, coords.split('-'))
        return chrom, start, end
    except ValueError:
        raise ValueError("Region must be in format 'chr:start-end'")

def discover_junctions(args):
    """
    MODE: DISCOVER
    Scans the gene for ANY gap in alignment. Useful for finding the exact 
    coordinates of a CRISPR deletion or a novel splice variant.
    """
    print("--- Running in Discovery Mode ---")
    chrom, start, end = parse_region(args.region)
    target_barcodes = load_barcodes_from_file(args.barcodes) if args.barcodes else None

    print(f"Scanning region {args.region} for junctions with a gap size >= {args.min_gap_size} bp...")
    
    junctions = []
    # Open the BAM file. 'rb' means Read Binary.
    with pysam.AlignmentFile(args.bam, "rb") as bamfile:
        # fetch() jump-reads only the relevant part (the coordinates being queried) of the indexed BAM, saving time.
        for read in bamfile.fetch(chrom, start, end):
            # Skip low-quality, duplicate, or unmapped reads.
            if read.is_unmapped or read.is_duplicate or read.is_secondary:
                continue
            
            # If a barcode list was provided, ignore reads from cells not in that list.
            if target_barcodes:
                if not read.has_tag('CB') or read.get_tag('CB').split('-')[0] not in target_barcodes:
                    continue
            
            # get_blocks() returns the physical segments of the read that aligned to the genome.
            # A 'spliced' read will have 2 or more blocks with a gap between them.
            blocks = read.get_blocks()
            if len(blocks) < 2:
                continue

            # Check the gaps between the aligned blocks.
            for i in range(len(blocks) - 1):
                gap_start = blocks[i][1] # End of the current block
                gap_end = blocks[i+1][0] # Start of the next block
                gap_size = gap_end - gap_start

                # If the gap is large enough (users can define the minimum gap size), we classify it as an intron or deletion.
                if gap_size >= args.min_gap_size:
                    junctions.append((chrom, gap_start, gap_end, gap_size))
    
    if not junctions:
        print("No junctions found matching the criteria.")
        return

    # Count how many times each unique junction appeared across all reads.
    junction_counts = Counter(junctions)
    print(f"\nFound {len(junction_counts)} unique junctions.")
    
    # Save results to a TSV.
    with open(args.output_file, 'w') as f:
        f.write("Count\tChromosome\tJunction_Start\tJunction_End\tGap_Size\n")
        # .most_common() sorts by frequency (highest counts first).
        for (chrom, j_start, j_end, gap_size), count in junction_counts.most_common():
            f.write(f"{count}\t{chrom}\t{j_start}\t{j_end}\t{gap_size}\n")
    print(f"Results saved to {args.output_file}")

def quantify_junctions(args):
    """
    MODE: QUANTIFY
    Looks for specific junctions (e.g., Wild-Type vs Deletion) and determines
    which alleles are present in every individual cell.
    """
    print("--- Running in Quantification Mode ---")
    chrom, start, end = parse_region(args.region)
    target_barcodes = load_barcodes_from_file(args.barcodes) if args.barcodes else None

    # Maps coordinates to a human-readable name (e.g., (100, 500) -> "Exon2_3_Del")
    junction_definitions = {}
    junction_map = {}
    for j in args.junctions:
        try:
            name, coords = j.split(':')
            j_start, j_end = map(int, coords.split('-'))
            junction_tuple = (j_start, j_end)
            junction_definitions[name] = junction_tuple
            junction_map[junction_tuple] = name
        except ValueError:
            raise ValueError("Junctions must be in format 'NAME:START-END'")
    
    print("Analyzing the following junctions:")
    for name, coords in junction_definitions.items():
        print(f"  - {name}: {coords[0]}-{coords[1]}")

    # Key = Cell Barcode, Value = Set of allele names found in that cell.
    cell_alleles = defaultdict(set)
    
    with pysam.AlignmentFile(args.bam, "rb") as bamfile:
        for read in bamfile.fetch(chrom, start, end):
            # 10x BAMs store the cell barcode in the 'CB' tag.
            if not read.has_tag('CB'):
                continue

            simple_barcode = read.get_tag('CB').split('-')[0]
            if target_barcodes and simple_barcode not in target_barcodes:
                continue

            blocks = read.get_blocks()
            if len(blocks) < 2:
                continue

            for i in range(len(blocks) - 1):
                junction_tuple = (blocks[i][1], blocks[i+1][0])
                # If this read spans one of our target junctions, record it for this cell.
                if junction_tuple in junction_map:
                    junction_name = junction_map[junction_tuple]
                    cell_alleles[simple_barcode].add(junction_name)

    if not cell_alleles:
        print("No cells found expressing any of the target junctions.")
        return

    # Count the 'States'. e.g. How many cells have 'WT' only vs 'WT + Del'?
    report = Counter()
    for cell, alleles_found in cell_alleles.items():
        # Sort names so 'WT + Del' is the same as 'Del + WT'
        state_key = " + ".join(sorted(list(alleles_found)))
        report[state_key] += 1

    total_cells = len(cell_alleles)
    print("\n--- Zygosity Report ---")
    print(f"Total cells expressing any target allele: {total_cells}")
    print("----------------------------------------------------------")
    print("Count\t\tAllele State(s) Detected in Cell")
    print("-----\t\t--------------------------------")
    for state, count in report.most_common():
        print(f"{count}\t\t{state}")
    print("----------------------------------------------------------")

def main():
    """
    Defines the Command Line Interface (CLI).
    """
    parser = argparse.ArgumentParser(description="Splice junction analyzer for scRNA-seq.")
    parser.add_argument("--bam", required=True, help="Input BAM file.")
    parser.add_argument("--region", required=True, help="Region e.g. 'chr17:29032000-29042500'.")
    parser.add_argument("--barcodes", help="Optional: List of high-quality barcodes.")
    
    # Subparsers allow us to have two different 'modes' with different options.
    subparsers = parser.add_subparsers(dest="mode", required=True, help="Operating mode")

    # Options for 'discover'
    parser_discover = subparsers.add_parser("discover", help="Find all junctions.")
    parser_discover.add_argument("--min-gap-size", type=int, default=200, help="Min intron length.")
    parser_discover.add_argument("--output-file", default="discovered_junctions.tsv")
    parser_discover.set_defaults(func=discover_junctions)

    # Options for 'quantify'
    parser_quantify = subparsers.add_parser("quantify", help="Count specific junctions.")
    parser_quantify.add_argument("--junctions", required=True, nargs='+', help="Format 'NAME:START-END'")
    parser_quantify.set_defaults(func=quantify_junctions)

    args = parser.parse_args()
    # Execute the function associated with the chosen mode.
    args.func(args)

if __name__ == "__main__":
    main()