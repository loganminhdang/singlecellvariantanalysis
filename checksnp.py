#!/usr/bin/env python3
import pysam
import argparse
from collections import Counter

def parse_position(pos_str):
    """Parses a genomic position string like 'chr17:29032660'."""
    try:
        chrom, pos = pos_str.split(':')
        return chrom, int(pos)
    except ValueError:
        raise ValueError("Position must be in format 'chr:position'")

def check_snp_pileup(args):
    """Performs a pileup at a specific genomic coordinate to check for SNP evidence."""
    chrom, pos = parse_position(args.position)
    # pysam uses 0-based indexing, so we subtract 1
    pos_0_based = pos - 1 

    print(f"--- Checking SNP at {chrom}:{pos} ---")
    
    base_counts = Counter()
    total_reads_at_site = 0
    
    with pysam.AlignmentFile(args.bam, "rb") as bamfile:
        # mpileup iterates over columns of the alignment
        for pileupcolumn in bamfile.pileup(chrom, pos_0_based, pos_0_based + 1, truncate=True):
            # We are only interested in our specific single position
            if pileupcolumn.pos == pos_0_based:
                total_reads_at_site = pileupcolumn.nsegments
                # Iterate through each read in the pileup at this position
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        # Get the base from the read
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                        base_counts[base] += 1

    if total_reads_at_site == 0:
        print("No reads found covering this specific position.")
        return

    print(f"\nFound {total_reads_at_site} reads covering the site.")
    print("Base counts from reads:")
    print("Base\tCount\tPercentage")
    print("----\t-----\t----------")
    for base, count in base_counts.most_common():
        percentage = (count / total_reads_at_site) * 100
        print(f"{base}\t{count}\t{percentage:.2f}%")

def main():
    parser = argparse.ArgumentParser(description="Perform a simple base pileup at a specific genomic position to investigate a known SNP.")
    parser.add_argument("--bam", required=True, help="Path to the indexed input BAM file.")
    parser.add_argument("--position", required=True, help="Genomic position of the SNP to check, e.g., 'chr17:29032660'.")
    args = parser.parse_args()
    check_snp_pileup(args)

if __name__ == "__main__":
    main()