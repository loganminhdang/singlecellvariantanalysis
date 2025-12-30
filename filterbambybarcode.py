#!/usr/bin/env python3
import pysam
import os
import argparse
import time
import subprocess
import shutil
from concurrent.futures import ProcessPoolExecutor

def load_barcodes_from_file(file_path):
    """
    Reads a list of target barcodes from a file. 
    It standardizes them by stripping 10x suffixes (e.g., -1 or _1) 
    so that barcodes from different software (CellRanger vs Seurat) can be matched.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Error: Barcode file not found at '{file_path}'")
    with open(file_path, 'r') as f:
        # Create a set of barcode 'stems' for O(1) lightning-fast lookup.
        barcodes = {line.strip().split('-')[0].split('_')[0] for line in f if line.strip()}
    print(f"Loaded {len(barcodes)} unique barcode stems to filter for.")
    return barcodes

def filter_chunk(args_tuple):
    """
    WORKER FUNCTION: Each CPU core runs this function on a different chromosome.
    This is what makes the script fast—it processes the genome in parallel chunks.
    """
    input_bam, temp_dir, chromosome, target_barcode_stems = args_tuple
    
    temp_output_bam = os.path.join(temp_dir, f"{chromosome}.bam")
    read_count = 0

    try:
        # Open the master BAM and create a temporary mini-BAM for this chromosome
        with pysam.AlignmentFile(input_bam, "rb") as infile, \
             pysam.AlignmentFile(temp_output_bam, "wb", header=infile.header) as outfile:
            
            # Fetch reads only for the assigned chromosome
            for read in infile.fetch(chromosome):
                if read.has_tag('CB'):
                    # Extract the barcode from the read and standardize it
                    full_barcode_from_bam = read.get_tag('CB')
                    simple_barcode = full_barcode_from_bam.split('-')[0].split('_')[0]
                    
                    # If this read belongs to one of our target cells, save it to the mini-BAM
                    if simple_barcode in target_barcode_stems:
                        outfile.write(read)
                        read_count += 1
                        
        # If the chromosome actually had relevant reads, index the mini-BAM
        if read_count > 0:
            pysam.index(temp_output_bam)
            return temp_output_bam
        else:
            # Otherwise, delete the empty file to save space
            if os.path.exists(temp_output_bam):
                os.remove(temp_output_bam)
            return None

    except Exception as e:
        print(f"Error processing chromosome {chromosome}: {e}")
        return None

def main():
    """
    ORCHESTRATOR: Manages the parallel execution and the final file merge.
    """
    parser = argparse.ArgumentParser(
        description="Filter a 10x BAM file in parallel to keep only reads from a specified list of cell barcodes.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    # Define required user inputs
    parser.add_argument("--input-bam", "-i", required=True, help="Path to the original BAM.")
    parser.add_argument("--barcodes", "-b", required=True, help="Path to text file with cell barcodes.")
    parser.add_argument("--output-bam", "-o", required=True, help="Path for the filtered output BAM.")
    parser.add_argument("--processes", "-p", type=int, default=os.cpu_count(), help="Number of CPU cores to use.")
    
    args = parser.parse_args()
    start_time = time.time()
    
    # Ensure 'samtools' is available for the final high-speed merge
    if not shutil.which("samtools"):
        print("Error: samtools is required for the final merge step.")
        return

    print(f"Using {args.processes} processes for parallel execution.")
    target_barcodes = load_barcodes_from_file(args.barcodes)

    # Create a workspace for the temporary per-chromosome files
    temp_dir = "temp_bam_filter_chunks"
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    os.makedirs(temp_dir)

    print("\n--- Phase 1: Filtering chromosomes in parallel... ---")
    
    # Identify all chromosomes (references) available in the BAM file
    with pysam.AlignmentFile(args.input_bam, "rb") as infile:
        chromosomes = infile.references

    # Distribute the work across multiple CPU processes
    tasks = [(args.input_bam, temp_dir, chrom, target_barcodes) for chrom in chromosomes]
    temp_bam_files = []
    with ProcessPoolExecutor(max_workers=args.processes) as executor:
        results = executor.map(filter_chunk, tasks)
        for result in results:
            if result:
                temp_bam_files.append(result)

    print(f"Parallel filtering completed in {time.time() - start_time:.2f} seconds.")

    if not temp_bam_files:
        print("\nNo reads found for the given barcodes.")
        shutil.rmtree(temp_dir)
        return

    # --- PHASE 2: HIGH-SPEED MERGE ---
    # Instead of Python, we call the optimized 'samtools merge' via system subprocess
    print("\n--- Phase 2: Merging temporary BAM files using multiple threads... ---")
    merge_start_time = time.time()
    threads_for_merge = str(args.processes)
    merge_command = ["samtools", "merge", "-f", "-@", threads_for_merge, args.output_bam] + temp_bam_files
    
    try:
        subprocess.run(merge_command, check=True, capture_output=True, text=True)
        print(f"Merging completed in {time.time() - merge_start_time:.2f} seconds.")
    except subprocess.CalledProcessError as e:
        print(f"Error during samtools merge: {e.stderr}")
        shutil.rmtree(temp_dir)
        return
        
    print("\n--- Phase 3: Cleaning up temporary files... ---")
    shutil.rmtree(temp_dir)

    total_time = time.time() - start_time
    print(f"\nFinished. Total process took {total_time:.2f} seconds.")
    print("\nIMPORTANT: You must now index the new BAM file before using it in IGV:")
    print(f"samtools index {args.output_bam}")

if __name__ == "__main__":
    main()