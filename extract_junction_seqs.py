#!/usr/bin/env python3
"""
extract_junction_seqs.py

Given a CSV of junctions with columns "junction" and "strand", this script will:
1) Drop all junctions with missing strand information (".")
2) Extract 250bp of flanking context + the site itself around both donor and acceptor
3) Use strand information to extract sequences in the correct orientation
4) Concatenate donor and acceptor sequences (1002bp total)
5) Add the result as a new column "sequence" after the "strand" column

Usage:
    pip install pyfaidx pandas
    ./extract_junction_seqs.py \
        --input   psi_ddliu.csv \
        --fasta   Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
        --output  with_sequences.csv
"""

import argparse
import pandas as pd
from pyfaidx import Fasta

def normalize_chrom(fasta, chrom):
    """
    Normalize contig names to match the FASTA:
      1) if there's a '|' (e.g. 'chrX|Y'), take the first part ('chrX')
      2) try raw chrom
      3) strip leading "chr"
      4) add leading "chr"
    """
    if "|" in chrom:
        chrom = chrom.split("|")[0]

    if chrom in fasta:
        return chrom
    no_chr = chrom[3:] if chrom.startswith("chr") else chrom
    if no_chr in fasta:
        return no_chr
    with_chr = "chr" + chrom if not chrom.startswith("chr") else chrom
    if with_chr in fasta:
        return with_chr

    raise KeyError(f"Chromosome '{chrom}' not found (tried {chrom}, {no_chr}, {with_chr}).")

def get_context_sequence(fasta, chrom, pos, flank=250):
    """
    Fetch flank bp upstream + site + flank bp downstream.
    Pads with 'N' if we run off the end of the chromosome.
    """
    contig = normalize_chrom(fasta, chrom)

    # pyfaidx uses 0-based indexing, end-exclusive
    start_idx = pos - 1 - flank
    end_idx   = pos - 1 + flank + 1

    # pad left
    if start_idx < 0:
        pad_left = "N" * (-start_idx)
        start_idx = 0
    else:
        pad_left = ""

    seq = fasta[contig][start_idx:end_idx].seq

    # pad right if needed
    desired_len = 2 * flank + 1
    if len(seq) < desired_len:
        seq = seq + "N" * (desired_len - len(seq))

    return pad_left + seq

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 
                  'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    return ''.join(complement[base] for base in reversed(seq))

def main():
    # File paths
    input_csv = '17samples_imputed_min10coverage_with_strand.csv'
    genome_fasta = 'references/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa'
    output_csv = '17samples_imputed_min10coverage_with_sequences.csv'

    # load
    print("Loading input CSV...")
    df = pd.read_csv(input_csv, dtype=str)
    print(f"Loaded {len(df)} total junctions")
    
    # Check for required columns
    if 'junction' not in df.columns:
        raise ValueError("Input CSV must have a 'junction' column")
    if 'strand' not in df.columns:
        raise ValueError("Input CSV must have a 'strand' column")
    
    # Drop junctions with missing strand information
    print("Filtering junctions with missing strand information...")
    original_count = len(df)
    df = df[df['strand'].notna() & (df['strand'] != '') & (df['strand'] != '.')]
    dropped_count = original_count - len(df)
    print(f"Dropped {dropped_count} junctions with missing strand information")
    print(f"Remaining junctions: {len(df)}")
    
    if len(df) == 0:
        print("No junctions with strand information found. Exiting.")
        return
    
    # Load FASTA
    print("Loading reference genome...")
    fasta = Fasta(genome_fasta, rebuild=False)

    # parse junctions: chrom:donor-acceptor
    print("Parsing junction coordinates...")
    parts = df['junction'].str.split('[:\-]', expand=True)
    df['chrom']    = parts[0]
    df['donor']    = parts[1].astype(int)
    df['acceptor'] = parts[2].astype(int)

    # extract sequences with strand-aware orientation
    print("Extracting sequences with strand-aware orientation...")
    def build_seq(row):
        # Get donor and acceptor sequences
        d_seq = get_context_sequence(fasta, row.chrom, row.donor, flank=250)
        a_seq = get_context_sequence(fasta, row.chrom, row.acceptor, flank=250)
        
        # For negative strand, we need to reverse complement the sequences
        # and swap donor/acceptor order to maintain biological orientation
        if row.strand == '-':
            d_seq = reverse_complement(d_seq)
            a_seq = reverse_complement(a_seq)
            # For negative strand, concatenate acceptor first, then donor
            return (a_seq + d_seq).upper()
        else:
            # For positive strand, concatenate donor first, then acceptor
            return (d_seq + a_seq).upper()

    df['sequence'] = df.apply(build_seq, axis=1)

    # clean up: drop temp cols, move 'sequence' right after 'strand'
    df = df.drop(columns=['chrom','donor','acceptor'])
    cols = list(df.columns)
    cols.remove('sequence')
    idx = cols.index('strand') + 1
    new_order = cols[:idx] + ['sequence'] + cols[idx:]
    df = df[new_order]

    # write out
    print(f"Writing {len(df)} junctions with sequences to {output_csv}")
    df.to_csv(output_csv, index=False)
    print(f"Successfully wrote {len(df)} junctions with uppercase sequences to {output_csv}")
    
    # Print strand distribution
    strand_counts = df['strand'].value_counts()
    print(f"\nStrand distribution in output:")
    for strand, count in strand_counts.items():
        print(f"  {strand}: {count} junctions")

if __name__ == "__main__":
    main()