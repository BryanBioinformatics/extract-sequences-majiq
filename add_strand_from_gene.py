#!/usr/bin/env python3

import pandas as pd
import gzip
from collections import defaultdict

def parse_gtf_for_gene_strands(gtf_file):
    """Parse GTF file to extract gene symbol -> strand mapping"""
    gene_strands = {}
    
    print(f"Parsing GTF file: {gtf_file}")
    
    # Handle gzipped or plain text files
    if gtf_file.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open
    
    with opener(gtf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            feature_type = parts[2]
            if feature_type != 'gene':
                continue
            
            # Parse attributes
            attributes = parts[8]
            gene_id = None
            gene_name = None
            
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr.startswith('gene_id'):
                    gene_id = attr.split('"')[1] if '"' in attr else attr.split(' ')[1]
                elif attr.startswith('gene_name'):
                    gene_name = attr.split('"')[1] if '"' in attr else attr.split(' ')[1]
            
            strand = parts[6]
            
            # Store both gene_id and gene_name mappings
            if gene_id:
                gene_strands[gene_id] = strand
            if gene_name:
                gene_strands[gene_name] = strand
    
    print(f"Found {len(gene_strands)} gene strand mappings")
    return gene_strands

def extract_gene_from_lsv(lsv_id):
    """Extract gene symbol from LSV ID format: GENE:s:coordinates:coordinates"""
    if ':' not in lsv_id:
        return None
    
    # Split by ':' and take the first part as gene symbol
    parts = lsv_id.split(':')
    if len(parts) >= 1:
        return parts[0]
    return None

def main():
    # File paths
    input_csv = '17samples_imputed_min10coverage.csv'
    gtf_file = 'references/Homo_sapiens.GRCh38.113.gtf'
    output_csv = '17samples_imputed_min10coverage_with_strand.csv'
    
    print("=== Adding Strand Information from Gene Symbols ===")
    
    # Parse GTF file to get gene strand information
    gene_strands = parse_gtf_for_gene_strands(gtf_file)
    
    # Read the input CSV
    print(f"Reading input CSV: {input_csv}")
    df = pd.read_csv(input_csv)
    
    print(f"Input CSV shape: {df.shape}")
    print(f"Sample LSV IDs:")
    print(df.iloc[:5, 0].tolist())
    
    # Extract gene symbols from LSV IDs
    print("Extracting gene symbols from LSV IDs...")
    df['gene_symbol'] = df.iloc[:, 0].apply(extract_gene_from_lsv)
    
    # Count unique gene symbols
    unique_genes = df['gene_symbol'].dropna().unique()
    print(f"Found {len(unique_genes)} unique gene symbols")
    print(f"Sample gene symbols: {unique_genes[:10].tolist()}")
    
    # Look up strand information
    print("Looking up strand information...")
    df['strand'] = df['gene_symbol'].map(gene_strands)
    
    # Count strand assignments
    strand_counts = df['strand'].value_counts()
    print(f"Strand distribution:")
    print(strand_counts)
    
    # Count missing strands
    missing_strands = df['strand'].isna().sum()
    print(f"Junctions with missing strand information: {missing_strands}")
    
    # Show some examples of missing strands
    if missing_strands > 0:
        print("\nSample junctions with missing strand information:")
        missing_examples = df[df['strand'].isna()].head(10)
        for _, row in missing_examples.iterrows():
            print(f"  LSV: {row.iloc[0]}")
            print(f"  Gene: {row['gene_symbol']}")
            print(f"  Junction: {row.iloc[1]}")
            print()
    
    # Save the output
    print(f"Saving output to: {output_csv}")
    df.to_csv(output_csv, index=False)
    
    print("=== Summary ===")
    print(f"Total junctions: {len(df)}")
    print(f"Junctions with strand information: {len(df) - missing_strands}")
    print(f"Junctions missing strand information: {missing_strands}")
    print(f"Success rate: {((len(df) - missing_strands) / len(df) * 100):.1f}%")

if __name__ == "__main__":
    main() 