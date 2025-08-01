#!/bin/bash

# Junction Annotation and Sequence Extraction Pipeline
# Simplified approach using gene symbols from LSV IDs

set -e  # Exit on any error

echo "=========================================="
echo "Junction Annotation and Sequence Pipeline"
echo "Simplified Gene-Based Approach"
echo "=========================================="

# Configuration
INPUT_CSV="17samples_imputed_min10coverage.csv"
GTF_FILE="references/Homo_sapiens.GRCh38.113.gtf"
GENOME_FASTA="references/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
OUTPUT_WITH_STRAND="17samples_imputed_min10coverage_with_strand.csv"
OUTPUT_WITH_SEQUENCES="17samples_imputed_min10coverage_with_sequences.csv"

# Check if input file exists
if [ ! -f "$INPUT_CSV" ]; then
    echo "Error: Input file $INPUT_CSV not found!"
    exit 1
fi

# Check if reference files exist
if [ ! -f "$GTF_FILE" ]; then
    echo "Error: GTF file $GTF_FILE not found!"
    exit 1
fi

if [ ! -f "$GENOME_FASTA" ]; then
    echo "Error: Genome FASTA file $GENOME_FASTA not found!"
    exit 1
fi

echo ""
echo "Step 1: Adding strand information from gene symbols..."
echo "Input: $INPUT_CSV"
echo "GTF: $GTF_FILE"
echo "Output: $OUTPUT_WITH_STRAND"
echo ""

python3 add_strand_from_gene.py

if [ ! -f "$OUTPUT_WITH_STRAND" ]; then
    echo "Error: Strand assignment failed - output file not created!"
    exit 1
fi

echo ""
echo "Step 2: Extracting junction sequences..."
echo "Input: $OUTPUT_WITH_STRAND"
echo "Genome: $GENOME_FASTA"
echo "Output: $OUTPUT_WITH_SEQUENCES"
echo ""

python3 extract_junction_seqs.py

if [ ! -f "$OUTPUT_WITH_SEQUENCES" ]; then
    echo "Error: Sequence extraction failed - output file not created!"
    exit 1
fi

echo ""
echo "=========================================="
echo "Pipeline completed successfully!"
echo "=========================================="
echo ""
echo "Final output: $OUTPUT_WITH_SEQUENCES"
echo ""

# Show final statistics
echo "Final statistics:"
echo "=================="
TOTAL_JUNCTIONS=$(wc -l < "$OUTPUT_WITH_SEQUENCES")
TOTAL_JUNCTIONS=$((TOTAL_JUNCTIONS - 1))  # Subtract header line
echo "Total junctions with sequences: $TOTAL_JUNCTIONS"

# Count strand distribution
echo ""
echo "Strand distribution:"
POSITIVE_STRAND=$(awk -F',' 'NR>1 && $(NF-1)=="+" {count++} END {print count+0}' "$OUTPUT_WITH_SEQUENCES")
NEGATIVE_STRAND=$(awk -F',' 'NR>1 && $(NF-1)=="-" {count++} END {print count+0}' "$OUTPUT_WITH_SEQUENCES")
echo "  Positive strand (+): $POSITIVE_STRAND junctions"
echo "  Negative strand (-): $NEGATIVE_STRAND junctions"

# Check sequence length
echo ""
echo "Sequence information:"
SEQUENCE_LENGTH=$(awk -F',' 'NR==2 {print length($NF)}' "$OUTPUT_WITH_SEQUENCES")
echo "  Sequence length: ${SEQUENCE_LENGTH}bp"

# File size
FILE_SIZE=$(ls -lh "$OUTPUT_WITH_SEQUENCES" | awk '{print $5}')
echo "  File size: $FILE_SIZE"

echo ""
echo "Pipeline completed successfully!"
echo "Output file: $OUTPUT_WITH_SEQUENCES" 