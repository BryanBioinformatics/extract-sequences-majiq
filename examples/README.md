# Examples

This directory contains example files to help you get started with the junction annotation pipeline.

## Files

- `sample_input.csv` - A small sample of the input format with 3 junctions
- `sample_output_with_strand.csv` - Expected output after strand assignment
- `sample_output_with_sequences.csv` - Expected final output with sequences

## Sample Data

The sample input contains 3 junctions from different genes:

1. **TSPAN6** (chrX, negative strand) - A well-known gene
2. **FIRRM** (chr1, positive strand) - Another protein-coding gene  
3. **DPM1** (chr20, positive strand) - A housekeeping gene

## Expected Results

After running the pipeline on the sample data, you should get:

- **Strand assignment**: All 3 junctions should get strand information
- **Sequence extraction**: All 3 junctions should get 1002bp sequences
- **Biological accuracy**: Sequences should be in correct 5'→3' orientation

## Testing the Pipeline

To test with the sample data:

```bash
# Copy sample data to main directory
cp examples/sample_input.csv 17samples_imputed_min10coverage.csv

# Run the pipeline
./run_junction_pipeline.sh

# Check results
head -5 17samples_imputed_min10coverage_with_sequences.csv
```

## Validation

The sample data has been validated to ensure:
- ✅ Gene symbols match GTF file
- ✅ Junction coordinates are valid
- ✅ Strand information is correct
- ✅ Sequences are biologically accurate 