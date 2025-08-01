#!/bin/bash

# Junction Annotation Pipeline Setup Script

echo "=========================================="
echo "Junction Annotation Pipeline Setup"
echo "=========================================="

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ Conda is not installed. Please install Miniconda or Anaconda first."
    echo "   Download from: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo "✅ Conda found"

# Create conda environment
echo ""
echo "Creating conda environment..."
conda env create -f environment.yml

if [ $? -eq 0 ]; then
    echo "✅ Conda environment created successfully!"
else
    echo "❌ Failed to create conda environment"
    exit 1
fi

# Make scripts executable
echo ""
echo "Making scripts executable..."
chmod +x run_junction_pipeline.sh

echo "✅ Scripts made executable"

# Create references directory
echo ""
echo "Creating references directory..."
mkdir -p references

echo "✅ References directory created"

echo ""
echo "=========================================="
echo "Setup completed successfully!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "1. Activate the environment: conda activate junction-annotation"
echo "2. Add your input files:"
echo "   - Put your CSV file in the current directory"
echo "   - Put GTF and FASTA files in references/ directory"
echo "3. Run the pipeline: ./run_junction_pipeline.sh"
echo ""
echo "For more information, see README.md" 