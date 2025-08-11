#!/bin/bash

# Script to rename and move FASTQ files from bamtofastq output
# Usage: rename_fastqs.sh <source_dir> <sample_name> <target_dir>

set -euo pipefail

# Check if correct number of arguments provided
if [ $# -ne 3 ]; then
    echo "Usage: $0 <sample_name> <source_dir> <target_dir>"
    echo "Example: $0 SRR30256650 /path/to/source /path/to/target"
    exit 1
fi

SAMPLE_NAME="$1"
SOURCE_DIR="$2"
TARGET_DIR="$3"

# Check if source directory exists
if [ ! -d "$SOURCE_DIR" ]; then
    echo "Error: Source directory '$SOURCE_DIR' does not exist"
    exit 1
fi

# Create target directory if it doesn't exist
mkdir -p "$TARGET_DIR"

echo "Renaming and moving FASTQ files..."
echo "Source: $SOURCE_DIR"
echo "Sample: $SAMPLE_NAME"
echo "Target: $TARGET_DIR"

# Ensure target directory exists
mkdir -p "$TARGET_DIR"

# Initialize sample counter
sample_counter=1

# Process each subdirectory in the source directory
for subdir in "$SOURCE_DIR"/*; do
    if [ -d "$subdir" ]; then
        subdir_name=$(basename "$subdir")
        echo "Processing directory: $subdir_name (S${sample_counter})"
        
        # Process all FASTQ files in the subdirectory
        for fastq_file in "$subdir"/*.fastq.gz; do
            if [ -f "$fastq_file" ]; then
                # Extract the original filename
                original_filename=$(basename "$fastq_file")
                
                # Extract the suffix after "bamtofastq_"
                # This will capture everything after "bamtofastq_" (e.g., "L001_R1_001.fastq.gz")
                if [[ "$original_filename" =~ bamtofastq_S1_(.+) ]]; then
                    original_suffix="${BASH_REMATCH[1]}"
                    
                    # Create new filename: ${sample_name}_S${sample_number}_${original_suffix}
                    new_filename="${SAMPLE_NAME}_S${sample_counter}_${original_suffix}"
                    
                    # Move and rename the file
                    echo "  $original_filename -> $new_filename"
                    mv "$fastq_file" "$TARGET_DIR/$new_filename"
                else
                    echo "  Warning: Could not parse filename pattern for $original_filename"
                fi
            fi
        done
        
        # Increment sample counter for next subdirectory
        ((sample_counter++))
    fi
done

echo "Completed processing. Files moved to: $TARGET_DIR"

# Count total files processed
total_files=$(find "$TARGET_DIR" -name "*.fastq.gz" | wc -l)
echo "Total FASTQ files processed: $total_files"