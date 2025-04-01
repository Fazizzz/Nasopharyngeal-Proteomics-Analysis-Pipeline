#!/bin/bash

# Usage function - prints to the terminal and not the logs
usage() {
    echo "Usage: $0 -i input_directory -r regions_bed [-t threads]"
    echo "Options:"
    echo "  -i: Path to input directory containing BAM files (e.g., *_Aligned.sortedByCoord.out.bam)"
    echo "  -r: Path to BED file with regions of interest"
    echo "  -t: Number of threads to use (default: 8)"
    exit 1
}

# Default number of threads
THREADS=8
# Default column to pick up for human GTF for Read Counts. Change to whaterver is appropriate for your file
COLUMN=10

# Parse command-line arguments
while getopts ":i:r:t:" opt; do
    case $opt in
        i) INPUT_DIR="$OPTARG" ;;
        r) REGIONS_BED="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
    esac
done

# Check if all required arguments are provided
if [ -z "$INPUT_DIR" ] || [ -z "$REGIONS_BED" ]; then
    usage
fi

# Move into the input directory
cd "$INPUT_DIR" || { echo "Error: Could not access input directory $INPUT_DIR."; exit 1; }

# Define and create the logs directory
LOG_DIR="${INPUT_DIR}/logs"
mkdir -p "$LOG_DIR" || { echo "Error: Could not create log directory."; exit 1; }

# Create Primary-alignments and Read-counts directories inside input directory
mkdir -p "Primary-alignments" "Read-counts" || { echo "Error: Could not create required directories."; exit 1; }

# Redirect stdout and stderr to a general log file
exec > >(tee -a "$LOG_DIR/Counts_general_log.out") 2> >(tee -a "$LOG_DIR/Counts_error.log" >&2)

# Check if required tools are installed
if ! command -v bedtools &> /dev/null; then
    echo "Bedtools could not be found. Please install it and try again." >&2
    exit 1
fi

# Step 1: Index BAM files
echo "Indexing BAM files..."
ls *_Aligned.sortedByCoord.out.bam | xargs -I{} bash -c '
    input={}; output=${input%_Aligned.sortedByCoord.out.bam};
    if [ ! -f "${output}.bai" ]; then
        samtools index "$input" && echo "Indexed $input successfully." || echo "Error indexing $input" >> "${LOG_DIR}/indexing_errors.log";
    else
        echo "$input already indexed. Skipping..."
    fi
'

# Step 2: Extract primary alignments using samtools
echo "Extracting primary alignments..."
ls *_Aligned.sortedByCoord.out.bam | xargs -I{} bash -c '
    input={}; output=${input%_Aligned.sortedByCoord.out.bam}_primary.bam;
    if [ ! -f "Primary-alignments/${output}" ]; then
        samtools view -@ '"$THREADS"' -b -f 3 -F 2816 "$input" > "Primary-alignments/${output}" && echo "Extracted primary alignments for $input successfully." || echo "Error extracting primary alignments for $input" >> "${LOG_DIR}/primary_alignment_errors.log";
    else
        echo "Primary alignment for $input already exists. Skipping..."
    fi
'

# Step 3: Index primary BAM files
echo "Indexing primary BAM files..."
ls Primary-alignments/*.bam | xargs -I{} bash -c '
    input={}; output=${input%.bam};
    if [ ! -f "${output}.bai" ]; then
        samtools index "$input" && echo "Indexed $input successfully." || echo "Error indexing $input" >> "${LOG_DIR}/primary_indexing_errors.log";
    else
        echo "$input already indexed. Skipping..."
    fi
'

# Step 4: Get read counts using bedtools
echo "Calculating read counts for targets"
cd Primary-alignments || { echo "Error: Could not access Primary-alignments directory."; exit 1; }
ls *.bam | xargs -I{} bash -c '
    input={}; output=${input%_primary.bam}_cov.all;
    bedtools coverage -F 0.1 -b "$input" -a '"$REGIONS_BED"' > "../Read-counts/${output}" && echo "Calculated read counts for $input successfully." || echo "Error calculating read counts for $input" >> "${LOG_DIR}/read_counts_errors.log";
'

# Step 5: Collate results into a single TSV file
echo "Collating results..."
cd ../Read-counts || { echo "Error: Could not access Read-counts directory."; exit 1; }

for file in *.all; do
    sample_name="${file%_cov.all}"
    echo "$sample_name" > "${sample_name}.txt"
    cut -f $COLUMN "$file" >> "${sample_name}.txt"
done

echo -e "Chromosome\tStart\tStop" > rows
cat $(ls *all | shuf -n 1) | cut -f 1,4,5 >> rows # change columns here to 1,2,3 for a standard bed. this is for a GTF file

paste -d "\t" rows *.txt > "Read_counts.tsv"

# Copy the Read_counts.tsv file back to the input directory
cp "Read_counts.tsv" "$INPUT_DIR/Read_counts.tsv"

# Final completion message
echo "Post-alignment processing and quantification complete. Results saved in ${INPUT_DIR}/Read_counts.tsv"

echo "Script completed."

