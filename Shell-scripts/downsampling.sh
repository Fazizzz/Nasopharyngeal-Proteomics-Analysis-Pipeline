#!/bin/bash

# Log setup: Direct all logs to the current directory instead of $HOME
LOG_DIR="${LOG_DIR:-"./Downsampling_logs"}"  # Changed log directory to the current directory
mkdir -p "$LOG_DIR"  # Create logs directory if it doesn't exist

# General log file will contain overall progress
exec > >(tee -a "$LOG_DIR/Downsampling_general_log.out") 2> >(tee -a "$LOG_DIR/Downsampling_error.log" >&2)  # Logs stdout and stderr

# Usage function to display help for the script
# This function will print directly to the terminal (stdout) and won't be piped to logs
usage() {
    echo "Usage: $0 -d input_directory -n number_of_reads -o output_directory [-t threads]"
    echo "-d: Input directory containing FASTQ files"
    echo "-n: Number of reads to sample"
    echo "-o: Output directory for downsampled FASTQ files"
    echo "-t: Optional, number of threads (default: 30)"
    exit 1
}

# Default settings
THREADS=30  # Default number of threads set to 30

# Parse command-line arguments
while getopts ":d:n:o:t:" opt; do  # Removed the -p option for seqtk path
    case $opt in
        d) FASTQ_DIR="$OPTARG" ;;  # Directory containing FASTQ files
        n) READS_COUNT="$OPTARG" ;;  # Number of reads to sample
        o) DOWNSAMPLED_FASTQ_DIR="$OPTARG" ;;  # Output directory
        t) THREADS="$OPTARG" ;;  # Optional: number of threads (default: 30)
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;  # Invalid option handler
    esac
done

# Check if all required arguments are provided
if [ -z "$FASTQ_DIR" ] || [ -z "$READS_COUNT" ] || [ -z "$DOWNSAMPLED_FASTQ_DIR" ]; then
    usage
fi

# Create the output directory if it doesn't exist
mkdir -p "$DOWNSAMPLED_FASTQ_DIR"

# Check if seqtk is installed
if ! command -v seqtk &> /dev/null; then  # Removed the -p seqtk_path handling, just call seqtk directly
    echo "Seqtk could not be found. Please install it and try again."
    exit 1
fi

echo "Starting downsampling of reads with Seqtk..."

# Loop through all filtered _R1.fastq.gz or _R1.fastq files and find their pairs (_R2.fastq.gz or _R2.fastq)
for FASTQ_FILE in "$FASTQ_DIR"/*_R1_filtered.fastq*; do
    BASENAME=$(basename "$FASTQ_FILE" _R1_filtered.fastq.gz)  # Get the base filename without extension
    [[ $FASTQ_FILE == *.fastq ]] && BASENAME=$(basename "$FASTQ_FILE" _R1_filtered.fastq)  # Handle uncompressed fastq files

    # Define the paired files for R1 and R2
    FASTQ_FILE_R1="${FASTQ_DIR}/${BASENAME}_R1_filtered.fastq.gz"
    FASTQ_FILE_R2="${FASTQ_DIR}/${BASENAME}_R2_filtered.fastq.gz"

    # Handle uncompressed .fastq files if compressed ones are not present
    if [ ! -f "$FASTQ_FILE_R1" ]; then
        FASTQ_FILE_R1="${FASTQ_DIR}/${BASENAME}_R1_filtered.fastq"
    fi
    if [ ! -f "$FASTQ_FILE_R2" ]; then
        FASTQ_FILE_R2="${FASTQ_DIR}/${BASENAME}_R2_filtered.fastq"
    fi

    if [[ -f "$FASTQ_FILE_R1" && -f "$FASTQ_FILE_R2" ]]; then
        echo "Processing $BASENAME..."

        # Skip if output files already exist
        if [ -f "$DOWNSAMPLED_FASTQ_DIR/${BASENAME}_R1_downsampled.fastq.gz" ] && [ -f "$DOWNSAMPLED_FASTQ_DIR/${BASENAME}_R2_downsampled.fastq.gz" ]; then
            echo "Downsampled files for $BASENAME already exist. Skipping..."
            continue
        fi

        # Log file for individual FASTQ file errors
        SAMPLE_ERROR_LOG="$LOG_DIR/${BASENAME}_error.log"

        # Downsample R1
        if seqtk sample -s100 "$FASTQ_FILE_R1" "$READS_COUNT" | gzip > "$DOWNSAMPLED_FASTQ_DIR/${BASENAME}_R1_downsampled.fastq.gz"; then
            echo "$BASENAME R1 downsampled successfully."
        else
            echo "Error: Downsampling R1 failed for $BASENAME" | tee -a "$SAMPLE_ERROR_LOG"
            continue
        fi

        # Downsample R2
        if seqtk sample -s100 "$FASTQ_FILE_R2" "$READS_COUNT" | gzip > "$DOWNSAMPLED_FASTQ_DIR/${BASENAME}_R2_downsampled.fastq.gz"; then
            echo "$BASENAME R2 downsampled successfully."
        else
            echo "Error: Downsampling R2 failed for $BASENAME" | tee -a "$SAMPLE_ERROR_LOG"
            continue
        fi
    else
        echo "Error: Missing paired file for $BASENAME" | tee -a "$LOG_DIR/error.log"
    fi
done

echo "Downsampling completed."


