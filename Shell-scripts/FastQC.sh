#!/bin/bash

# Combined log (stderr and stdout) setup
# ----------------------------------------
# LOG_DIR: This variable sets the default directory for logs if not already defined.
# mkdir -p: This command creates the directory if it doesn't exist.
LOG_DIR="${LOG_DIR:-"./FastQC_logs"}"
mkdir -p "$LOG_DIR"

# GENERAL_LOG: This will store all combined output (general log + error log) for the entire script.
GENERAL_LOG="$LOG_DIR/FastQC_General.log"

# Redirect both stdout and stderr to the general log file
# --------------------------------------------------------
# This redirection will ensure that all echo statements, outputs, and errors are logged into general_log.out
exec > >(tee -a "$GENERAL_LOG") 2>&1

# Trap to reset file descriptors on exit
# ---------------------------------------
# This trap ensures that original file descriptors are restored when the script exits.
trap 'exec 2>&4 1>&3' 0 1 2 3
exec 3>&1 4>&2

# Function to display script usage (help message)
# ------------------------------------------------
# usage(): This function will be called if the script is not provided with the required arguments.
usage() {
    echo "Usage: $0 -d input_directory -o output_directory [-t threads]" | tee -a "$GENERAL_LOG"
    exit 1
}

# Default number of threads
THREADS=8

# Parse command-line arguments
# ---------------------------------
# This block reads the input options and assigns the values to the respective variables.
while getopts ":d:o:t:" opt; do
    case $opt in
        d) TRIMMED_FASTQ_DIR="$OPTARG" ;;  # Input directory for trimmed FASTQ files
        o) TRIMMED_QC_DIR="$OPTARG" ;;     # Output directory for FastQC results
        t) THREADS="$OPTARG" ;;            # Number of threads
        \?) echo "Invalid option -$OPTARG" | tee -a "$GENERAL_LOG"  # Log invalid options
            usage ;;
    esac
done

# Check if all required arguments are provided
# ----------------------------------------------
# This block ensures that both input and output directories are provided, otherwise, it triggers the usage() function.
if [ -z "$TRIMMED_FASTQ_DIR" ] || [ -z "$TRIMMED_QC_DIR" ]; then
    usage
fi

# Create output directory if it doesn't exist
# ----------------------------------------------
# This block ensures that the output directory for FastQC results exists. If not, it creates the directory.
if [ ! -d "$TRIMMED_QC_DIR" ]; then
    echo "Output directory does not exist. Creating: $TRIMMED_QC_DIR" | tee -a "$GENERAL_LOG"
    mkdir -p "$TRIMMED_QC_DIR"
    if [[ $? -ne 0 ]]; then
        echo "Error: Could not create output directory $TRIMMED_QC_DIR" | tee -a "$GENERAL_LOG"
        exit 1
    fi
fi

# Check if FastQC is installed
# -------------------------------
# This checks if the "fastqc" command is available on the system. If it's not installed, the script will exit.
if ! command -v fastqc &> /dev/null; then
    echo "FastQC could not be found. Please install it and try again." | tee -a "$GENERAL_LOG"
    exit 1
fi

# Check if input directory contains FASTQ files
# ------------------------------------------------
# If the input directory doesn't contain any files with the ".fastq.gz" extension, the script exits with a message.
if [ -z "$(ls -A "$TRIMMED_FASTQ_DIR"/*.fastq.gz 2>/dev/null)" ]; then
    echo "No trimmed FASTQ files found in the input directory: $TRIMMED_FASTQ_DIR" | tee -a "$GENERAL_LOG"
    exit 1
fi

# Running FastQC on trimmed FASTQ files
# -----------------------------------------
# This is the main loop that processes each FASTQ file in the directory.
echo "Running FastQC on trimmed FASTQ files..." | tee -a "$GENERAL_LOG"

for TRIMMED_FASTQ_FILE in "$TRIMMED_FASTQ_DIR"/*.fastq.gz; do
    BASENAME=$(basename "$TRIMMED_FASTQ_FILE" .fastq.gz)

    # Individual log for each file
    # ------------------------------
    # Each sample will have its own log in the log directory to capture FastQC-specific output and errors.
    SAMPLE_LOG="$LOG_DIR/${BASENAME}_fastqc_trimmed.log"

    # Skip if output already exists
    # ------------------------------
    # If FastQC output (HTML) already exists for the current sample, the script skips reprocessing that sample.
    if [ -f "$TRIMMED_QC_DIR/${BASENAME}_fastqc.html" ]; then
        echo "FastQC output for $BASENAME already exists. Skipping..." | tee -a "$GENERAL_LOG"
        continue
    fi

    # Process each FASTQ file using FastQC
    # --------------------------------------
    # FastQC command is run, and its output is redirected to an individual sample log.
    echo "Processing $TRIMMED_FASTQ_FILE..." | tee -a "$GENERAL_LOG"
    fastqc -t $THREADS -o "$TRIMMED_QC_DIR" "$TRIMMED_FASTQ_FILE" > "$SAMPLE_LOG" 2>&1

    # Check if FastQC ran successfully
    # ---------------------------------
    # $? checks the exit status of the last command (FastQC). If FastQC fails, an error message is logged.
    if [[ $? -ne 0 ]]; then
        echo "Error: FastQC failed for $TRIMMED_FASTQ_FILE" | tee -a "$GENERAL_LOG"
    else
        echo "FastQC for $BASENAME completed." | tee -a "$GENERAL_LOG"
    fi
done

# Final completion message
# --------------------------
# A final message is logged when the quality control process completes.
echo "Quality control on trimmed data completed." | tee -a "$GENERAL_LOG"

