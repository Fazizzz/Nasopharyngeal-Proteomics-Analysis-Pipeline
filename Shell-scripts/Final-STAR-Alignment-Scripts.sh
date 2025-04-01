#!/bin/bash

# Define and create the logs directory
LOG_DIR="./logs"  
mkdir -p "$LOG_DIR" || { echo "Error: Could not create log directory."; exit 1; }

# Redirect stdout and stderr to the general log file, except for usage which prints to terminal
exec > >(tee -a "$LOG_DIR/Alignment_general_log.out") 2> >(tee -a "$LOG_DIR/Alignment_error.log" >&2)

# Usage function - prints to the terminal and not the logs
usage() {
    echo "Usage: $0 -i input_fastq -o output_dir [-f fasta_file | -g genome_index] [-t threads] [-s gtf_file]"
    echo "-i: Input FASTQ file or directory"
    echo "-o: Output directory"
    echo "-f: FASTA file to build genome index"
    echo "-g: STAR Genome index directory (if already built)"
    echo "-t: Optional, number of threads (default: 8)"
    echo "-s: Optional, GTF file for genome index generation"
    exit 1
}

# Default settings
THREADS=8
READ_LENGTH=99
FASTA_FILE=""
GENOME_INDEX=""
GTF_FILE=""

# Parse command-line arguments
while getopts ":r:o:f:g:t:s:" opt; do
    case $opt in
        r) INPUT_FASTQ="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        f) FASTA_FILE="$OPTARG" ;;
        g) GENOME_INDEX="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        s) GTF_FILE="$OPTARG" ;;  # Added GTF file option
        \?) echo "Invalid option: -$OPTARG" >&2; usage ;;
    esac
done

# Check if required options are provided
if [[ -z "$INPUT_FASTQ" || -z "$OUTPUT_DIR" ]]; then
    usage
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR" || { echo "Error: Could not create output directory $OUTPUT_DIR"; exit 1; }

# Generate genome index if FASTA is provided
if [[ -n "$FASTA_FILE" ]]; then
    echo "FASTA file provided. Generating genome index in $OUTPUT_DIR/genome_index..."
    mkdir -p "$OUTPUT_DIR/genome_index" || { echo "Error: Could not create genome index directory."; exit 1; }

    # Add GTF file option if provided
    if [[ -n "$GTF_FILE" ]]; then
        STAR --runThreadN "$THREADS" --runMode genomeGenerate --genomeDir "$OUTPUT_DIR/genome_index" \
             --genomeFastaFiles "$FASTA_FILE" --sjdbOverhang $READ_LENGTH --sjdbGTFfile "$GTF_FILE" > "$LOG_DIR/genome_index.log" 2>> "$LOG_DIR/error.log"
    else
        STAR --runThreadN "$THREADS" --runMode genomeGenerate --genomeDir "$OUTPUT_DIR/genome_index" \
             --genomeFastaFiles "$FASTA_FILE"  > "$LOG_DIR/genome_index.log" 2>> "$LOG_DIR/error.log"
    fi

    if [[ $? -ne 0 ]]; then
        echo "Error: Genome index generation failed." >> "$LOG_DIR/error.log"
        exit 1
    fi
    GENOME_INDEX="$OUTPUT_DIR/genome_index"
    echo "Genome index successfully generated in $OUTPUT_DIR/genome_index"
elif [[ -n "$GENOME_INDEX" ]]; then
    echo "Using provided genome index: $GENOME_INDEX"
else
    echo "Error: No FASTA file or genome index provided."
    usage
fi

# Check if STAR genome index was generated successfully
if [[ ! -d "$GENOME_INDEX" ]]; then
    echo "Error: Genome index directory not found."
    exit 1
fi

# Run STAR alignment
echo "Running STAR alignment for $INPUT_FASTQ..."
for FASTQ_FILE in "$INPUT_FASTQ"/*_R1_downsampled.fastq.gz; do
    SAMPLE_NAME=$(basename "$FASTQ_FILE" "_R1_downsampled.fastq.gz")
    FASTQ_R2="${INPUT_FASTQ}/${SAMPLE_NAME}_R2_downsampled.fastq.gz"
    
    # Skip if no paired FASTQ_R2 is found
    if [[ ! -f "$FASTQ_R2" ]]; then
        echo "Error: Paired FASTQ file for $FASTQ_FILE not found. Skipping..." | tee -a "$LOG_DIR/error.log"
        continue
    fi

    # Log file for individual sample alignment
    SAMPLE_LOG="$LOG_DIR/${SAMPLE_NAME}_alignment.log"
    
    STAR --runThreadN "$THREADS" --genomeDir "$GENOME_INDEX" --readFilesIn "$FASTQ_FILE" "$FASTQ_R2" \
         --readFilesCommand zcat --outFileNamePrefix "$OUTPUT_DIR/${SAMPLE_NAME}_" \
         --outSAMtype BAM SortedByCoordinate > "$SAMPLE_LOG" 2>> "$LOG_DIR/error.log"
         
    # Check if STAR ran successfully
    if [[ $? -ne 0 ]]; then
        echo "Error: STAR alignment failed for $SAMPLE_NAME. Check $SAMPLE_LOG for details." | tee -a "$LOG_DIR/error.log"
    else
        echo "STAR alignment for $SAMPLE_NAME completed successfully."
    fi
done

# Check for STAR alignment summaries
echo "Generating Alignment Summary..."
echo "Sample" > "$LOG_DIR/samplename.txt"
echo "Number of input reads" > "$LOG_DIR/inputReads.txt"
echo "Average input read length" > "$LOG_DIR/readlength.txt"
echo "% Uniquely mapped reads" > "$LOG_DIR/unique.txt"
echo "Average mapped length" > "$LOG_DIR/mapLen.txt"
echo "% of reads mapped to multiple loci" > "$LOG_DIR/loci.txt"
echo "% of reads mapped to too many loci" > "$LOG_DIR/2loci.txt"
echo "% reads unmapped by mismatches" > "$LOG_DIR/missMatch.txt"
echo "% reads unmapped by length" > "$LOG_DIR/lenny.txt"
echo "% reads unmapped by other" > "$LOG_DIR/other.txt"

for file in "$OUTPUT_DIR"/*Log.final.out; do
    SAMPLE=$(basename "$file" Log.final.out)  # use file basename
    echo "$SAMPLE" >> "$LOG_DIR/samplename.txt"
    grep "Number of input reads" "$file" | cut -f 2 >> "$LOG_DIR/inputReads.txt"
    grep "Average input read length" "$file" | cut -f 2 >> "$LOG_DIR/readlength.txt"
    grep "Uniquely mapped reads %" "$file" | cut -f 2 >> "$LOG_DIR/unique.txt"
    grep "Average mapped length" "$file" | cut -f 2 >> "$LOG_DIR/mapLen.txt"
    grep "% of reads mapped to multiple loci" "$file" | cut -f 2 >> "$LOG_DIR/loci.txt"
    grep "% of reads mapped to too many loci" "$file" | cut -f 2 >> "$LOG_DIR/2loci.txt"
    grep "% of reads unmapped: too many mismatches" "$file" | cut -f 2 >> "$LOG_DIR/missMatch.txt"
    grep "% of reads unmapped: too short" "$file" | cut -f 2 >> "$LOG_DIR/lenny.txt"
    grep "% of reads unmapped: other" "$file" | cut -f 2 >> "$LOG_DIR/other.txt"
done

paste "$LOG_DIR/samplename.txt" "$LOG_DIR/inputReads.txt" "$LOG_DIR/readlength.txt" "$LOG_DIR/unique.txt" \
      "$LOG_DIR/mapLen.txt" "$LOG_DIR/loci.txt" "$LOG_DIR/2loci.txt" "$LOG_DIR/missMatch.txt" "$LOG_DIR/lenny.txt" \
      "$LOG_DIR/other.txt" > "$LOG_DIR/combined-starLogFinal.txt"

echo "Alignment Summary generated at $LOG_DIR/combined-starLogFinal.txt"

# Cleanup temporary log files
rm "$LOG_DIR/samplename.txt" "$LOG_DIR/inputReads.txt" "$LOG_DIR/readlength.txt" "$LOG_DIR/unique.txt" \
   "$LOG_DIR/mapLen.txt" "$LOG_DIR/loci.txt" "$LOG_DIR/2loci.txt" "$LOG_DIR/missMatch.txt" "$LOG_DIR/lenny.txt" \
   "$LOG_DIR/other.txt"

echo "Script completed."

