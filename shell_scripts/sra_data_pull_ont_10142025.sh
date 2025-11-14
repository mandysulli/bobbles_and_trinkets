#!/bin/bash
#$ -o sra_ont.$JOB_ID.out
#$ -e sra_ont.$JOB_ID.err
#$ -N SRA_pull_ont
#$ -pe smp 4
#$ -l h_rt=12:00:00
#$ -l h_vmem=12G
#$ -q flu.q
#$ -cwd
#$ -V

module load sra-toolkit/3.2.1

# Define the working directory
WORKDIR="/scicomp/groups-pure/OID/NCIRD/ID-OD/VSDB/BIA/MHS/large_flu_test_data/large_flu_wgs_ont"
SRA_DIR="$WORKDIR/sra_ont_pull"
FASTQ_DIR="$WORKDIR/fastq_pass"

# Ensure the directories exist
mkdir -p "$SRA_DIR"
mkdir -p "$FASTQ_DIR"

# Define the list of samples to download
SAMPLES=(
    "SRR35746125"
    "SRR32055894"
    "SRR32055466"
    "SRR23852495"
    "SRR33435879"
    "SRR31649752"
    "ERR15108850"
)

# Loop through each sample
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing sample: $SAMPLE"

    # Prefetch the SRA file
    prefetch --output-directory "$SRA_DIR" "$SAMPLE"

    # Check if the prefetch was successful
    if [ $? -ne 0 ]; then
        echo "Error: Failed to prefetch $SAMPLE"
        continue
    fi

    # Construct the path to the .sra file
    SRA_FILE="$SRA_DIR/$SAMPLE/$SAMPLE.sra"

    # Verify the .sra file exists
    if [ ! -f "$SRA_FILE" ]; then
        echo "Error: SRA file not found for $SAMPLE at $SRA_FILE"
        continue
    fi

    # Run fasterq-dump on the downloaded SRA file without splitting
    fasterq-dump "$SRA_FILE" --split-3 -O "$FASTQ_DIR"

    # Check if fasterq-dump was successful
    if [ $? -ne 0 ]; then
        echo "Error: Failed to convert $SAMPLE to FASTQ"
        continue
    fi

    echo "Successfully processed $SAMPLE"
done

# Gzip all FASTQ files in the fastqs folder
echo "Compressing FASTQ files in $FASTQ_DIR..."
find "$FASTQ_DIR" -type f -name "*.fastq" -exec gzip {} \;

# Check if gzip was successful
if [ $? -ne 0 ]; then
    echo "Error: Failed to compress FASTQ files"
else
    echo "Successfully compressed all FASTQ files in $FASTQ_DIR"
fi
