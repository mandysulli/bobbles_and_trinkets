#!/bin/bash
#$ -o mira-oxide.$JOB_ID.out
#$ -e mira-oxide.$JOB_ID.err
#$ -N mira-oxide_singularity
#$ -pe smp 4
#$ -l h_rt=24:00:00
#$ -l h_vmem=200G
#$ -q highmem.q
#$ -cwd
#$ -V

module load singularity

WORKDIR=/scicomp/groups/OID/NCIRD/ID-OD/VSDB/BIA/FLU_SC2_SEQUENCING/sra-2026-01-15

cd "$WORKDIR" || exit 1

singularity pull docker:cdcgov/mira-oxide:v1.3.1

cp /scicomp/groups/OID/NCIRD/ID-OD/VSDB/BIA/FLU_SC2_SEQUENCING/sra-2026-01-15/mira-output/aggregate_outputs/dais-ribosome/DAIS_ribosome.seq mira-output/manual_mira_oxide_run

echo "Running MIRA-Oxide prepare-mira-reports command"
echo "Working directory: $PWD"

singularity exec \
    --bind "$WORKDIR":"$WORKDIR":rw \
    --bind /scicomp/groups/OID/NCIRD/ID-OD/VSDB/BIA/MIRA-NF:/scicomp/groups/OID/NCIRD/ID-OD/VSDB/BIA/MIRA-NF:ro \
    mira-oxide_v1.3.1.sif \
    mira-oxide prepare-mira-reports \
    -s "$WORKDIR/samplesheet.csv" \
    -i "$WORKDIR/mira-output" \
    -o "$WORKDIR/mira-output/manual_mira_oxide_run" \
    -q /scicomp/groups/OID/NCIRD/ID-OD/VSDB/BIA/MIRA-NF/bin/irma_config/qc_pass_fail_settings.yaml \
    -p illumina \
    -w /scicomp/groups/OID/NCIRD/ID-OD/VSDB/BIA/MIRA-NF \
    -r sra-2026-01-15 \
    -v flu \
    -f
