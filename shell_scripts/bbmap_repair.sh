#!/bin/bash
source /etc/profile

#$ -N fix_files
#$ -q all.q
#$ -cwd
#$ -m abe
#$ -o fix_files.out
#$ -e fix_files.err
#$ -l h_rt=4:00:00
#$ -l h_vmem=32G

module load BBMap/38.90

for i in SRR28525021 SRR28525065 SRR28567179 SRR30131960 SRR30131961 SRR30131962 SRR30131963 SRR30131964 SRR30131965 SRR30131966 SRR30131967 SRR30131968 SRR30131969 SRR30131970 SRR30131971 SRR30131972 SRR28958948 SRR28958949 SRR28958950 SRR28958951 SRR28958952 SRR28958953 SRR28958954 SRR28958955 SRR28958956 SRR28958957 SRR28958958 SRR28958959 SRR28958960 SRR28958961; do
    repair.sh in=/scicomp/home-pure/xpa3/mira_cli_test_data/large_flu_illumina/unfixed_reads/$i\_R1.fastq in2=/scicomp/home-pure/xpa3/mira_cli_test_data/large_flu_illumina/unfixed_reads/$i\_R2.fastq out=/scicomp/home-pure/xpa3/mira_cli_test_data/large_flu_illumina/fastqs/$i\_R1.fastq out2=/scicomp/home-pure/xpa3/mira_cli_test_data/large_flu_illumina/fastqs/$i\_R2.fastq outs=/scicomp/home-pure/xpa3/mira_cli_test_data/large_flu_illumina/singletons/$i\_S.fastq repair=t
done
