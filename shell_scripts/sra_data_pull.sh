#!/bin/bash
source /etc/profile

#$ -N SRA_data_pull
#$ -q all.q
#$ -cwd
#$ -m abe
#$ -o SRA_data_pull.out
#$ -e SRA_data_pull.err
#$ -l h_rt=00:30:00
#$ -l h_vmem=20G

module load sratoolkit/2.10.9

for i in SRR28525021 SRR28525065 SRR28567179 SRR28567180 SRR30131960 SRR30131961 SRR30131962 SRR30131963 SRR30131964 SRR30131965 SRR30131966 SRR30131967 SRR30131968 SRR30131969 SRR30131970 SRR30131971 SRR30131972 SRR28958948 SRR28958949 SRR28958950 SRR28958951 SRR28958952 SRR28958953 SRR28958954 SRR28958955 SRR28958956 SRR28958957 SRR28958958 SRR28958959 SRR28958960 SRR28958961; do
    #prefetch -O /scicomp/home-pure/xpa3/mira_cli_test_data/large_flu_illumina/sra_files $i
    cd /scicomp/home-pure/xpa3/mira_cli_test_data/large_flu_illumina/sra_files/$i
    fasterq-dump --split-files $i\.sra -O /scicomp/home-pure/xpa3/mira_cli_test_data/large_flu_illumina/fastqs_redo/
done
