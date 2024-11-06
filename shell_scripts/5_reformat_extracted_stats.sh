#!/bin/bash
source /etc/profile

#$ -N reformat_extracts
#$ -q all.q
#$ -cwd
#$ -M xpa3@cdc.gov
#$ -m abe
#$ -o reformat_extracts.out
#$ -e reformat_extracts.err
#$ -l h_rt=00:30:00
#$ -l h_vmem=1G

cd /scicomp/scratch/xpa3/bam_check

sed -e 's/'\ SN\ average\ length\:\ '/ /' extract_stats_d.out >extract_stats_d_2.txt
sed -e 's/'\ SN\ average\ length\:\ '/ /' extract_stats_g.out >extract_stats_g_2.txt

sed -e 's/'\ SN\ average\ quality\:\ '/ /' extract_stats_d.out >extract_stats_d.csv
sed -e 's/'\ SN\ average\ quality\:\ '/ /' extract_stats_g.out >extract_stats_g.csv
