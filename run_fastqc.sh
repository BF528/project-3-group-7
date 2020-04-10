#!/bin/bash -l

#$ -P bf528
# Specify hard time limit for the job. 
#   The job will be aborted if it runs longer than this time.
#   The default time is 12 hours
#$ -pe mpi_16_tasks_per_node 16
#$ -l h_rt=48:00:00

# Send an email when the job finishes or if it is aborted (by default no email is sent).
#$ -m ea

# Give job a name
#$ -N fastqc_p3

# Combine output and error files into a single file
#$ -j y

# Specify the output file name
#$ -o fastqc.qlog

echo "Job started: $(date +%F)"
fastqc *.fastq.gz -o qc_data

echo "Job finished: $(date +%F)"