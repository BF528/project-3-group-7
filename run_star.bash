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
#$ -N star_p3

# Combine output and error files into a single file
#$ -j y

# Specify the output file name
#$ -o star.qlog
module load gcc
module load star/2.6.0c
for i in $(ls *.fastq.gz | sed 's/_[12].fastq.gz//' | uniq); do STAR --genomeDir /projectnb2/bf528/project_3/reference/rn4_STAR/ \
--readFilesIn ${i}_1.fastq.gz ${i}_2.fastq.gz \
--runThreadN 16  \
--outFileNamePrefix aligned/$i. \
--outSAMtype BAM SortedByCoordinate \
--readFilesCommand zcat ; done
