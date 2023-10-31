#!/bin/bash
# Set path to trimmomatic


# Loop through each sample
for sample in SRR10998570 SRR10998569 SRR10998566 SRR10998565 SRR10998564 SRR10998563 SRR10998562 SRR10998561 SRR10998560 SRR10998559 SRR10998568 SRR10998567
do
  # Set input and output file paths
  input1=~/assignment/assignment_data/${sample}.fastq.gz 
  output1=~/assignment/assignment_data_trimmed/${sample}_trimmed.fastq.gz 
  
  # Run trimmomatic
  java -jar ~/assignment/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 3 $input1 $output1 \
  ILLUMINACLIP:./Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 HEADCROP:10 MINLEN:65
done