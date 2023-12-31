#!/bin/bash

# copy adapter to the current folder
cp ~/assignment/Trimmomatic-0.39/adapters/TruSeq3-SE.fa .

# Loop through each sample
for sample in SRR10998570 SRR10998569 SRR10998566 SRR10998565 SRR10998564 SRR10998563 SRR10998562 SRR10998561 SRR10998560 SRR10998559 SRR10998568 SRR10998567
do
  # Set input and output file paths
  input1=~/assignment/fastq/raw/${sample}.fastq.gz 
  output1=~/assignment/fastq/trimmed/${sample}_trimmed.fastq.gz 
  
  # Run trimmomatic
  java -jar ~/assignment/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 8 $input1 $output1 \
  ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 HEADCROP:10 MINLEN:65
done

# run quality check again on the trimmed reads
fastqc ~/assignment/fastq/trimmed/*.gz -t 8 -o ~/assignment/qc/trimmed
multiqc -d ~/assignment/qc/trimmed -o ~/assignment/qc/trimmed

echo "Finished Step 2"

# if the trimmed reads seem to pass the quality check, proceed to read mapping.
