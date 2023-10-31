# This script checks the qualitiy of our fastq files and performs an alignment to the human cDNA transcriptome reference with Kallisto.
# To run this 'shell script' you will need to open your terminal and navigate to the directory where this script resides on your computer.
# This should be the same directory where you fastq files and reference fasta file are found.
# Change permissions on your computer so that you can run a shell script by typing: 'chmod 777 readMapping.sh' (without the quotes) at the terminal prompt 
# Then type './readMapping.sh' (without the quotes) at the prompt.  
# This will begin the process of running each line of code in the shell script.

# first use fastqc to check the quality of our fastq files:
###fastqc *.gz -t 4

# next, we want to build an index from our reference fasta file 
# I get my reference mammalian transcriptome files from here: https://useast.ensembl.org/info/data/ftp/index.html
###kallisto index -i ~/assignment/assignment_data/Solanum_lycopersicum.SL3.0.cdna.all.index ~/assignment/assignment_data/Solanum_lycopersicum.SL3.0.cdna.all.fa.gz

# now map reads to the indexed reference host transcriptome
# use as many 'threads' as your machine will allow in order to speed up the read mapping process.
# note that we're also including the '&>' at the end of each line
# this takes the information that would've been printed to our terminal, and outputs this in a log file that is saved in /data/course_data

# first the control subjects (CS)
kallisto quant -i ~/assignment/assignment_data/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/assignment_readmap/CS01 -t 4 --single -l 250 -s 30 ~/assignment/assignment_data/SRR10998560.fastq.gz &> ~/assignment/assignment_readmap/CS01.log
kallisto quant -i ~/assignment/assignment_data/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/assignment_readmap/CS02 -t 4 --single -l 250 -s 30 ~/assignment/assignment_data/SRR10998561.fastq.gz &> ~/assignment/assignment_readmap/CS02.log
kallisto quant -i ~/assignment/assignment_data/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/assignment_readmap/CS03 -t 4 --single -l 250 -s 30 ~/assignment/assignment_data/SRR10998562.fastq.gz &> ~/assignment/assignment_readmap/CS03.log
kallisto quant -i ~/assignment/assignment_data/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/assignment_readmap/CS04 -t 4 --single -l 250 -s 30 ~/assignment/assignment_data/SRR10998566.fastq.gz &> ~/assignment/assignment_readmap/CS04.log
kallisto quant -i ~/assignment/assignment_data/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/assignment_readmap/CS05 -t 4 --single -l 250 -s 30 ~/assignment/assignment_data/SRR10998569.fastq.gz &> ~/assignment/assignment_readmap/CS05.log
kallisto quant -i ~/assignment/assignment_data/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/assignment_readmap/CS06 -t 4 --single -l 250 -s 30 ~/assignment/assignment_data/SRR10998570.fastq.gz &> ~/assignment/assignment_readmap/CS06.log

# then the heat treatment (HT)
kallisto quant -i ~/assignment/assignment_data/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/assignment_readmap/HT07 -t 4 --single -l 250 -s 30 ~/assignment/assignment_data/SRR10998559.fastq.gz &> ~/assignment/assignment_readmap/HT07.log
kallisto quant -i ~/assignment/assignment_data/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/assignment_readmap/HT08 -t 4 --single -l 250 -s 30 ~/assignment/assignment_data/SRR10998563.fastq.gz &> ~/assignment/assignment_readmap/HT08.log
kallisto quant -i ~/assignment/assignment_data/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/assignment_readmap/HT09 -t 4 --single -l 250 -s 30 ~/assignment/assignment_data/SRR10998564.fastq.gz &> ~/assignment/assignment_readmap/HT09.log
kallisto quant -i ~/assignment/assignment_data/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/assignment_readmap/HT10 -t 4 --single -l 250 -s 30 ~/assignment/assignment_data/SRR10998565.fastq.gz &> ~/assignment/assignment_readmap/HT10.log
kallisto quant -i ~/assignment/assignment_data/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/assignment_readmap/HT11 -t 4 --single -l 250 -s 30 ~/assignment/assignment_data/SRR10998567.fastq.gz &> ~/assignment/assignment_readmap/HT11.log
kallisto quant -i ~/assignment/assignment_data/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/assignment_readmap/HT12 -t 4 --single -l 250 -s 30 ~/assignment/assignment_data/SRR10998568.fastq.gz &> ~/assignment/assignment_readmap/HT12.log

# summarize fastqc and kallisto mapping results in a single summary html using MultiQC
###multiqc -d . 

echo "Finished"

