#!/bin/bash
# Change permissions on your computer so that you can run a shell script by typing: 'chmod 777 step1.sh' (without the quotes) at the terminal prompt 
# Then type './step1.sh' (without the quotes) at the prompt.  

# make assignment directory with the necessary subdirectories
mkdir -p ~/assignment/{fastq/{raw,trimmed},fasta,qc/{raw,trimmed},read_map/log}

# download fastq files
wget -nc -P ~/assignment/fastq/raw ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR109/063/SRR10998563/SRR10998563.fastq.gz
wget -nc -P ~/assignment/fastq/raw ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR109/064/SRR10998564/SRR10998564.fastq.gz
wget -nc -P ~/assignment/fastq/raw ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR109/062/SRR10998562/SRR10998562.fastq.gz
wget -nc -P ~/assignment/fastq/raw ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR109/059/SRR10998559/SRR10998559.fastq.gz
wget -nc -P ~/assignment/fastq/raw ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR109/061/SRR10998561/SRR10998561.fastq.gz
wget -nc -P ~/assignment/fastq/raw ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR109/065/SRR10998565/SRR10998565.fastq.gz
wget -nc -P ~/assignment/fastq/raw ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR109/069/SRR10998569/SRR10998569.fastq.gz
wget -nc -P ~/assignment/fastq/raw ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR109/067/SRR10998567/SRR10998567.fastq.gz
wget -nc -P ~/assignment/fastq/raw ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR109/068/SRR10998568/SRR10998568.fastq.gz
wget -nc -P ~/assignment/fastq/raw ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR109/060/SRR10998560/SRR10998560.fastq.gz
wget -nc -P ~/assignment/fastq/raw ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR109/066/SRR10998566/SRR10998566.fastq.gz
wget -nc -P ~/assignment/fastq/raw ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR109/070/SRR10998570/SRR10998570.fastq.gz

# download md5sums for fastq files and make sure reads are downloaded correctly
wget -nc -P ~/assignment/fastq/raw raw.githubusercontent.com/ziokim1/QBio304/main/md5sums.txt

# download reference cDNA for Solanum lycopersicum
wget -nc -P ~/assignment/fasta ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/solanum_lycopersicum/cdna/Solanum_lycopersicum.SL3.0.cdna.all.fa.gz

# download Trimmomatic
wget -nc -P ~/assignment http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
unzip ~/assignment/Trimmomatic-0.39.zip -d ~/assignment

# download study design
wget -nc -P ~/assignment/read_map raw.githubusercontent.com/ziokim1/QBio304/main/studydesign_assignment.txt

# run fastqc on fastq files
fastqc ~/assignment/fastq/raw/*.gz -t 8 -o ~/assignment/qc/raw
multiqc -d ~/assignment/qc/raw -o ~/assignment/qc/raw

# looking at the multiqc file and individual fastqc files, make the judgment whether to trim reads or not.

echo "Finished Step 1"
