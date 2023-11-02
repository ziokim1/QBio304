# build an index from Solanum lycopersicum reference fasta file 
kallisto index -i ~/assignment/read_map/Solanum_lycopersicum.SL3.0.cdna.all.index ~/assignment/fasta/Solanum_lycopersicum.SL3.0.cdna.all.fa.gz

# map reads to the indexed reference transcriptome

# first the control subjects (CS)
kallisto quant -i ~/assignment/read_map/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/read_map/CS01 -t 8 --single -l 250 -s 30 ~/assignment/fastq/trimmed/SRR10998560_trimmed.fastq.gz &> ~/assignment/read_map/log/CS01.log
kallisto quant -i ~/assignment/read_map/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/read_map/CS02 -t 8 --single -l 250 -s 30 ~/assignment/fastq/trimmed/SRR10998561_trimmed.fastq.gz &> ~/assignment/read_map/log/CS02.log
kallisto quant -i ~/assignment/read_map/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/read_map/CS03 -t 8 --single -l 250 -s 30 ~/assignment/fastq/trimmed/SRR10998562_trimmed.fastq.gz &> ~/assignment/read_map/log/CS03.log
kallisto quant -i ~/assignment/read_map/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/read_map/CS04 -t 8 --single -l 250 -s 30 ~/assignment/fastq/trimmed/SRR10998566_trimmed.fastq.gz &> ~/assignment/read_map/log/CS04.log
kallisto quant -i ~/assignment/read_map/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/read_map/CS05 -t 8 --single -l 250 -s 30 ~/assignment/fastq/trimmed/SRR10998569_trimmed.fastq.gz &> ~/assignment/read_map/log/CS05.log
kallisto quant -i ~/assignment/read_map/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/read_map/CS06 -t 8 --single -l 250 -s 30 ~/assignment/fastq/trimmed/SRR10998570_trimmed.fastq.gz &> ~/assignment/read_map/log/CS06.log

# then the heat treatment (HT)
kallisto quant -i ~/assignment/read_map/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/read_map/HT07 -t 8 --single -l 250 -s 30 ~/assignment/fastq/trimmed/SRR10998559_trimmed.fastq.gz &> ~/assignment/read_map/log/HT07.log
kallisto quant -i ~/assignment/read_map/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/read_map/HT08 -t 8 --single -l 250 -s 30 ~/assignment/fastq/trimmed/SRR10998563_trimmed.fastq.gz &> ~/assignment/read_map/log/HT08.log
kallisto quant -i ~/assignment/read_map/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/read_map/HT09 -t 8 --single -l 250 -s 30 ~/assignment/fastq/trimmed/SRR10998564_trimmed.fastq.gz &> ~/assignment/read_map/log/HT09.log
kallisto quant -i ~/assignment/read_map/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/read_map/HT10 -t 8 --single -l 250 -s 30 ~/assignment/fastq/trimmed/SRR10998565_trimmed.fastq.gz &> ~/assignment/read_map/log/HT10.log
kallisto quant -i ~/assignment/read_map/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/read_map/HT11 -t 8 --single -l 250 -s 30 ~/assignment/fastq/trimmed/SRR10998567_trimmed.fastq.gz &> ~/assignment/read_map/log/HT11.log
kallisto quant -i ~/assignment/read_map/Solanum_lycopersicum.SL3.0.cdna.all.index -o ~/assignment/read_map/HT12 -t 8 --single -l 250 -s 30 ~/assignment/fastq/trimmed/SRR10998568_trimmed.fastq.gz &> ~/assignment/read_map/log/HT12.log

echo "Finished"

# proceed to R scripts
