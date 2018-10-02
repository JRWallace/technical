#!/bin/bash

#Download SRA toolkit
#wget "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.2/sratoolkit.2.9.2-mac64.tar.gz"

#install it
#tar xvzf sratoolkit.2.9.2-mac64.tar.gz

#use it to download the .fastq files
#./sratoolkit.2.9.2-mac64/bin/fastq-dump --outdir . --split-files SRR961514

#download the genome fasta file from EBI 
#wget -O sequence.fasta "https://www.ebi.ac.uk/ena/data/view/K03455&display=fasta"

#set the options cases -- still need to figure out how to include default behaviors and error messages for entries that are outside the expected inputs

while getopts "f:r:x:q:b:" opt; do
	case "$opt" in
		f) forward="$OPTARG";;
		r) reverse="$OPTARG";;
		x) reference="$OPTARG";;
		q) quality="$OPTARG";;
		b) base_out="$OPTARG";;	
	esac
done

#install trimmomatic
#conda install -c bioconda trimmomatic

#call trimmomatic
#trimmomatic PE -threads 8 "$forward" "$reverse" -baseout "$base_out" AVGQUAL:"$quality"

#install bwa, samtools, bedtools
#conda install -c bioconda bwa
#conda install -c bioconda samtools
#conda install -c bioconda bedtools

#make BWA index of reference
#bwa index "$reference"

#align reads
#bwa mem "$reference" "${base_out}_1P" "${base_out}_2P" > "${base_out}_mem.sam"

#sam to bam
#samtools view -b "${base_out}_mem.sam" > "${base_out}_mem.bam"

#sort bam
#samtools sort "${base_out}_mem.bam" > "${base_out}_mem_sorted.bam"

#index bam
#samtools index "${base_out}_mem_sorted.bam"

#mpileup to get the coverage information per base
#samtools mpileup "${base_out}_mem_sorted.bam" > "${base_out}_mem_sorted_pileup.txt"

#make faidx index
#samtools faidx "$reference"

#make windows
bedtools makewindows -g "${reference}.fai" -w 25 -s 1 > "${reference}_windows.bed"

#convert bam alignment to bed
bedtools bamtobed -i "${base_out}_mem_sorted.bam" > "${base_out}_mem_sorted.bed"

#run bedtools coverage
bedtools coverage -a "${reference}_windows.bed" -b "${base_out}_mem_sorted.bed" -mean > "${base_out}_mem_sorted_bedcoverage.txt"

#run bedtools nuc
bedtools nuc -fi "$reference" -bed "${reference}_windows.bed" > "${reference}_windows_nuc.txt" 







