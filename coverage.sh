#!/bin/bash

#Download SRA toolkit
wget "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.2/sratoolkit.2.9.2-mac64.tar.gz"

#install it
tar xvzf sratoolkit.2.9.2-mac64.tar.gz

#use it to download the .fastq files
./sratoolkit.2.9.2-mac64/bin/fastq-dump --outdir . --split-files SRR961514

#download the genome fasta file from EBI 
wget -O sequence.fasta "https://www.ebi.ac.uk/ena/data/view/K03455&display=fasta"

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
conda install -c bioconda trimmomatic

#call trimmomatic
trimmomatic PE -threads 8 "$forward" "$reverse" -baseout "$base_out" AVGQUAL:"$quality"

