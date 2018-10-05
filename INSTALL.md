# Installation information

## download SRA toolkit

```
wget "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.2/sratoolkit.2.9.2-mac64.tar.gz"
```

## uncompress it
```
tar xvzf sratoolkit.2.9.2-mac64.tar.gz
```
## use it to download the .fastq files
```
./sratoolkit.2.9.2-mac64/bin/fastq-dump --outdir . --split-files SRR961514
```
## download the genome fasta file from EBI
```
wget -O sequence.fasta "https://www.ebi.ac.uk/ena/data/view/K03455&display=fasta"
```

## install bioinformatics packages trimmomatic, bwa, samtools, bedtools
```
conda install -c bioconda trimmomatic
conda install -c bioconda bwa
conda install -c bioconda samtools
conda install -c bioconda bedtools
```
## install other python packages
```
conda install -c conda-forge matplotlib
conda install -c conda-forge pandas
conda install -c conda-forge glob2     
```
