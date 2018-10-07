#!/bin/bash

usage="Usage: [-f <fastq file of forward or 1 read>] [-r <fastq file of reverse or 2 read>] [-x <reference fasta file>] [-q <mean quality score threshold of reads to keep for mapping to reference, options are integers between 0 and 41, defaults to 0>] [-b <this is the base file name trimmomatic will use for the trimming output files>]"

while getopts ":f:r:x:q:b:" opt; do
	case "$opt" in
		f) forward="$OPTARG"
		echo "forward read file is $OPTARG" ;;
		r) reverse="$OPTARG"
		echo "reverse read file is $OPTARG" ;;
		x) reference="$OPTARG"
		echo "reference file is $OPTARG" ;;
		q) quality="$OPTARG"
		echo "quality threshold is $OPTARG" ;;
		b) base_out="$OPTARG"
		echo "trimmomatic output base file name is $OPTARG" ;;
		:) echo $usage
		exit 1
		;;
		\?) echo $usage
		exit 1
		;;
	esac
done

if [[ ! $quality =~ [0-9]{1,2} ]] || [[ $quality -gt 41 ]]; then quality=0
echo "quality value must be an integer between 0 and 41, using 0 default value"
fi

if [[ $OPTIND -lt 11 ]]; then
echo "You have not used all required arguments. ${usage}"
fi

if [[ $OPTIND -lt 11 ]]; then
exit
fi	


#call trimmomatic
trimmomatic PE -threads 8 "$forward" "$reverse" -baseout "$base_out" AVGQUAL:"$quality"

#make BWA index of reference
bwa index "$reference"

#align reads
bwa mem "$reference" "${base_out}_1P" "${base_out}_2P" > "${base_out}_mem.sam"

#sam to bam
samtools view -b "${base_out}_mem.sam" > "${base_out}_mem.bam"

#sort bam
samtools sort "${base_out}_mem.bam" > "${base_out}_mem_sorted.bam"

#index bam
samtools index "${base_out}_mem_sorted.bam"

#mpileup to get the coverage information per base
samtools mpileup "${base_out}_mem_sorted.bam" > "${base_out}_mem_sorted_pileup.txt"

#make faidx index
samtools faidx "$reference"

#make windows
bedtools makewindows -g "${reference}.fai" -w 25 -s 1 > "${reference}_windows.bed"

#convert bam alignment to bed
bedtools bamtobed -i "${base_out}_mem_sorted.bam" > "${base_out}_mem_sorted.bed"

#run bedtools coverage
bedtools coverage -a "${reference}_windows.bed" -b "${base_out}_mem_sorted.bed" -mean > "${base_out}_mem_sorted_bedcoverage.txt"

#run bedtools nuc
bedtools nuc -fi "$reference" -bed "${reference}_windows.bed" > "${reference}_windows_nuc.txt"

#print mpileup columns I want
awk 'BEGIN {OFS="\t"} {print $1, $2, $4}' "${base_out}_mem_sorted_pileup.txt" > "${base_out}_mem_sorted_pileup_coverage_only.txt"

#make header line for the mpileup output
echo $'chromosome\tcoordinate\tcoverage' > "${base_out}_mem_sorted_pileup_coverage_only_header.txt"

#add header to the mpileup
cat "${base_out}_mem_sorted_pileup_coverage_only_header.txt" "${base_out}_mem_sorted_pileup_coverage_only.txt" > "${base_out}_mem_sorted_pileup_coverage_only_for_plotting.txt"

cat "${base_out}_mem_sorted_pileup_coverage_only_for_plotting.txt" > "${quality}_coverage.tsv"

#Make headers for the bedcoverage file, use same conventions as the nuc file
echo $'#1_usercol\t2_usercol\t3_usercol\tcoverage' > "${base_out}_mem_sorted_bedcoverage_header.txt"

#add the header to the bedcoverage file
cat "${base_out}_mem_sorted_bedcoverage_header.txt" "${base_out}_mem_sorted_bedcoverage.txt" > "${base_out}_mem_sorted_bedcoverage_with_header.txt"

python python_report.py 
mv Report.pdf "${quality}_Report.pdf"
#clean up the folder
#make a dir
mkdir "${quality}""_""${base_out}"

#remove raw reads
rm "${base_out}_1P" 
rm "${base_out}_2P" 
rm "${base_out}_1U" 
rm "${base_out}_2U" 

#remove sam file
rm *.sam 

#move files to the new dir
mv "${reference}.bwt" ./"${quality}""_""${base_out}"
mv "${reference}.pac" ./"${quality}""_""${base_out}"
mv "${reference}.ann" ./"${quality}""_""${base_out}"
mv "${reference}.amb" ./"${quality}""_""${base_out}"
mv "${reference}.sa" ./"${quality}""_""${base_out}"
mv "${reference}_windows_nuc.txt" ./"${quality}""_""${base_out}"
mv *.bam* ./"${quality}""_""${base_out}"
mv *mem_sorted* ./"${quality}""_""${base_out}"
mv *.fai ./"${quality}""_""${base_out}"
mv *_windows.bed ./"${quality}""_""${base_out}"
mv *_Report.pdf ./"${quality}""_""${base_out}"
mv *_coverage.tsv ./"${quality}""_""${base_out}"


