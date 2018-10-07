### Pipeline overview

##### This is a pipeline to map a set of paired-end MiSeq reads to a reference genome. It requires the forward and reverse fastq files and a reference fasta file. It also requires arguments to be passed when the script is called in order to set an average read quality threshold and an output file name, both for trimmomatic. The reads, reference, coverage.sh, and python_report.py should all be in the same directory. The following example would run the pipeline with SRR961514_1.fastq as the forward read, SRR961514_2.fastq as the reverse read, sequence.fasta as the reference file, a quality threshold of 30, and SRR961514 as a base file name for trimmomatic output files. 

```
bash coverage.sh -f SRR961514_1.fastq -r SRR961514_2.fastq -x sequence.fasta -q 30 -b SRR961514
```
##### The quality threshold defaults to 0 and will accept values between 0 and 41. All arguments are required for the pipeline to proceed.

### Pipeline details

1. The first part of the script uses trimmomatic to trim the reads. In this particular case, it was requested for the reads to be filtered based on a reads average quality score. To my knowledge, re-trimming the reads is the best way to meet this requirement. BWA and SAMtools will allow filtering based on *mapping* quality or the quality of a *base*, but not the average quality over the length of the entire read. There was no other trimming or read quality control performed.

2. The next step is indexing the reference and mapping the reads using BWA. I used BWA mem, as it is designed to align 70-1Mbp queries (as per the BWA documentation). This particular set of MiSeq reads was from a 2 x 250 sequencing run, so BWA mem seemed appropriate.

3. After running the mapping, I used SAMtools to convert the SAM to BAM, sort and index it, and used mpileup to get the read coverage. You could also probably use the SAMtools depth option, but from what I understand the mpileup output is more conservative when considering what counts as a mapped read. The mpileup file did not print out with file headers, so I added those separately.

4. I also used SAMtools to index the reference sequence. I used bedtools to set a sliding window of 25 nt and a 1nt step size (this is the windows_.bed file). I converted the alignment BAM file to a BED file and ran bedtools nuc and coverage to get information about the nucleotide content and coverage in each of the sliding windows. The coverage file did not print out with file headers, so I added those separately.

5. The next step is to run the python_report.py script. After it finishes, I change the file names to reflect the quality filtering that is reflected in the report. This seemed more straight forward than passing a bash variable to a python script.

6. The last step of the script does a clean up of the directory. It makes a directory (named for the quality threshold and trimmomatic output file name), and puts the BAM file into that folder. I remove all of the other outputs because they are quite large, but the BAM file seems like something that is useful and should be kept.

### Python_report details
1. The first part of this script imports the files using glob2, since the exact file names will change depending on the quality thresholds used.

2. The next part of the script merges the bedtools nuc and coverage files and renames the file headers to remove special characters.

3. I run the correlation analyses (kendall, spearman, pearson) on the whole dataframe and write a function to bin the output of the correlation. If the absolute value of the correlation is less than 0.3, it is considered weak, if it is greater than 0.7 it is considered strong, and if it is between 0.3 and 0.7 it is considered moderate.

4. The next few lines deal with formatting the dataframes. Since we are particularly interested in how nucleotide composition correlates with sequencing depth, I kept the column containing the correlation score for this variable only.

5. The last part of this script generates the report.pdf and moves the outputs into a new folder named for the quality cutoff and the trimmomatic output file name. The first statement generates the plot showing coverage over the length of the reference, and the next 3 add tables for the correlation analyses. The coverage_corr column is the correlation score from the particular test that was run, and the third column references the strength of the correlation. The variables in the variable column are as follows:

  - w_start: the starting position of the sliding window
  - w_end: the ending position of the sliding window
  - pct_at: percent AT in window
  - pct_gc: percent GC in window
  - num_A: number of A in window
  - num_C: number of C in window
  - num_G: number of G in window
  - num_T: number of T in window
  - num_N: number of N in window
  - num_OTH: number of bases other than A, T, G, C, or N in window
  - seq_len: length of the sliding window
  - coverage: average depth of sequencing in the window
