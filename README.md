# Lough Feeagh sedaDNA analysis

These commands and pipelines are the bioinformatic steps which accompany the Tighe et al. 2025 article, showing how both the metabarcoding data and shotgun sequences were analysed.
All analyses were run on a Linux Server running Ubuntu 24.04 LTS, with 64 cores and 472 Gb memory.

## Analysing the metabarcoding data

### Firstly install cutadapt, R, RDP classifier (which requires java) and local BLAST.
```
pip3 install --user cutadapt

sudo apt install r-base

wget https://downloads.sourceforge.net/project/rdp-classifier/rdp-classifier/2.13/rdp_classifier_2.13.zip
unzip rdp_classifier_2.13

sudo apt install ncbi-blast+
```
### Download the relevent databases for RDP classifier and local BLAST
```
wget https://github.com/terrimporter/CO1Classifier/releases/download/RDP-COI-v5.1.0/RDP_COIv5.1.0.zip
wget https://github.com/terrimporter/12SvertebrateClassifier/releases/tag/v3.0.0-ref/12SvertebrateNA_v3.0.0_ref.zip

wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt.*.tar.gz
```
### Run cutadapt to demultiplex the samples
```
mkdir 1_demux

cutadapt -j 50 -e 0.02 --no-indels -g file:forward_primer_tags.txt -G file:reverse_primer_tags.txt -o 1_demux/{name1}-{name2}_R1_001.fastq.gz -p 1_demux/{name1}-{name2}_R2_001.fastq.gz fwdsequences.fq revseequences.fq --minimum-length 100 > fwd_orient.cutadapt.stat
```
### Merge and denoise the forward and reverse sequence files, and then generate actual sequence variants (ASVs), using DADA2 pipeline 
Rename your demultiplexed files to actual sample names.
```
mkdir true_hits
mv	r1-f1_R1_001.fastq	true_hits/Sample1_1.fastq
mv	r1-f1_R2_001.fastq	true_hits/Sample1_2.fastq
...
```
Run dada scripts specific to the primers and gene used.
```
mkdir fwd_dada
Rscript dada_1.R
```
### Phylogenetically place each ASV using RDP classifer, and merge results to sample information using custom R scripts.
```
java -Xmx20g -jar /path_to_classifier/RDP/rdp_classifier_213/dist/classifier.jar classify -t /path_to_reference_database/DBs/specific_directory/rRNAClassifier.properties -o fwd_dada/RDP_classified_fwd.txt fwd_dada/ASVs_1.fasta

Rscript Fwd_abundance_table.R
```
### Double check each phylogenetic placement using local BLAST.

First pull out the ASV sequences from the output table from Fwd_abundance_table.R, and put them in a fasta format suitable for BLAST. Command for comma separated csv, 'NR!=1' skips the first line (the header), gsub swaps any " for nothing, print ">"$1 prints the first column with an ">" added in front, and "\n" means move the second column to a new line.
```
awk -F "," 'NR!=1 {gsub(/"/, ""); print ">"$1, "\n"$2}' Fwd_abundance_table.csv > ASV_fwd_for_blast.fasta
```
BLASTed against the entire GenBank database locally.
```
blastn -query ASV_fwd_for_blast.fasta -db /pathway_to_local_database/nt/2nt -num_threads 50 -max_target_seqs 10 -outfmt "6 qseqid stitle pident length evalue" -out ASV_fwd_blast_results.txt
```
Print top hit of each result only.
```
sort -k1,1 -u ASV_fwd_blast_results.txt > ASV_fwd_blast_top_hit.csv
```
Merge the blast results with the RDP results using custom R script.
```
Rscript Fwd_combine_blast.R
```
## Analysing the shotgun sequence data

### Install and run PANDASEQ to assemble paired-end reads.
```
sudo apt install pandaseq
pandaseq -f forward.fastq -r reverse.fastq
```
### BLAST all assembled and unassembled raw reads.
```
./blastn -query rawreads.fasta -db /pathway_to_local_database/nt_2024/nt -num_threads 50 -max_target_seqs 10 -outfmt "6 qseqid stitle pident length evalue" -out Raw_blast_results.txt
```
### Custom BASH scripts to remove all human, bacterial and fungal reads from the original database. It searches using awk in the first field (tab delimited) for the Genus name and then pulls out the second field, which is the query ID.
Example of first three commands of BASH script:
```
#!/bin/bash

awk  -F '\t' '$1 ~ /A.radiobacter/ {print $2}' 2nt_nucl_test6_myown_95_w44_assembled > bad_accesions.txt | echo "Pulling out A.radiobacter accessions.."
awk  -F '\t' '$1 ~ /Abiotrophia/ {print $2}' 2nt_nucl_test6_myown_95_w44_assembled >> bad_accesions.txt | echo "Pulling out Abiotrophia accessions.."
awk  -F '\t' '$1 ~ /Achromobacter/ {print $2}' 2nt_nucl_test6_myown_95_w44_assembled >> bad_accesions.txt | echo "Pulling out Achromobacter accessions.."
...
```
Next used awk to make a file of unique accessions that blasted to human or microbial.
```
awk '{ a[$1]++ } END { for (b in a) { print b } }' bad_accesions.txt > Unique_bad_accesions.txt
```
Then used awk to print the accessions in field 2 of "Raw_blast_results.txt" which are not in the Unique bad accession list, and the output is piped into Good_accessions.txt.
```
awk -F '\t' 'NR == FNR { list[tolower($1)]=1; next } { if (! list[tolower($2)]) { print $2 } }' Unique_bad_accesions.txt 2nt_nucl_test6_myown_95_w44_assembled > Good_accessions.txt
```
Get the unique query accessions and put them in a new file.
```
awk '{ a[$1]++ } END { for (b in a) { print b } }' Good_accessions.txt > Unique_Good_accessions.txt
```
Pull the good sequences out of the original fastq file.
```
grep -A3 -Ff Unique_Good_accessions.txt Raw_reads.assembled.fastq > Reads_for_reblasting.fastq
```
The output from above has an extra line per sequence with --, so this is removed using reverse grep.
```
grep -v -e "--" Reads_for_reblasting.fastq > Reads_for_reblasting_nodash.fastq
```

### Preprocessing, using cutadapt to trim polyG tails, remove any reads of low quality (Q20) and any low complexity reads below 30 bp.
```
cutadapt --quality-cutoff 20 -a "G{10}" --minimum-length 30 -o trimmed_Reads_for_reblasting.fastq Reads_for_reblasting_nodash.fastq
```
### Running BLAST locally on the processed sequences, with the standard output (-oufmt 6) which can be read by MEGAN.
```
blastn -query Trimmed_reads_for_reblasting.fastq -db /data04/nt/2nt -num_threads 50 -max_target_seqs 10 -outfmt 6 -out Clean_blast_results.txt
```




