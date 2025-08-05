# Lough Feeagh sedaDNA analysis

These commands and pipelines are the bioinformatic steps which accompany the Tighe et al. 2025 article, showing how both the metabarcoding data and shotgun sequences were analysed.
All analyses were run on a Linux Server running Ubuntu 24.04.2 LTS, with 64 cores and 472 Gb memory.

## Analysing the metabarcoding data

### Firstly install cutadapt, R, RDP classifier (which requires java) and local BLAST.
```
sudo apt install cutadapt
R

wget https://sourceforge.net/projects/rdp-classifier/files/rdp-classifier/rdp_classifier_2.13.zip
unzip rdp_classifier_2.13

BLAST
```
### Download the relevent databases for RDP classifier
```
wget https://github.com/terrimporter/CO1Classifier/releases/download/RDP-COI-v5.1.0/RDP_COIv5.1.0.zip
wget https://github.com/terrimporter/12SvertebrateClassifier/releases/tag/v3.0.0-ref/12SvertebrateNA_v3.0.0_ref.zip
```
### Run cutadapt to demultiplex the samples
```
mkdir 1_demux

cutadapt -j 50 -e 0.02 --no-indels -g file:$forwardtag -G file:$reversetag -o 1_demux/{name1}-{name2}_R1_001.fastq.gz -p 1_demux/{name1}-{name2}_R2_001.fastq.gz $fwdfile $revfile --minimum-length 100 > fwd_orient.cutadapt.stat
```




## Analysing the shotgun sequence data

### Install panda to assemble paired-end reads.
```
run panda
```
### BLAST all assembled and unassembled raw reads.
```
./blastn -query rawreads.fasta -db /data01/nt_2024/nt -num_threads 50 -max_target_seqs 10 -outfmt "6 qseqid stitle pident length evalue" -out Raw_blast_results.txt
```
### Custom BASH scripts to remove all human, bacterial and fungal reads from the original database. It searches using awk in the first field (tab delimited) for the Genus name and then pulls out the second field, which is the query ID.
Example of first three commands of BASH script:
```
#!/bin/bash

awk  -F '\t' '$1 ~ /A.radiobacter/ {print $2}' 2nt_nucl_test6_myown_95_w44_assembled > bad_accesions.txt | echo "Pulling out A.radiobacter accessions.."
awk  -F '\t' '$1 ~ /Abiotrophia/ {print $2}' 2nt_nucl_test6_myown_95_w44_assembled >> bad_accesions.txt | echo "Pulling out Abiotrophia accessions.."
awk  -F '\t' '$1 ~ /Achromobacter/ {print $2}' 2nt_nucl_test6_myown_95_w44_assembled >> bad_accesions.txt | echo "Pulling out Achromobacter accessions.."
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




