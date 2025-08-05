# Lough Feeagh sedaDNA analysys

These commands and pipelines are the bioinformatic steps which accompany the Tighe et al. 2025 article, showing how both the metabarcoding data and shotgun sequences were analysed.
All analyses were run on XX Linux Server running 11 (?), with 60 cores and 1Tb memory.

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




## Analysing the metabarcoding data

### Install panda to assemble paired-end reads
```
run panda
```
### BLAST all assmebled and unassembled raw reads
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





