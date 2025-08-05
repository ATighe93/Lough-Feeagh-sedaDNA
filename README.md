# Lough-Feagh-sedaDNA

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
