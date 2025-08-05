# dada2 pipeline Rscript #

library("dada2")

directory = getwd()
setwd(file.path(directory,"fwd_dada"))			# Where you want the output files to save
path <- file.path(directory,"1_demux/true_hits")	 		# Where the fastq.gz demuxd files are located

fnFs <- sort(list.files(path, pattern="_1.fastq.gz", full.names = TRUE)) 	   # list forward files 
fnRs <- sort(list.files(path,pattern="_2.fastq.gz", full.names = TRUE)) 	   # list reverse files
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Quality good both ways (checked with FastQC), will truncate F at , reverse at  

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(218,218), minLen=70, maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)

exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]

sink("filterAndTrim.txt")
head(out,20)
sink()

sink("leanrErrors.txt")
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
sink()

make.monotone.decreasing <- function(v) sapply(seq_along(v), function(i) max(v[i:length(v)])) # enforce monotonicity in error estimates for Novaseq data
errF.md <- t(apply(getErrors(errF), 1, make.monotone.decreasing))
dimnames(errF.md) <- dimnames(errF$err_out)
errR.md <- t(apply(getErrors(errR), 1, make.monotone.decreasing))
dimnames(errR.md) <- dimnames(errR$err_out)

save.image("dada2.RData")

errF$err_out <- errF.md
pdf("errFmd.pdf")
plotErrors(errF, nominalQ=TRUE)
dev.off()

errR$err_out <- errR.md
pdf("errRmd.pdf")
plotErrors(errR, nominalQ=TRUE)
dev.off()

sink("dadaFs.txt")
dadaFs <- dada(filtFs, err=errF, multithread=TRUE, pool=TRUE)
sink()
sink("dadaRs.txt")
dadaRs <- dada(filtRs, err=errR, multithread=TRUE, pool=TRUE)
sink()

save.image("dada2.RData")

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, minOverlap = 40, verbose=TRUE)
# Inspect the merger data.frame from the first sample
sink("head_mergers.txt")
head(mergers[[1]])
sink()

save.image("dada2.RData")

seqtab <- makeSequenceTable(mergers)		# Make ASV table

write.csv(table(nchar(getSequences(seqtab))), "seq_lengths.csv")
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 303:323] # lengths set as 303-323 bp
write.csv(table(nchar(getSequences(seqtab2))), "seq_lengths_filt.csv")
write.csv(seqtab2, "seqtab2.csv")

save.image("dada2.RData")

sink("chim_rem.txt")
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab2)
sink()

save.image("dada2.RData")

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, "track_fwd_orient.csv")
save.image("dada2.RData")

uniquesToFasta(seqtab.nochim, "ASVs_1.fasta", ids=colnames(seqtab.nochim))
save.image("dada2.RData")
saveRDS(seqtab.nochim, "seqtab_nochim.rds")
