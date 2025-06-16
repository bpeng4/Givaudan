setwd("/work/benson/bpeng4/Givaudan/Givaudan_All_16S")
suppressMessages({
library("tidyverse")
library("dada2")
library("gridExtra")
library("devtools")
})
no_of_cores = 36
metadata = "Meta.csv"
metadata_df = read.csv(metadata)
metadata_df = metadata_df[,-1]
fnFs <- metadata_df$fq1
fnRs <- metadata_df$fq2

sample.names <- metadata_df$Sample_Name

filt_path <- paste0(getwd(), '/part1filtered') # don't change this 
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

FORWARD_TRUNC <- 250 # determine from quality plots
REVERSE_TRUNC <- 160 # determine from quality plots

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(FORWARD_TRUNC,REVERSE_TRUNC), 
                     trimLeft=c(20, 20), maxEE=c(2,2), 
                     multithread=no_of_cores,
                     matchIDs=TRUE, compress=TRUE, 
                     verbose=TRUE)

derepFs <- derepFastq(filtFs, n = 1e+06, verbose = TRUE)
derepRs <- derepFastq(filtRs, n = 1e+06, verbose = TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

errF <- learnErrors(filtFs, verbose=TRUE, multithread=no_of_cores)
errR <- learnErrors(filtRs, verbose=TRUE, multithread=no_of_cores)

dadaFs <- dada(derepFs, err=errF, pool=TRUE, multithread=no_of_cores, 
               verbose=TRUE)
dadaRs <- dada(derepRs, err=errR, pool=TRUE, multithread=no_of_cores, 
               verbose=TRUE)

cat("dada-class: object describing DADA2 denoising results", sep="\n")

cat(paste(length(dadaFs[[1]]$denoised), 
          "sequence variants were inferred from", 
          length(derepFs[[1]]$uniques), 
          "input unique sequences."), sep="\n")

cat(paste("Key parameters: OMEGA_A =", dadaFs[[1]]$opts$OMEGA_A, 
          "OMEGA_C =", 
          dadaFs[[1]]$opts$OMEGA_C, "BAND_SIZE =", 
          dadaFs[[1]]$opts$BAND_SIZE), sep="\n")

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

table(nchar(getSequences(seqtab)))

# barplot(table(nchar(getSequences(seqtab))))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=no_of_cores, verbose=TRUE)

#View(seqtab.nochim)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x)) # getUniques() gets abundance of unique sequences
# Calculate the number of reads at different steps
track <- cbind(out, 
               sapply(dadaFs, getN), 
               sapply(dadaRs, getN), 
               sapply(mergers, getN), 
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track

track.long <- as.data.frame(track, 
                            row.names = row.names(track)) %>% 
  gather(., key = steps, value = counts, input:nonchim, factor_key = TRUE)

Datatracking<-ggplot(track.long, aes(x = steps, y = counts, color = steps)) +
  theme_classic() +
  geom_boxplot() +
  geom_jitter(shape = 16, position = position_jitter(0.3)) +
  scale_y_continuous(labels = scales::comma)
ggsave("/work/benson/bpeng4/Givaudan/Plot/Datatracking.png", plot = Datatracking, width = 6, height = 4, dpi = 300)

dir.create("./intermediate/", showWarnings = FALSE)
save(metadata_df, seqtab.nochim, file = "./intermediate/Givaudian.rda")
