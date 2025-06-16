setwd("/work/benson/bpeng4/76Seeds_Peptides/data")
suppressMessages({
library("tidyverse")
library("dada2")
library("gridExtra")
library("devtools")
})
no_of_cores = 16
metadata = "SampleMeta_Majumder_Peptide.csv"
metadata_df = read.csv(metadata)
metadata_df = metadata_df[,-1]
fnFs <- metadata_df$fq1

sample.names <- metadata_df$sample

filt_path <- paste0(getwd(), '/filtered') # don't change this 
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))

FORWARD_TRUNC <- 250 # determine from quality plots

out <- filterAndTrim(fnFs, filtFs,
                     truncLen=FORWARD_TRUNC, 
                     trimLeft=20, maxEE=2, 
                     multithread=FALSE,
                     matchIDs=TRUE, compress=TRUE, 
                     verbose=TRUE)

derepFs <- derepFastq(filtFs, n = 1e+06, verbose = TRUE)

names(derepFs) <- sample.names

errF <- learnErrors(filtFs, verbose=TRUE, multithread=no_of_cores)

dadaFs <- dada(derepFs, err=errF, pool=TRUE, multithread=no_of_cores, 
               verbose=TRUE)

seqtab <- makeSequenceTable(dadaFs)


seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread=no_of_cores, verbose=TRUE)

save(metadata_df, seqtab.nochim, file = "./intermediate/PiptideForward.rda")
