rm(list = ls())
gc()
setwd( "/work/benson/bpeng4/Givaudan/Givaudan_All_16S")
suppressMessages({
  library("tidyverse")
  library("dada2")
  library("gridExtra")
  library("phyloseq")
  library("Biostrings")
})
no_of_cores = 36
load("./intermediate/Givaudian.rda")
ls()
gc(full=TRUE) 
#silva <- "silva_nr99_v138.1_wSpecies_train_set.fa.gz" # don't change this
silva <- "silva_nr99_v138.2_toSpecies_trainset.fa.gz"  #Use the New Version of Silva
taxa <- assignTaxonomy(seqtab.nochim, silva, multithread=no_of_cores, verbose=TRUE)
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
# Name the sequence file
sequence_outfile <- "dada2_nochim.fa"
print(paste0('Writing FASTA file ', sequence_outfile, ' with ', ncol(seqtab.nochim), ' sequences.'))
file.remove(sequence_outfile) # Make sure the file does not exist
for (i in 1:ncol(seqtab.nochim)){
  cat(paste0(">ASV_", i), file = sequence_outfile, sep="\n", append = TRUE)
  cat(colnames(seqtab.nochim)[i], file = sequence_outfile, sep="\n", append= TRUE)
}
system("head dada2_nochim.fa", intern=TRUE) # use system() to run bash command in R

# We need to make the external packages MAFFT and FastTree available for our use;
# this is what the next three lines do.
source(file.path(Sys.getenv("LMOD_PKG"), "init/R"))
module("load", "mafft")
module("load", "fasttree")

# Multiple sequence alignment with MAFFT
system("mafft --auto --thread -1 dada2_nochim.fa > dada2_nochim_mafft_msa.fa", intern=TRUE)

# Phylogenetic tree reconstruction with FastTree
system("fasttree -gtr -nt < dada2_nochim_mafft_msa.fa > dada2_nochim.tree", intern=TRUE)

list.files(pattern = "*dada2*")
seqfile <- "dada2_nochim.fa"
treefile <- "dada2_nochim.tree"
ls()
seqtab.nochim[1:3,1:3]
taxa[1:3,] 
metadata_df[1:3,] 
seqtab.nochim.0 <- seqtab.nochim
taxa.0 <- taxa
metadata_df.0 <- metadata_df
colnames(seqtab.nochim) <- paste0("ASV_", seq(1:ncol(seqtab.nochim)))
rownames(taxa) <- paste0("ASV_", seq(1:nrow(taxa)))
seqtab.nochim[1:3,1:3]
taxa[1:3,]
metadata_df[1:3,]
rownames(metadata_df)<-metadata_df$Sample_Name

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(metadata_df),
               tax_table(taxa))
ps
# Add rep seqs
repseqs <- readDNAStringSet(seqfile)
ps <- merge_phyloseq(ps, repseqs)
ps
# Add tree
tree <- read_tree(treefile)
ps <- merge_phyloseq(ps, tree)
ps

save(ps, file = "intermediate/Givaudian_ps.rda")
