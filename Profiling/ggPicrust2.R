rm(list = ls())
gc()
suppressMessages({
  library("tidyverse")
  library("phyloseq")
  library("rstatix")
  library("vegan")
  library("picante")
  library("ggpubr")
  library(vegan)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(poppr)  # For AMOVA if genetic analysis is needed
  library(microbiomeutilities)
  library(viridis)
  library(patchwork)
  library(pheatmap)
  library(RColorBrewer)
})
no_of_cores = 16
setwd("/work/benson/bpeng4/Givaudan/Givaudan_All_16S")
load("intermediate/Givaudian_ps.rda")
ls()
ps

###Pre-arrangement on  the phyloseq file
################
ps.clean <- subset_taxa(ps, Kingdom == "Bacteria") %>%
  subset_taxa(!is.na(Phylum)) %>%
  subset_taxa(!Class %in% c("Chloroplast")) %>%
  subset_taxa(!Family %in% c("Mitochondria"))
ps.clean
#Filter Out Taxa Exist in At Least 25% samples
ps.clean.p0 <- filter_taxa(ps.clean, function (x) {sum(x > 0) >= 255}, prune=TRUE)
ps.clean.p0

#Difine Functinos
phyloseq_counts_to_df <- function(ps) {
  df <- as.data.frame(t(otu_table(ps)))
  df <- cbind("#OTU ID" = rownames(df), df)  # Cheeky way to add a comment to first line
  
  # Note: we use OTU ID above since at the moment this function is meanly used
  # to pipe files into picrust, which needs the 'OTU ID'.
  return(df)
}

write_counts_to_file <- function(ps, filename) {
  df <- phyloseq_counts_to_df(ps)
  
  write.table(df,
              file=filename,
              sep = "\t",
              row.names = FALSE)
}


write_seqs_to_file <- function(ps, filename) {
  df <- as.data.frame(refseq(ps))
  
  unlink(filename)  # make sure we delete before we concatenate below
  
  names <- rownames(df)
  for (i in 1:nrow(df)) {
    cat(paste0(">", names[i]), file=filename, sep="\n", append=TRUE)
    cat(df[i,1], file=filename, sep="\n", append=TRUE)
  }
  message(paste0("Wrote ", nrow(df), " ASV sequences to file ", filename))
}
#Prepare the inputs file Picrust2
dir.create("picrust2", showWarnings = FALSE)
write_counts_to_file(ps.clean.p0, filename = "picrust2/raw_counts.tsv")
write_seqs_to_file(ps.clean.p0, filename = "picrust2/seqs.fna")

###
#{bash}
module load picrust2/2.4
# Use an extra cd to go to your work directory (not shown)
cd picrust2

picrust2_pipeline.py -s seqs.fna -i raw_counts.tsv -p 2 -o results \
                     --stratified --verbose

cd /work/benson/bpeng4/Kemin/Kemin_All/picrust2/results

add_descriptions.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
-o EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
-o KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

add_descriptions.py -i pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
-o pathways_out/path_abun_unstrat_descrip.tsv.gz

# Decompress tsv files
find . -name "*desc*" -exec pigz -dkv \{\} \;
####








