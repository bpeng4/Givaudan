#
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
  library(ANCOMBC)
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

#Check Total Reads to find the Thresholds for Rarefication
OTU<-otu_table(ps.clean.p0) |>as.data.frame()
OTU$Total<- rowSums(OTU)

#Subset for each Sequencing run and output the otu table
ps.run1<- subset_samples(ps.clean.p0, Sample_Project=="Givaudan Part One" )
ps.run2<- subset_samples(ps.clean.p0, Sample_Project=="Givaudan Part Two" )
ps.run3<- subset_samples(ps.clean.p0, Sample_Project=="Givaudan Part Three" )

#Subset phyloseq files 
ps.clean.p0 <- subset_samples(ps.clean.p0, Microbiome %in% c("S1","S2","S3","S4","S5",
                                                             "S6","S12","S13","S14","S16"))
#Rarefication
ps.rare<-rarefy_even_depth(ps.clean.p0, sample.size = 5000,
                           rngseed = 111, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
#Change to relative abundance
ps.clean.re <- transform_sample_counts(ps.rare, function(x) x / sum(x))
ps.sample.re<-subset_samples(ps.clean.re, !ps.clean.re@sam_data$Sample %in% c("FBB0","FBB16"))
######################


#####
##Ancombc(SS filter)
#Set "FBB0" as the reference level for the dunnet comparison
sample_data(ps.clean.p0)$Category <- factor(
  sample_data(ps.clean.p0)$Category,
  levels = c("FBB0", "Synergistic effect (Fiber)",   "Postbiotic",                  
             "Prebiotic control",  "Prebiotic fiber" ,"Vitamin + fiber",             
             "Vitamin", "Vitamin (Synergistic)",       
             "Vitamin (Synergistic effect)", "Polyphenols" )
)

output <- ancombc2(
  data = ps.clean.p0,
  assay_name = "counts",
  tax_level = "Genus",
  fix_formula = "Category",
  rand_formula = "Microbiome + Rep",
  p_adj_method = "holm",
  pseudo_sens = TRUE,
  prv_cut = 0.10,
  lib_cut = 1000,
  s0_perc = 0.05,
  group = "Category",  # <- now safe
  struc_zero = FALSE,
  neg_lb = FALSE,
  alpha = 0.05,
  n_cl = 1,
  verbose = FALSE,
  global = TRUE,
  pairwise = TRUE,
  dunnet = TRUE,
  trend = FALSE,
  iter_control = list(tol = 1e-5, max_iter = 20, verbose = FALSE),
  em_control = list(tol = 1e-5, max_iter = 100),
  lme_control = NULL,
  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
  trend_control = NULL
)

#####Test
output_simple <- ancombc2(
  data = ps.clean.p0,
  assay_name = "counts",
  tax_level = "Genus",
  fix_formula = "Category",
  rand_formula = NULL,
  p_adj_method = "holm",
  pseudo_sens = TRUE,
  prv_cut = 0.05,
  lib_cut = 500,
  s0_perc = 0.05,
  group = "Category",
  struc_zero = TRUE,
  neg_lb = FALSE,
  alpha = 0.05,
  n_cl = 1,
  verbose = TRUE,
  global = TRUE,
  pairwise = TRUE,
  dunnet = TRUE,
  trend = FALSE
)

dim(output_simple$feature_table)
head(output_simple$res_pair)
table(output_simple$res_pair$diff_abn)
summary(output_simple$res_pair$lfc)

######


res_pair = output$res_pair %>%
  rowwise() %>%
  mutate(species = strsplit(taxon, "_")[[1]][7]) %>%
  filter(species != "NA") %>%
  ungroup()


res_pair = output$res_pair %>%
  separate(taxon, into = paste0("part", 1:10), sep = "_", fill = "right") %>%
  mutate(species = part7) %>%
  filter(!is.na(species) & species != "NA")

df_zero = output$zero_ind %>%
  rowwise() %>%
  mutate(species = strsplit(taxon, "_")[[1]][7]) %>%
  filter(species != "NA") %>%
  mutate(idx = sum(across(matches("structural")))) %>%
  filter(idx %in% 1:2) %>%
  transmute(Taxon = taxon, Species = species,
            None = ifelse(`structural_zero (surgery_type = none)` == TRUE, "Absence", "Presence"),
            Ileocolonic = ifelse(`structural_zero (surgery_type = ileocolonic)` == TRUE, "Absence", "Presence"),
            Colectomy = ifelse(`structural_zero (surgery_type = colectomy)` == TRUE, "Absence", "Presence")) %>%
  arrange(None, Ileocolonic)

## Visualization
df_fig1 = res_pair %>%
  filter(p_surgery_typeileocolonic < 0.05 | 
           p_surgery_typecolectomy < 0.05 |
           p_surgery_typecolectomy_surgery_typeileocolonic < 0.05) %>%
  transmute(species,
            lfc1 = ifelse(p_surgery_typeileocolonic < 0.05, 
                          round(lfc_surgery_typeileocolonic, 2), 0),
            lfc2 = ifelse(p_surgery_typecolectomy < 0.05, 
                          round(lfc_surgery_typecolectomy, 2), 0),
            lfc3 = ifelse(p_surgery_typecolectomy_surgery_typeileocolonic < 0.05, 
                          round(lfc_surgery_typecolectomy_surgery_typeileocolonic, 2), 0)) %>%
  pivot_longer(lfc1:lfc3, names_to = "group", values_to = "value") %>%
  arrange(species)

df_fig2 = res_pair %>%
  filter(p_surgery_typeileocolonic < 0.05 | 
           p_surgery_typecolectomy < 0.05 |
           p_surgery_typecolectomy_surgery_typeileocolonic < 0.05) %>%
  transmute(species,
            lfc1 = ifelse(q_surgery_typeileocolonic < 0.05, "aquamarine3", "black"),
            lfc2 = ifelse(q_surgery_typecolectomy < 0.05, "aquamarine3", "black"),
            lfc3 = ifelse(q_surgery_typecolectomy_surgery_typeileocolonic < 0.05, "aquamarine3", "black")) %>%
  pivot_longer(lfc1:lfc3, names_to = "group", values_to = "color") %>%
  arrange(species)

df_fig3 = res_pair %>%
  filter(p_surgery_typeileocolonic < 0.05 | 
           p_surgery_typecolectomy < 0.05 |
           p_surgery_typecolectomy_surgery_typeileocolonic < 0.05) %>%
  transmute(species,
            lfc1 = ifelse(passed_ss_surgery_typeileocolonic == TRUE, "yes", "no"),
            lfc2 = ifelse(passed_ss_surgery_typecolectomy == TRUE, "yes", "no"),
            lfc3 = ifelse(passed_ss_surgery_typecolectomy_surgery_typeileocolonic == TRUE, "yes", "no")) %>%
  pivot_longer(lfc1:lfc3, names_to = "group", values_to = "asterisk") %>%
  arrange(species)

df_fig = df_fig1 %>%
  dplyr::left_join(df_fig2, by = c("species", "group")) %>%
  dplyr::left_join(df_fig3, by = c("species", "group")) %>%
  dplyr::mutate(txt = ifelse(color == "aquamarine3" & asterisk == "yes",
                             paste0("*", value), value))

df_fig$group = recode(df_fig$group, 
                      `lfc1` = "Ileocolonic - None",
                      `lfc2` = "Colectomy - None", 
                      `lfc3` = "Colectomy - Ileocolonic")
df_fig$group = factor(df_fig$group, 
                      levels = c("Ileocolonic - None", 
                                 "Colectomy - None", 
                                 "Colectomy - Ileocolonic"))

fig_heat_pair = df_fig %>%
  ggplot(aes(x = group, y = species, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = 0, limit = c(-6, 6),
                       name = NULL) +
  geom_text(aes(group, species, label = txt, color = color), size = 2) +
  scale_color_identity(guide = FALSE) +
  labs(x = NULL, y = NULL, title = NULL) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text = element_text(size = 4.5),
        legend.text = element_text(size = 5),
        axis.text.y = element_text(face = "italic"))





























#ANCOM-BC2 pattern analysis
res_trend = output$res_trend

df_fig_trend = res_trend %>%
  dplyr::filter(diff_abn == 1) %>%
  dplyr::mutate(lfc1 = round(lfc_bmioverweight, 2),
                lfc2 = round(lfc_bmilean, 2),
                color = ifelse(passed_ss == 1, "aquamarine3", "black")) %>%
  tidyr::pivot_longer(cols = lfc1:lfc2, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)

df_fig_trend$group = recode(df_fig_trend$group, 
                            `lfc1` = "Overweight - Obese",
                            `lfc2` = "Lean - Obese")
df_fig_trend$group = factor(df_fig_trend$group, 
                            levels = c("Overweight - Obese", 
                                       "Lean - Obese"))

lo = floor(min(df_fig_trend$value))
up = ceiling(max(df_fig_trend$value))
mid = (lo + up)/2
fig_trend = df_fig_trend %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to obese subjects") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(color = df_fig_trend %>%
                                     dplyr::distinct(taxon, color) %>%
                                     .$color))
fig_trend














