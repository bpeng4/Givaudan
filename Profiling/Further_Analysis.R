#Load Package
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

#Remove Fs
ps <- subset_samples(ps , Category != "FS")
#Remove Spaces
# Remove leading/trailing whitespace from all levels in the Category column
ps@sam_data$Category <- trimws(ps@sam_data$Category)

###Pre-arrangement on  the phyloseq file
################
ps.clean <- subset_taxa(ps, Kingdom == "Bacteria") %>%
  subset_taxa(!is.na(Phylum)) %>%
  subset_taxa(!Class %in% c("Chloroplast")) %>%
  subset_taxa(!Family %in% c("Mitochondria"))
ps.clean
#ps.clean.re <- transform_sample_counts(ps.clean, function(x) x / sum(x))
#Filter Out Taxa Exist in At Least 25% samples
ps.clean.p0 <- filter_taxa(ps.clean, function (x) {sum(x > 0) >= 251}, prune=TRUE)
ps.clean.p0

#Check Total Reads to find the Thresholds for Rarefication
OTU<-otu_table(ps.clean.p0) |>as.data.frame()
OTU$Total<- rowSums(OTU)

#Rarefication
ps.rare<-rarefy_even_depth(ps.clean.p0, sample.size = 6138,
                           rngseed = 111, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
#Change to relative abundance
ps.clean.re <- transform_sample_counts(ps.clean.p0, function(x) x / sum(x))
ps.sample.re<-subset_samples(ps.clean.re, !ps.clean.re@sam_data$Sample %in% c("FBB0","FBB16"))
######################

###Plot Stack Bar chart for Baseline relative abundance
####################
ps.re.FBB<-subset_samples(ps.clean.re, ps.clean.re@sam_data$Sample %in% c("FBB0","FBB16"))

#Set up mycolor
mycolor<-c(
  '#F781BF', '#CCEBC5', '#FFFF99', '#FF7F00', '#8DD3C7', '#7FC97F', '#999999',  
  '#CAB2D6', '#FDBF6F', '#D95F02', '#666666', '#A6761D', '#B2DF8A',  
  '#6A3D9A', '#FB9A99', '#E41A1C', '#FFFFB3', '#FDB462', '#F0027F', '#D9D9D9', '#4DAF4A', 
  '#FFFF33', '#B15928', '#FDC086', '#80B1D3', '#A6CEE3', '#BC80BD', '#1B9E77', '#E31A1C',  
  '#FCCDE5', '#386CB0', '#377EB8', '#984EA3', '#FFED6F', '#66A61E', '#E6AB02',  
  '#BF5B17', '#1F78B4', '#A65628', '#B3DE69', '#7570B3', '#E7298A', '#33A02C','#FB8072',
  '#f7754f', '#ce9032', '#97a431', '#32b166', '#35ad9c', '#38a9c5', '#3ca2f4','#a48cf4', 
  '#f45cf2', '#f66bad')

sample_order <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8",
                  "S9", "S10", "S11", "S12", "S13", "S14", "S15", "S16")

ps.re.FBB@sam_data$Microbiome <- factor(ps.re.FBB@sam_data$Microbiome, levels = sample_order)


#Plot for Microbiome Baseline Profile
phyloseq::plot_bar(ps.re.FBB, fill = "Family") + 
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Microbiome, scales = "free_x",ncol = 4) +
  theme(axis.text.x = element_text(angle = 30, hjust = 0.5, size = 8),
        legend.key.size = unit(17, "pt"))+
  scale_fill_manual(values=mycolor)+
  guides(fill = guide_legend(ncol = 1))

########################
####Alpha Diversity Analysis and Plotting
## alpha diversity should be calculated before filtering on abundance and prevalence
#Remove NCs
ps.nanc<-subset_samples(ps, !ps@sam_data$Category %in% c("FBB0", "FBB16") )
tree = phyloseq::phy_tree(ps.nanc)
samp = data.frame(phyloseq::otu_table(ps.nanc))

##Alpha Diversity Among the Categories
adiv <- data.frame(
  phyloseq::estimate_richness(ps.nanc, measures = c( "Shannon", "Chao1", "Simpson", "InvSimpson", "Fisher")),
  "PD" = picante::pd(samp, tree, include.root=FALSE)[,1],
  dplyr::select(as_tibble(phyloseq::sample_data(ps.nanc)), Sample_Number, Sample_Name, 
                Sample_Project, Microbiome, Rep, Sample_Abbr, Sample,Category)) %>%
  dplyr::select(-se.chao1)
#Glance the median values
adiv %>%
  group_by( Category ) %>%
  dplyr::summarise(median_Chao1 = median(Chao1),
                   median_Shannon = median(Shannon),
                   median_Simpson = median(Simpson),
                   median_InvSimpson = median(InvSimpson),
                   median_Fisher = median(Fisher),
                   median_PD = median(PD))
#Plot
adivglm<-glm(Simpson ~ Category + Microbiome + Rep, data = adiv)
Anova(adivglm, type = "III")

# Plotting
adiv %>%
  pivot_longer(cols = c("Shannon", "Chao1", "Simpson", "PD"), 
               names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("Shannon", "Chao1", "Simpson", "PD"))) %>%
  ggplot(aes(x = Category, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Category), height = 0, width = .2) +
  labs(
    title = "Alpha Diversity Measures",
    x = "", 
    y = "Alpha Diversity Measures"
  )+
  facet_wrap(~ metric, scales = "free", ncol = 1) +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 14),
        title = element_text(size = 18)) +
  stat_anova_test()



AdivPlot<-function(adiv, title_text){
  adiv %>%
    pivot_longer(cols = c("Shannon", "Chao1", "Simpson", "PD"), 
                 names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = c("Shannon", "Chao1", "Simpson", "PD"))) %>%
    ggplot(aes(x = Category, y = value)) +
    geom_boxplot(outlier.color = NA) +
    geom_jitter(aes(color = Category), height = 0, width = .2) +
    labs(
      title = title_text,
      x = "", 
      y = "Alpha Diversity Measures"
    )+
    labs(x = "", y = "Alpha Diversity Measures") +
    facet_wrap(~ metric, scales = "free", ncol = 1) +
    theme(legend.position = "none", 
          axis.text.x = element_text(size = 10), 
          axis.text.y = element_text(size = 10),
          strip.text = element_text(size = 14),
          title = element_text(size = 18)) +
    stat_anova_test()
}

# Loop with dynamic titles
for (i in c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8",
            "S9", "S10", "S11", "S12", "S13", "S14", "S15", "S16")) {
  
  subset_data <- adiv[adiv$Microbiome == i, ]
  plot_title <- paste(i, "Alpha Diversity by Category")
  print(AdivPlot(subset_data, plot_title))
}
 


####Beta Diversity Analysis and Plotting
########################
ps.beta= ps.clean.re

###Principle Component Analysis
ordBC <- ordinate(ps.beta, "PCoA", "bray")
ordJC <- ordinate(ps.beta, "PCoA", "jaccard")
ordUF <- ordinate(ps.beta, "PCoA", "unifrac")
ordwUF <- ordinate(ps.beta, "PCoA", "wunifrac")
smpID <- ps.beta@sam_data$Sample_Name

# add sample_data info
BC <- merge(data.frame(ordBC$vectors[,1:2], sample = smpID, method = 'BC'), data.frame(sample_data(ps.beta)), by = "row.names")
rownames(BC) <- BC$Row.names
BC$Row.names <- NULL  # Remove the extra column created by merge

Jaccard <- merge(data.frame(ordJC$vectors[,1:2], sample = smpID,method = 'Jaccard'), data.frame(sample_data(ps.beta)), by = "row.names")
rownames(Jaccard) <- Jaccard$Row.names
Jaccard$Row.names <- NULL  # Remove the extra column created by merge

unifrac <- merge(data.frame(ordUF$vectors[,1:2], sample = smpID,method = 'unifrac'), data.frame(sample_data(ps.beta)), by = "row.names")
rownames(unifrac) <- unifrac$Row.names
unifrac$Row.names <- NULL  # Remove the extra column created by merge

wunifrac <- merge(data.frame(ordwUF$vectors[,1:2], sample = smpID,method = 'wunifrac'), data.frame(sample_data(ps.beta)), by = "row.names")
rownames(wunifrac) <- wunifrac$Row.names
wunifrac$Row.names <- NULL  # Remove the extra column created by merge

# Keep first 2 vectors (latent variables, PCs) of each distance matrix
df <- rbind(BC,Jaccard,unifrac,wunifrac)

#Form a column named MicroCatg
df$MicroCatg <- paste(df$Microbiome, df$Category, sep = "_")
# Calculate the mean for each level of MicroCatg
summary_df <- df %>%
  group_by(method,Microbiome,Category, MicroCatg) %>%
  summarise(
    Axis.1 = mean(Axis.1, na.rm = TRUE),
    Axis.2 = mean(Axis.2, na.rm = TRUE)
  )
#Order Category. levels
summary_df$Category <- factor(summary_df$Category, 
                               levels = c( "Polyphenols","Postbiotic",
                                           "Prebiotic control", "Prebiotic fiber",
                                           "Synergistic effect (Fiber)", "Vitamin",
                                           "Vitamin (Synergistic effect)","Vitamin (Synergistic)",
                                           "Vitamin + fiber", "FBB16","FBB0"))
#Plot without arrows
yourcolor<-c(
  '#F781BF', '#FF7F00', '#8DD3C7', '#7FC97F',  
  '#A6761D',   
  '#6A3D9A', '#FB9A99', '#E41A1C', '#FDB462', '#F0027F', '#D9D9D9', '#4DAF4A', 
  '#FFFF33', '#B15928', '#FDC086', '#80B1D3', '#A6CEE3', '#BC80BD', '#1B9E77', '#E31A1C',  
  '#FCCDE5', '#386CB0', '#377EB8', '#984EA3', '#FFED6F', '#66A61E', '#E6AB02',  
  '#BF5B17', '#1F78B4', '#A65628', '#B3DE69', '#7570B3', '#E7298A', '#33A02C','#FB8072',
  '#f7754f', '#ce9032', '#97a431', '#32b166', '#35ad9c', '#38a9c5', '#3ca2f4','#a48cf4', 
  '#f45cf2', '#f66bad')

ggplot(data = summary_df, aes(Axis.1, Axis.2, color = Microbiome, shape = Category)) + 
  geom_point(size=3) + 
  facet_wrap(~method, scales = 'free') +
  scale_color_manual(values = yourcolor) +
  scale_shape_manual(values = c(2:11,19))  # Define shape values

#Set up FBB0 as the starting points
summary_df0<-summary_df
summary_df0<-summary_df0[summary_df0$Category!="FBB16",]
fbb0_df <- summary_df0 %>%
  filter(Category == "FBB0") %>%
  select(method, Microbiome, Axis.1, Axis.2) %>%
  mutate(Axis.1_start = Axis.1, Axis.2_start = Axis.2) %>%
  select(-Axis.1, -Axis.2)  # Remove old columns to avoid duplicates

summary_df0 <- summary_df0 %>%
  left_join(fbb0_df, by = c("Microbiome", "method"))  # Match on method & Microbiome

# Order the levels
sample_order <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8",
                  "S9", "S10", "S11", "S12", "S13", "S14", "S15", "S16")

summary_df0$Microbiome <- factor(summary_df0$Microbiome, levels = sample_order)


#Plot with arrows
ggplot(summary_df0, aes(Axis.1, Axis.2, color = Microbiome, shape = Category.x)) + 
  geom_point(size = 3) + 
  geom_segment(aes(x = Axis.1_start, y = Axis.2_start, 
                   xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +  # Arrows from FBB0 to other categories
  facet_wrap(~method, scales = 'free') +
  scale_color_manual(values = yourcolor) +
  scale_shape_manual(values = c(2:10,19)) +  # Adjust shapes as needed
  labs(shape = "Category", color = "Microbiome") +  # Rename legends
  theme_minimal()
#############

#Beta Diversity Plot with arrows for each microbiome
####################
# Define color palette
yourcolor <- c(
  '#F781BF', '#FF7F00', '#8DD3C7', '#7FC97F', '#A6761D', '#6A3D9A', '#FB9A99',
  '#E41A1C', '#FDB462', '#F0027F', '#D9D9D9', '#4DAF4A', '#FFFF33', '#B15928',
  '#FDC086', '#80B1D3', '#A6CEE3', '#BC80BD', '#1B9E77', '#E31A1C', '#FCCDE5',
  '#386CB0', '#377EB8', '#984EA3', '#FFED6F', '#66A61E', '#E6AB02', '#BF5B17',
  '#1F78B4', '#A65628', '#B3DE69', '#7570B3', '#E7298A', '#33A02C', '#FB8072',
  '#f7754f', '#ce9032', '#97a431', '#32b166', '#35ad9c', '#38a9c5', '#3ca2f4',
  '#a48cf4', '#f45cf2', '#f66bad'
)

# Define plotting function
BetaPlot <- function(sub_df, micro_name, color_palette) {
  ggplot(sub_df, aes(Axis.1, Axis.2, color = Microbiome, shape = Category.x)) + 
    geom_point(size = 3) + 
    geom_segment(aes(x = Axis.1_start, y = Axis.2_start, xend = Axis.1, yend = Axis.2),
                 arrow = arrow(length = unit(0.2, "cm")), color = "black") +
    facet_wrap(~method, scales = 'free') +
    scale_color_manual(values = yourcolor_named) +
    scale_shape_manual(values = c(2:10, 19)) +
    labs(shape = "Category", color = "Microbiome") +
    theme_minimal() +
    ggtitle(paste(micro_name, "Beta Diversity Shift from the Baseline"))
}

# Loop over each microbiome and generate plot
microbiomes <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8",
                 "S9", "S10", "S11", "S12", "S13", "S14", "S15", "S16")

yourcolor_named <- setNames(yourcolor[1:length(microbiomes)], microbiomes)


for (i in seq_along(microbiomes)) {
  micro_name <- microbiomes[i]
  sub_df <- summary_df0[summary_df0$Microbiome == micro_name, ]
  
  if (nrow(sub_df) > 0) {
    print(BetaPlot(sub_df, micro_name, yourcolor))
  } else {
    message(paste("No data for", micro_name, "- skipping plot."))
  }
}


#Beta Diversity Plot with arrows for each group of health and disease microbiomes
####################
summary_df0_K3K5<-summary_df0[summary_df0$Microbiome%in%c("K3","K5"),]
summary_df0_K3K7<-summary_df0[summary_df0$Microbiome%in%c("K3","K7"),]
summary_df0_K10K3<-summary_df0[summary_df0$Microbiome%in%c("K10","K3"),]
summary_df0_K10K5<-summary_df0[summary_df0$Microbiome%in%c("K10","K5"),]
summary_df0_K10K7<-summary_df0[summary_df0$Microbiome%in%c("K10","K7"),]
#Get the first four colors
first_four <- brewer.pal(n = 8, name = "Set2")[1:4]
ggplot(summary_df0_K3K5, aes(Axis.1, Axis.2, color = Microbiome, shape = Category..x)) + 
  geom_point(size = 3) + 
  geom_segment(aes(x = Axis.1_start, y = Axis.2_start, 
                   xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +  # Arrows from FBB0 to other categories
  facet_wrap(~method, scales = 'free') +
  scale_color_manual(values = first_four[c(2,3)]) +
  scale_shape_manual(values = c(2:7,19)) +  # Adjust shapes as needed
  labs(shape = "Category", color = "Microbiome") +  # Rename legends
  theme_minimal() +
  ggtitle(paste("K3&K5", "Beta Diversity Shift from the Baseline"))

ggplot(summary_df0_K3K7, aes(Axis.1, Axis.2, color = Microbiome, shape = Category..x)) + 
  geom_point(size = 3) + 
  geom_segment(aes(x = Axis.1_start, y = Axis.2_start, 
                   xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +  # Arrows from FBB0 to other categories
  facet_wrap(~method, scales = 'free') +
  scale_color_manual(values = first_four[c(2,4)]) +
  scale_shape_manual(values = c(2:7,19)) +  # Adjust shapes as needed
  labs(shape = "Category", color = "Microbiome") +  # Rename legends
  theme_minimal() +
  ggtitle(paste("K3&K7", "Beta Diversity Shift from the Baseline"))

ggplot(summary_df0_K10K3, aes(Axis.1, Axis.2, color = Microbiome, shape = Category..x)) + 
  geom_point(size = 3) + 
  geom_segment(aes(x = Axis.1_start, y = Axis.2_start, 
                   xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +  # Arrows from FBB0 to other categories
  facet_wrap(~method, scales = 'free') +
  scale_color_manual(values = first_four[c(2,4)]) +
  scale_shape_manual(values = c(2:7,19)) +  # Adjust shapes as needed
  labs(shape = "Category", color = "Microbiome") +  # Rename legends
  theme_minimal() +
  ggtitle(paste("K10&K3", "Beta Diversity Shift from the Baseline"))

ggplot(summary_df0_K10K5, aes(Axis.1, Axis.2, color = Microbiome, shape = Category..x)) + 
  geom_point(size = 3) + 
  geom_segment(aes(x = Axis.1_start, y = Axis.2_start, 
                   xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +  # Arrows from FBB0 to other categories
  facet_wrap(~method, scales = 'free') +
  scale_color_manual(values = first_four[c(2,4)]) +
  scale_shape_manual(values = c(2:7,19)) +  # Adjust shapes as needed
  labs(shape = "Category", color = "Microbiome") +  # Rename legends
  theme_minimal() +
  ggtitle(paste("K10&K5", "Beta Diversity Shift from the Baseline"))

ggplot(summary_df0_K10K7, aes(Axis.1, Axis.2, color = Microbiome, shape = Category..x)) + 
  geom_point(size = 3) + 
  geom_segment(aes(x = Axis.1_start, y = Axis.2_start, 
                   xend = Axis.1, yend = Axis.2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +  # Arrows from FBB0 to other categories
  facet_wrap(~method, scales = 'free') +
  scale_color_manual(values = first_four[c(2,4)]) +
  scale_shape_manual(values = c(2:7,19)) +  # Adjust shapes as needed
  labs(shape = "Category", color = "Microbiome") +  # Rename legends
  theme_minimal() +
  ggtitle(paste("K10&K7", "Beta Diversity Shift from the Baseline"))
#################

#Subset data frame with only Weighted Unifrace Method
df0<-df[,c("Axis.1","Axis.2","Microbiome","MicroCatg","Category","Sample","method")]
df0<-df0[df0$method=="wunifrac",]

####Boxplot of Distances for Each Microbiome
#Set a Function to Plot the Distances to FBB0
####################
process_microbiome <- function(df0, Si, mycolor) {

  # Subset for FBB0
  df_FBB0 <- df0 %>% filter(Sample == "FBB0" & Microbiome == Si)
  
  # Calculate the mean Axis.1 & Axis.2 values for FBB0
  FBB0_Axis1 <- mean(df_FBB0$Axis.1, na.rm = TRUE)
  FBB0_Axis2 <- mean(df_FBB0$Axis.2, na.rm = TRUE)
  
  # Remove FBB0 from the dataset
  df_Si <- df0 %>% filter(Sample != "FBB0" & Microbiome == Si)
  
  # Reset Axis.1 and Axis.2 based on the new coordinate axis origin (FBB0_Axis1, FBB0_Axis2)
  df_Si <- df_Si %>%
    mutate(Axis.1 = Axis.1 - FBB0_Axis1,
           Axis.2 = Axis.2 - FBB0_Axis2)
  
  # Calculate distances based on Axis.1 and Axis.2 values
  df_Si <- df_Si %>%
    mutate(distance = sqrt(Axis.1^2 + Axis.2^2))
  
  # Order Category levels
#  df_Si$Category. <- factor(df_Si$Category, 
#                            levels = c("Prebiotic fiber","Vitamin + fiber",             
#                                       "Vitamin","Vitamin (Synergistic)",       
#                                       "Vitamin (Synergistic effect)", "Polyphenols",                 
#                                       "Synergistic effect (Fiber)","Postbiotic",                  
#                                       "Prebiotic control", "FBB16" ))
  
  # Box plot: Distance by Category
  p1 <- ggplot(df_Si, aes(x = Category, y = distance, fill = Category)) + 
    geom_boxplot(outlier.shape = 16,   
                 outlier.color = "black", 
                 outlier.size = 4) +
    labs(x = "Category", y = "Distance", 
         title = paste(Si, "Category-Wise Distribution of Weighted UniFrac Distances to FBB0")) +
    theme_minimal() +
    scale_fill_manual(values = mycolor)+
    theme(title = element_text(size=15))
  
  print(p1)
  return(df_Si)
}


#Plot the Category and Sample Distances to FBB0
df_S1 <- process_microbiome(df0, "S1", mycolor)
for (Subject in microbiomes){
  process_microbiome(df0, Subject, mycolor)
}
####################


###Write a function to do Amova test to calculate the Fst and p value for overall and each pair of sample
#################
perform_amova_for_categories <- function(df_Ki) {
  
  # Compute overall AMOVA for all categories
  distance_matrix <- dist(df_Ki$distance)
  overall_amova <- adonis2(distance_matrix ~ Category., data = df_Ki, permutations = 999)
  
  overall_result <- data.frame(
    Category1 = "Overall",
    Category2 = "",
    Fs_value = overall_amova$F[1],
    p_value = overall_amova$`Pr(>F)`[1]
  )
  
  # Generate all unique pairs of categories
  category_pairs <- expand.grid(Category1 = unique(df_Ki$Category.), 
                                Category2 = unique(df_Ki$Category.)) %>%
    mutate(across(everything(), as.character)) %>%  # Convert factor to character
    filter(Category1 < Category2)  # Remove duplicate pairs (A,B) vs (B,A)
  
  # Initialize an empty results dataframe
  amova_results <- data.frame(Category1 = character(), Category2 = character(), 
                              Fs_value = numeric(), p_value = numeric())
  
  # Loop through each pair and perform AMOVA
  for (i in 1:nrow(category_pairs)) {
    c1 <- category_pairs$Category1[i]
    c2 <- category_pairs$Category2[i]
    
    # Subset only the two categories being compared
    df_pair <- df_Ki %>% filter(Category. %in% c(c1, c2))
    
    # Compute distance matrix for the pair
    distance_matrix_pair <- dist(df_pair$distance)
    
    # Run AMOVA for the pair
    amova_result <- adonis2(distance_matrix_pair ~ Category., data = df_pair, permutations = 999)
    
    # Extract Fs and p-value
    Fs_value <- amova_result$F[1]
    p_value <- amova_result$`Pr(>F)`[1]
    
    # Append results to dataframe
    amova_results <- rbind(amova_results, data.frame(Category1 = c1, 
                                                     Category2 = c2, 
                                                     Fs_value = Fs_value, 
                                                     p_value = p_value))
  }
  
  # Combine overall and pairwise results
  final_results <- rbind(overall_result, amova_results)
  
  return(final_results)
}


#Calcualte Fst and p value for each pair of categories for each microbiome
df_K3_amova_categories <- perform_amova_for_categories(df_K3)
write.csv(df_K3_amova_categories, file = "/work/benson/bpeng4/Kemin/Kemin_All/Plot/K3_amova_categories.csv")
df_K5_amova_categories <- perform_amova_for_categories(df_K5)
write.csv(df_K5_amova_categories, file = "/work/benson/bpeng4/Kemin/Kemin_All/Plot/K5_amova_categories.csv")
df_K7_amova_categories <- perform_amova_for_categories(df_K7)
write.csv(df_K7_amova_categories, file = "/work/benson/bpeng4/Kemin/Kemin_All/Plot/K7_amova_categories.csv")
df_K10_amova_categories <- perform_amova_for_categories(df_K10)
write.csv(df_K10_amova_categories, file = "/work/benson/bpeng4/Kemin/Kemin_All/Plot/K10_amova_categories.csv")
#################


###Write a function to do Amova test to calculate the Fst and p value for overall and each pair of sample
#################
perform_amova_for_samples <- function(df_Ki) {
  
  # Compute overall AMOVA for all samples
  distance_matrix <- dist(df_Ki$distance)
  overall_amova <- adonis2(distance_matrix ~ Sample, data = df_Ki, permutations = 999)
  
  overall_result <- data.frame(
    Sample1 = "Overall",
    Sample2 = "",
    Fs_value = overall_amova$F[1],
    p_value = overall_amova$`Pr(>F)`[1]
  )
  
  # Generate all unique pairs of samples
  sample_pairs <- expand.grid(Sample1 = unique(df_Ki$Sample), 
                              Sample2 = unique(df_Ki$Sample)) %>%
    mutate(across(everything(), as.character)) %>%  # Convert factor to character
    filter(Sample1 < Sample2)  # Remove duplicate pairs (A,B) vs (B,A)
  
  # Initialize an empty results dataframe
  amova_results <- data.frame(Sample1 = character(), Sample2 = character(), 
                              Fs_value = numeric(), p_value = numeric())
  
  # Loop through each pair and perform AMOVA
  for (i in 1:nrow(sample_pairs)) {
    s1 <- sample_pairs$Sample1[i]
    s2 <- sample_pairs$Sample2[i]
    
    # Subset only the two samples being compared
    df_pair <- df_Ki %>% filter(Sample %in% c(s1, s2))
    
    # Compute distance matrix for the pair
    distance_matrix_pair <- dist(df_pair$distance)
    
    # Run AMOVA for the pair
    amova_result <- adonis2(distance_matrix_pair ~ Sample, data = df_pair, permutations = 999)
    
    # Extract Fs and p-value
    Fs_value <- amova_result$F[1]
    p_value <- amova_result$`Pr(>F)`[1]
    
    # Append results to dataframe
    amova_results <- rbind(amova_results, data.frame(Sample1 = s1, 
                                                     Sample2 = s2, 
                                                     Fs_value = Fs_value, 
                                                     p_value = p_value))
  }
  
  # Combine overall and pairwise results
  final_results <- rbind(overall_result, amova_results)
  
  return(final_results)
}

#Calcualte Fst and p value for each pair of samples for each microbiome
df_K3_amova_sample <- perform_amova_for_samples(df_K3)
df_K5_amova_sample <- perform_amova_for_samples(df_K5)
df_K7_amova_sample <- perform_amova_for_samples(df_K7)
df_K10_amova_sample <- perform_amova_for_samples(df_K10)

write.csv(df_K3_amova_sample, file = "/work/benson/bpeng4/Kemin/Kemin_All/Plot/K3_amova_sample.csv")
write.csv(df_K5_amova_sample, file = "/work/benson/bpeng4/Kemin/Kemin_All/Plot/K5_amova_sample.csv")
write.csv(df_K7_amova_sample, file = "/work/benson/bpeng4/Kemin/Kemin_All/Plot/K7_amova_sample.csv")
write.csv(df_K10_amova_sample, file = "/work/benson/bpeng4/Kemin/Kemin_All/Plot/K10_amova_sample.csv")

#################


###Weighted Unifrac Distance and Weight Association 
###################
# Extract unique samples
unique_samples <- unique(df0$Sample)
# Remove specific values
Sample <- unique_samples[!unique_samples %in% c("FBB0", "FBB16")]

#Write the function to plot
WeiDis <- function(df_K) {
  plot_list <- list()  # Store all plots
  
  # Determine global y-axis range
  y_min <- min(df_K$distance, na.rm = TRUE)
  y_max <- max(df_K$distance, na.rm = TRUE)
  
  for (i in Sample) {
    sub_sample <- df_K[df_K$Sample == i, ]  # Filter data for each sample
    
    # Fit the linear model
    model <- lm(distance ~ Amount, data = sub_sample)
    
    # Extract coefficients
    intercept <- round(coef(model)[1], 3)
    slope <- round(coef(model)[2], 3)
    
    # Create the regression formula text
    formula_text <- paste0("y = ", slope, "x + ", intercept)
    
    # Generate plot
    plot <- ggplot(sub_sample, aes(x = Amount, y = distance)) +
      geom_jitter(width = 0.2, alpha = 0.6, size = 0.7) +  # Jitter to prevent overlap
      geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Linear regression line
      annotate("text", x = max(sub_sample$Amount, na.rm = TRUE), 
               y = y_max, label = formula_text, hjust = 1, vjust = 1, 
               size = 2, color = "red") +  # Add regression formula
      theme_minimal() +
      labs(title = paste(i),
           x = "Amount (mg)",
           y = "Distance") +
      theme(plot.title = element_text(hjust = 0.5, size =10),  # Decreased title font size
            axis.title.x = element_text(size = 8),  # Decreased x-axis label font size
            axis.title.y = element_text(size = 8)   # Decreased y-axis label font size
            ) +
      ylim(y_min, y_max)  # Fix y-axis range
    
    plot_list[[i]] <- plot  # Store the plot
  }
  
  # Arrange plots in a 4-row, 6-column grid
  final_plot <- wrap_plots(plot_list, nrow = 4, ncol = 6)
  
  print(final_plot)  # Print the combined plot
}

WeiDis(df_K3)
WeiDis(df_K5)
WeiDis(df_K7)
WeiDis(df_K10)



###################

###Single taxon and Weight Association 
###################
# Set a meta data 
ps.genus = phyloseq::tax_glom(ps.rare, taxrank = rank_names(ps.rare)[6])
genus_meta <- sample_data(ps.genus)
genus_otu <- otu_table(ps.genus) |> as.data.frame()
genus_taxa <- tax_table(ps.genus) |> as.data.frame()

select_genus_taxa <- c("Fusobacterium", "Bacillus", "Bifidobacterium", "Faecalibacterium", 
                       "Blautia", "Lachnospiraceae UCG-010", "Clostridium") #Clostridium perfringens is the only species in Clostridium in the data set
select_genus_ASV <- genus_taxa[genus_taxa$Genus %in% select_genus_taxa,]
select_genus_otu <- genus_otu[, rownames(select_genus_ASV)]
colnames(select_genus_otu) <- select_genus_ASV$Genus
genus_meta <- cbind(genus_meta, select_genus_otu)
genus_meta <- genus_meta[!genus_meta$Sample %in% c("FBB0", "FBB16"),]

# Extract unique samples
Sample <- unique(genus_meta$Sample)

#Write the function to plot
# Determine global y-axis range
WeiTaxa <- function(genus_meta, taxon) {
  plot_list <- list()  # Store all plots
  
  y_min <- min(genus_meta[[taxon]], na.rm = TRUE)
  y_max <- max(genus_meta[[taxon]], na.rm = TRUE)
  
  for (i in Sample) {
    sub_sample <- genus_meta[genus_meta$Sample == i, ]  # Filter data for each sample
    
    # Fit the linear model
    model <- lm(sub_sample[[taxon]] ~ sub_sample$Amount, data = sub_sample)
    
    # Extract coefficients
    intercept <- round(coef(model)[1], 3)
    slope <- round(coef(model)[2], 3)
    
    # Create the regression formula text
    formula_text <- paste0("y = ", slope, "x + ", intercept)
    
    # Generate plot
    plot <- ggplot(sub_sample, aes(x = Amount, y = !!sym(taxon))) +
      geom_jitter(width = 0.2, alpha = 0.6, size = 0.7) +  # Jitter to prevent overlap
      geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Linear regression line
      annotate("text", x = max(sub_sample$Amount, na.rm = TRUE), 
               y = y_max, label = formula_text, hjust = 1, vjust = 1, 
               size = 3, color = "red") +  # Add regression formula
      theme_minimal() +
      labs(title = paste(i),
           x = "Amount (mg)",
           y = "Abundance") +
      theme(plot.title = element_text(hjust = 0.5, size = 10),  # Decreased title font size
            axis.title.x = element_text(size = 8),  # Decreased x-axis label font size
            axis.title.y = element_text(size = 8)   # Decreased y-axis label font size
      ) +
      ylim(y_min, y_max)  # Fix y-axis range
    
    plot_list[[i]] <- plot  # Store the plot
  }
  
  # Arrange plots in a 4-row, 6-column grid
  final_plot <- wrap_plots(plot_list, nrow = 4, ncol = 6) + 
    plot_annotation(title = paste("K10 Dose-Response Relationship of", taxon),
                    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
  
  print(final_plot)  # Print the combined plot
}

#Subset for each Microbiome
genus_meta_K3 <- genus_meta[genus_meta$Microbiome=="K3",]
genus_meta_K5 <- genus_meta[genus_meta$Microbiome=="K5",]
genus_meta_K7 <- genus_meta[genus_meta$Microbiome=="K7",]
genus_meta_K10 <- genus_meta[genus_meta$Microbiome=="K10",]

# Apply the function for each taxon
for (taxon in select_genus_taxa) {
  WeiTaxa(genus_meta_K3, taxon)
}
for (taxon in select_genus_taxa) {
  WeiTaxa(genus_meta_K5, taxon)
}
for (taxon in select_genus_taxa) {
  WeiTaxa(genus_meta_K7, taxon)
}
for (taxon in select_genus_taxa) {
  WeiTaxa(genus_meta_K10, taxon)
}

##For Family Level
ps.family = phyloseq::tax_glom(ps.rare, taxrank = rank_names(ps.rare)[5])
family_meta <- sample_data(ps.family)
family_otu <- otu_table(ps.family) |> as.data.frame()
family_taxa <- tax_table(ps.family) |> as.data.frame()

select_family_taxa <- "Enterococcaceae"
select_family_ASV <- family_taxa[family_taxa$Family %in% select_family_taxa,]
select_family_otu <- family_otu[, rownames(select_family_ASV)] |> as.data.frame()
colnames(select_family_otu) <- select_family_ASV$Family
family_meta <- cbind(family_meta, select_family_otu)
family_meta <- family_meta[!family_meta$Sample %in% c("FBB0", "FBB16"),]

# Extract unique samples
Sample <- unique(family_meta$Sample)


WeiTaxafam <- function(genus_meta_K) {
  plot_list <- list()  # Store all plots
  
  y_min <- min(genus_meta_K$Enterococcaceae, na.rm = TRUE)
  y_max <- max(genus_meta_K$Enterococcaceae, na.rm = TRUE)
  
  for (i in Sample) {
    sub_sample <- genus_meta_K[genus_meta_K$Sample == i, ]  # Filter data for each sample
    
    # Fit the linear model
    model <- lm(Enterococcaceae ~ Amount, data = sub_sample)
    
    # Extract coefficients
    intercept <- round(coef(model)[1], 3)
    slope <- round(coef(model)[2], 3)
    
    # Create the regression formula text
    formula_text <- paste0("y = ", slope, "x + ", intercept)
    
    # Generate plot
    plot <- ggplot(sub_sample, aes(x = Amount, y = Enterococcaceae)) +
      geom_jitter(width = 0.2, alpha = 0.6, size = 0.7) +  # Jitter to prevent overlap
      geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Linear regression line
      annotate("text", x = max(sub_sample$Amount, na.rm = TRUE), 
               y = y_max, label = formula_text, hjust = 1, vjust = 1, 
               size = 3, color = "red") +  # Add regression formula
      theme_minimal() +
      labs(title = paste(i),
           x = "Amount (mg)",
           y = "Abundance") +
      theme(plot.title = element_text(hjust = 0.5, size =10),  # Decreased title font size
            axis.title.x = element_text(size = 8),  # Decreased x-axis label font size
            axis.title.y = element_text(size = 8)   # Decreased y-axis label font size
      ) +
      ylim(y_min, y_max)  # Fix y-axis range
    
    plot_list[[i]] <- plot  # Store the plot
  }
  
  # Arrange plots in a 4-row, 6-column grid and add the title
  final_plot <- wrap_plots(plot_list, nrow = 4, ncol = 6) +
    plot_annotation(title = "K10 Dose-Response Relationship of Enterococcaceae",
                    theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))
  
  print(final_plot)  # Print the combined plot
}


#Subset for each Microbiome
family_meta_K3 <- family_meta[family_meta$Microbiome=="K3",]
family_meta_K5 <- family_meta[family_meta$Microbiome=="K5",]
family_meta_K7 <- family_meta[family_meta$Microbiome=="K7",]
family_meta_K10 <- family_meta[family_meta$Microbiome=="K10",]

WeiTaxafam(family_meta_K3)
WeiTaxafam(family_meta_K5)
WeiTaxafam(family_meta_K7)
WeiTaxafam(family_meta_K10)


###########

###Parallel Coordinate plot
#############
# Create an empty dataframe to store slopes
sampleslopes <- data.frame(Microbiome = character(), Sample = character(), Taxon = character(), Slope = numeric(), stringsAsFactors = FALSE)

FourMicro <- c("K3", "K5", "K7", "K10")

# Loop through each microbiome separately
for (mb in FourMicro) {
  genus_meta0 <- genus_meta[genus_meta$Microbiome == mb,]  # Filter by microbiome
  
  # Extract unique samples for this microbiome
  Sample <- unique(genus_meta0$Sample_Abbrv)
  
  for (taxon in select_genus_taxa) {
    for (i in Sample) {
      sub_sample <- genus_meta0[genus_meta0$Sample_Abbrv == i, ]  # Filter data for each sample
      
      if (nrow(sub_sample) > 1) {  # Ensure there are enough points for a linear model
        model <- lm(sub_sample[[taxon]] ~ sub_sample$Amount, data = sub_sample)
        
        # Extract the slope coefficient
        slope <- round(coef(model)[2], 3)
        
        # Append results to the dataframe
        sampleslopes <- bind_rows(sampleslopes, data.frame(Microbiome = mb, Sample = i, Taxon = taxon, Slope = slope))
      }
    }
  }
}

# View the results
View(sampleslopes)


#Set up the clean table
sampleslopes_wide <- sampleslopes %>%
  pivot_wider(names_from = Sample, values_from = Slope)

write.csv(sampleslopes_wide, file ="/work/benson/bpeng4/Kemin/Kemin_All/sampleslopes_wide.csv" )

#Parallel Coordinate Plot
library(GGally)

ggparcoord(data = sampleslopes_wide,
           columns = 4:26,
           groupColumn = "Taxon",
           scale = "center",
           showPoints = TRUE) +
  scale_color_brewer(palette = "Set3") +
  theme(
    axis.title = element_text(size = 20, face = "bold"),  # Axis title font size
    axis.text = element_text(size = 10),                 # Axis values font size
    legend.title = element_text(size = 20, face = "bold"), # Legend title size
    legend.text = element_text(size = 17)                 # Legend text size
  )

#############



####Heatmap
#Heatmap for Category
####################
#glommate to Genus Level
ps.genus = phyloseq::tax_glom(ps.rare, taxrank = rank_names(ps.rare)[6])

#Plotting
heat.sample <- plot_taxa_heatmap(ps.genus,
                                 subset.top = 50,
                                 VariableA = "Category.",
                                 heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                                 transformation = "log10")

#############Calculate the average read counts for each category
# Extract OTU (abundance) table and metadata from phyloseq
otu_table_df <- as.data.frame(otu_table(ps.genus))
meta_df <- as.data.frame(sample_data(ps.genus))

# Add category information to the OTU table
otu_table_df$Category <- meta_df$Category.[match(rownames(otu_table_df), rownames(meta_df))]

# Aggregate by category (mean abundance per category)
otu_table_avg <- otu_table_df %>%
  group_by(Category) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

# Convert back to matrix and ensure row names are correct
otu_table_avg <- column_to_rownames(otu_table_avg, var = "Category") %>%
  as.matrix()

# Convert back to a phyloseq OTU table
otu_table_avg_phylo <- otu_table(otu_table_avg, taxa_are_rows = FALSE)

# Ensure tax_table matches the taxa in otu_table_avg_phylo
common_taxa <- intersect(rownames(tax_table(ps.genus)), colnames(otu_table_avg_phylo))
tax_table_filtered <- tax_table(ps.genus)[common_taxa, ]

# Create a new sample metadata table
meta_avg_phylo <- sample_data(data.frame(Category = rownames(otu_table_avg_phylo),
                                         row.names = rownames(otu_table_avg_phylo)))

# Set factor levels for Category in the desired order
meta_avg_phylo$Category <- factor(meta_avg_phylo$Category, 
                                  levels = c("Botanicals", "Complex fiber", "Gums",
                                             "Oligosaccharrides","Waste Stream",
                                             "b-glucan","FBB16","FBB0"))  # Replace with actual category names

# Create a new phyloseq object with matching OTU and tax tables
ps.genus.avg <- phyloseq(otu_table_avg_phylo, tax_table_filtered, meta_avg_phylo)

# Generate heatmap of top 50 taxa, now grouped by category
heat.category <- plot_taxa_heatmap(ps.genus.avg,
                                   subset.top = 50,
                                   VariableA = "Category",
                                   heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                                   transformation = "log10")

# Plot the heatmap
heat.category

####################

#Heatmap for Sample
##############
heatmapsample <- function(ps.genus) {
  library(phyloseq)
  library(dplyr)
  library(tibble)
  library(RColorBrewer)
  
  # Extract OTU (abundance) table and metadata from phyloseq
  otu_table_df <- as.data.frame(otu_table(ps.genus))
  meta_df <- data.frame(sample_data(ps.genus))
  
  # Add sample and category info to OTU table
  otu_table_df$Sample <- meta_df$Sample
  otu_table_df$Category <- meta_df$Category
  
  # Aggregate by sample and category (mean abundance per sample)
  otu_table_avg <- otu_table_df %>%
    group_by(Category, Sample) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
  
  # Set Sample as rownames before dropping it from columns
  otu_table_avg <- otu_table_avg %>%
    column_to_rownames(var = "Sample")
  
  # Remove 'Category' column, then convert to matrix
  otu_table_matrix <- otu_table_avg[, !colnames(otu_table_avg) %in% "Category"] %>%
    as.matrix()
  
  # Create a phyloseq OTU table object
  otu_table_avg_phylo <- otu_table(otu_table_matrix, taxa_are_rows = FALSE)
  
  # Prepare corresponding sample metadata
  meta_avg_phylo <- meta_df %>%
    distinct(Sample, .keep_all = TRUE)
  rownames(meta_avg_phylo) <- meta_avg_phylo$Sample
  meta_avg_phylo <- sample_data(meta_avg_phylo)
  
  # Match taxa between original and averaged OTU tables
  common_taxa <- intersect(rownames(tax_table(ps.genus)), colnames(otu_table_matrix))
  tax_table_filtered <- tax_table(ps.genus)[common_taxa, ]
  
  # Create a new phyloseq object
  ps.genus.avg <- phyloseq(otu_table_avg_phylo, tax_table_filtered, meta_avg_phylo)
  
  # Plot heatmap of top 50 taxa
  heat.sample <- plot_taxa_heatmap(
    ps.genus.avg,
    subset.top = 50,
    VariableA = c("Category"),
    heatcolors = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
    transformation = "log10"
  )
  
  return(heat.sample)
}


ps.genus.K3<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "K3")
heatmapsample(ps.genus.K3)
ps.genus.K5<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "K5")
heatmapsample(ps.genus.K5)
ps.genus.K7<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "K7")
heatmapsample(ps.genus.K7)
ps.genus.K10<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "K10")
heatmapsample(ps.genus.K10)
##############


###Permanova Calculation and Visualization
##############
#glommate to Genus Level
ps.genus = phyloseq::tax_glom(ps.rare, taxrank = rank_names(ps.rare)[6])

#####Write the function 
ps.filter<-function(ps.genus){
  #Extract OTU/ASV table and metadata:
  otu_table <- as.data.frame(otu_table(ps.genus))
  metadata <- data.frame(sample_data(ps.genus))  # Force conversion
  taxa_table<- as.data.frame(tax_table(ps.genus))

  #Use the adonis function (PERMANOVA) from vegan to analyze treatment(Sample) effects
  metadata$Sample <- as.factor(metadata$Sample)  # Ensure it's a factor
  adonis_results <- adonis2(otu_table ~ Sample, data = metadata, permutations = 999, method = "bray")
  
  #Extract p-values
  p_value <- adonis_results$`Pr(>F)`[1]  # Extract p-value for treatment
  
  #test each taxon individually:
  p_values <- apply(otu_table, 2, function(x) {
    df <- data.frame(x = x, Sample = metadata$Sample)
    fit <- aov(x ~ Sample, data = df)
    summary(fit)[[1]][["Pr(>F)"]][1]  # Extract p-value
  })
  
  #Adjust for multiple testing (FDR correction):
  p_values_adj <- p.adjust(p_values, method = "fdr")
  
  #Mark taxa as significant if p < 0.05
  signif_taxa <- names(p_values_adj[p_values_adj < 0.05])

  #Filter for only the significant taxa
  ps.genus.filtered <- prune_taxa(taxa_names(ps.genus) %in% signif_taxa, ps.genus)
  
  return(ps.genus.filtered)
}
###############

###Heatmap for significant taxa calculated by Permanova for each Subject 
#Subset for each Subject
ps.genus<-subset_samples(ps.genus, !ps.genus@sam_data$Category == "FBB16")

ps.genus.S1<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S1")
ps.genus.filtered.S1 <- ps.filter(ps.genus.S1)
heatmapsample(ps.genus.filtered.S1)

ps.genus.S2<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S2")
ps.genus.filtered.S2 <- ps.filter(ps.genus.S2)
heatmapsample(ps.genus.filtered.S2)

ps.genus.S3<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S3")
ps.genus.filtered.S3 <- ps.filter(ps.genus.S3)
heatmapsample(ps.genus.filtered.S3)

ps.genus.S4<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S4")
ps.genus.filtered.S4 <- ps.filter(ps.genus.S4)
heatmapsample(ps.genus.filtered.S4)

ps.genus.S5<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S5")
ps.genus.filtered.S5 <- ps.filter(ps.genus.S5)
heatmapsample(ps.genus.filtered.S5)

ps.genus.S6<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S6")
ps.genus.filtered.S6 <- ps.filter(ps.genus.S6)
heatmapsample(ps.genus.filtered.S6)

ps.genus.S7<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S7")
ps.genus.filtered.S7 <- ps.filter(ps.genus.S7)
heatmapsample(ps.genus.filtered.S7)

ps.genus.S8<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S8")
ps.genus.filtered.S8 <- ps.filter(ps.genus.S8)
heatmapsample(ps.genus.filtered.S8)

ps.genus.S9<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S9")
ps.genus.filtered.S9 <- ps.filter(ps.genus.S9)
heatmapsample(ps.genus.filtered.S9)

ps.genus.S10<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S10")
ps.genus.filtered.S10 <- ps.filter(ps.genus.S10)
heatmapsample(ps.genus.filtered.S10)

ps.genus.S11<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S11")
ps.genus.filtered.S11 <- ps.filter(ps.genus.S11)
heatmapsample(ps.genus.filtered.S11)

ps.genus.S12<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S12")
ps.genus.filtered.S12 <- ps.filter(ps.genus.S12)
heatmapsample(ps.genus.filtered.S12)

ps.genus.S13<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S13")
ps.genus.filtered.S13 <- ps.filter(ps.genus.S13)
heatmapsample(ps.genus.filtered.S13)

ps.genus.S14<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S14")
ps.genus.filtered.S14 <- ps.filter(ps.genus.S14)
heatmapsample(ps.genus.filtered.S14)

ps.genus.S15<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S15")
ps.genus.filtered.S15 <- ps.filter(ps.genus.S15)
heatmapsample(ps.genus.filtered.S15)

ps.genus.S16<-subset_samples(ps.genus, ps.genus@sam_data$Microbiome == "S16")
ps.genus.filtered.S16 <- ps.filter(ps.genus.S16)
heatmapsample(ps.genus.filtered.S16)


ps.genus.filtered <- ps.filter(ps.genus)
heatmapsample(ps.genus.filtered)

#ps.rare<-subset_samples(ps.rare, !ps.rare@sam_data$Category == "FBB16")

ps.family = phyloseq::tax_glom(ps.rare, taxrank = rank_names(ps.rare)[5])
ps.family.filtered <- ps.filter(ps.family)
heatmapsample(ps.family.filtered)

ps.order = phyloseq::tax_glom(ps.rare, taxrank = rank_names(ps.rare)[4])
ps.order.filtered <- ps.filter(ps.order)
heatmapsample(ps.order.filtered)

ps.class = phyloseq::tax_glom(ps.rare, taxrank = rank_names(ps.rare)[3])
ps.class.filtered <- ps.filter(ps.class)
heatmapsample(ps.class.filtered)


ps.genus.sample <- subset_samples(ps.genus, !ps.genus@sam_data$Category %in% c("FBB0","FS"))
#Modify mannually by add "1" to each value in ps.genus otu table
otu_table(ps.genus.sample) <- otu_table(ps.genus.sample) + 1
ps.genus.sample.filtered <- ps.filter(ps.genus.sample)
heatmapsample(ps.genus.sample.filtered)



##############
####Prebiotic Potential Index Calculation
##Define Beneficial Bacteria and Harmful Bacteria
Beneficial<-c("Akkermansia", "Anaerostipes", "Barnesiella", "Bifidobacterium",
              "Blautia", "Butyricicoccus", "Catenibacterium", "Coprococcus",
              "Eubacterium", "Faecalibacterium", "Fusicatenibacter", "Lactobacillus",
              "Megasphaera", "Oscillibacter", "Parabacteroides", "Prevotella",
              "Roseburia", "Ruminococcus")
Harmful <- c("Bilophila", "Desulfovibrio", "Escherichia-Shigella","Fusobacterium", 
             "Haemophilus", "Klebsiella", "Paraprevotella", "Parasutterella", 
             "Streptococcus", "Sutterella", "Veillonella")

####Glomerate the phyloseq files to Beneficial and Harmful Categories
########
#Write a function to classified taxa to Beneficial and Harmful Categories
TaxGlom<-function(ps.rare){
  ps.genus = phyloseq::tax_glom(ps.rare, taxrank = rank_names(ps.rare)[6])
  
  genus_level <- rank_names(ps.rare)[6]  # Usually "Genus"
  
  #Modify the tax_table to include a new "Group" column
  tax_df <- as.data.frame(tax_table(ps.genus))
  tax_df$Group <- ifelse(tax_df[[genus_level]] %in% Beneficial, "Beneficial",
                         ifelse(tax_df[[genus_level]] %in% Harmful, "Harmful", "Other"))
  
  #Reassign the modified taxonomy table back to the phyloseq object
  tax_table(ps.genus) <- tax_table(as.matrix(tax_df))
  
  # Rename the new column as a taxonomic rank (so tax_glom can use it)
  colnames(tax_table(ps.genus))[ncol(tax_table(ps.genus))] <- "Group"
  
  #Glom based on your new "Group" column
  ps.grouped <- tax_glom(ps.genus, taxrank = "Group")
  
  #(optional): Remove "Other" if you only want Beneficial + Harmful
  ps.grouped <- subset_taxa(ps.grouped, Group %in% c("Beneficial", "Harmful"))
  
  #Extract OTU (ASV) table and tax_table
  otu_df <- as.data.frame(t(otu_table(ps.grouped)))
  tax_df <- as.data.frame(tax_table(ps.grouped))
  
  #Add Group column to OTU table using tax_table
  otu_df$Group <- tax_df$Group
  
  #Sum abundance across samples by Group
  grouped_abund <- otu_df %>%
    group_by(Group) %>%
    summarise(across(everything(), sum))
  
  #Convert back to phyloseq components (optional)
  otu_grouped <- as.matrix(grouped_abund[, -1])
  rownames(otu_grouped) <- grouped_abund$Group
  
  otu_ps <- otu_table(otu_grouped, taxa_are_rows = TRUE)
  otu_ps <- t(otu_ps)
  tax_grouped <- tax_table(matrix(data = grouped_abund$Group, ncol = 1,
                                  dimnames = list(grouped_abund$Group, "Group")))
  
  #Create new phyloseq object
  ps.grouped_final <- phyloseq(otu_ps, tax_grouped, sample_data(ps.grouped))
  
  return(ps.grouped_final)
}

#########
ps.grouped_final<-TaxGlom(ps.rare)

ps.rareS1<-subset_samples(ps.rare, Microbiome == "S1")
ps.rareS2<-subset_samples(ps.rare, Microbiome == "S2")
ps.rareS3<-subset_samples(ps.rare, Microbiome == "S3")
ps.rareS4<-subset_samples(ps.rare, Microbiome == "S4")
ps.rareS5<-subset_samples(ps.rare, Microbiome == "S5")
ps.rareS6<-subset_samples(ps.rare, Microbiome == "S6")
ps.rareS7<-subset_samples(ps.rare, Microbiome == "S7")
ps.rareS8<-subset_samples(ps.rare, Microbiome == "S8")
ps.rareS9<-subset_samples(ps.rare, Microbiome == "S9")
ps.rareS10<-subset_samples(ps.rare, Microbiome == "S10")
ps.rareS11<-subset_samples(ps.rare, Microbiome == "S11")
ps.rareS12<-subset_samples(ps.rare, Microbiome == "S12")
ps.rareS13<-subset_samples(ps.rare, Microbiome == "S13")
ps.rareS14<-subset_samples(ps.rare, Microbiome == "S14")
ps.rareS15<-subset_samples(ps.rare, Microbiome == "S15")
ps.rareS16<-subset_samples(ps.rare, Microbiome == "S16")

ps.S1.grouped_final<-TaxGlom(ps.rareS1)
ps.S2.grouped_final<-TaxGlom(ps.rareS2)
ps.S3.grouped_final<-TaxGlom(ps.rareS3)
ps.S4.grouped_final<-TaxGlom(ps.rareS4)
ps.S5.grouped_final<-TaxGlom(ps.rareS5)
ps.S6.grouped_final<-TaxGlom(ps.rareS6)
ps.S7.grouped_final<-TaxGlom(ps.rareS7)
ps.S8.grouped_final<-TaxGlom(ps.rareS8)
ps.S9.grouped_final<-TaxGlom(ps.rareS9)
ps.S10.grouped_final<-TaxGlom(ps.rareS10)
ps.S11.grouped_final<-TaxGlom(ps.rareS11)
ps.S12.grouped_final<-TaxGlom(ps.rareS12)
ps.S13.grouped_final<-TaxGlom(ps.rareS13)
ps.S14.grouped_final<-TaxGlom(ps.rareS14)
ps.S15.grouped_final<-TaxGlom(ps.rareS15)
ps.S16.grouped_final<-TaxGlom(ps.rareS16)


#Calculate the average read counts for each category
#############
# Extract OTU (abundance) table and metadata from phyloseq
otu_table_df <- as.data.frame(otu_table(ps.grouped_final))
meta_df <- as.data.frame(sample_data(ps.grouped_final))

# Add category information to the OTU table
otu_table_df$Category <- meta_df$Category[match(rownames(otu_table_df), rownames(meta_df))]

# Aggregate by category (mean abundance per category)
otu_table_avg <- otu_table_df %>%
  group_by(Category) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))

# Convert back to matrix and ensure row names are correct
otu_table_avg <- column_to_rownames(otu_table_avg, var = "Category") %>%
  as.matrix()

# Convert back to a phyloseq OTU table
otu_table_avg_phylo <- otu_table(otu_table_avg, taxa_are_rows = FALSE)

# Ensure tax_table matches the taxa in otu_table_avg_phylo
common_taxa <- intersect(rownames(tax_table(ps.grouped_final)), colnames(otu_table_avg_phylo))
tax_table_filtered <- tax_table(ps.grouped_final)[common_taxa, ]

# Create a new sample metadata table
meta_avg_phylo <- sample_data(data.frame(Category = rownames(otu_table_avg_phylo),
                                         row.names = rownames(otu_table_avg_phylo)))

# Set factor levels for Category in the desired order
meta_avg_phylo$Category <- factor(meta_avg_phylo$Category, 
                                  levels = c("FBB0", "FBB16",                       
                                              "FS", "Postbiotic","Prebiotic control",           
                                              "Prebiotic fiber", "Synergistic effect (Fiber)",  
                                              "Vitamin","Vitamin (Synergistic effect)",
                                              "Vitamin (Synergistic)","Vitamin + fiber"
                                              ))  # Replace with actual category names

# Create a new phyloseq object with matching OTU and tax tables
ps.grouped.final.avg <- phyloseq(otu_table_avg_phylo, tax_table_filtered, meta_avg_phylo)
############
ps.grouped.final.avg

#####Calculate  PPIa for each Category level
OTU_Average<- ps.grouped.final.avg@otu_table |> as.data.frame()
OTU_Average$Total<-rowSums(OTU_Average)
OTU_Average$PPIa<-(OTU_Average$Beneficial- OTU_Average$Harmful)/OTU_Average$Total 

write.csv(OTU_Average, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_Category.csv")



#Write a funciton to calculate the average read counts for each sample
#############

SampleRead <- function(ps.grouped_final){
  # Extract OTU (abundance) table and metadata from phyloseq
  otu_table_df <- as.data.frame(otu_table(ps.grouped_final))
  meta_df <- as.data.frame(sample_data(ps.grouped_final))
  
  # Add sample information to the OTU table
  otu_table_df$Sample <- meta_df$Sample[match(rownames(otu_table_df), rownames(meta_df))]
  
  # Aggregate by sample (mean abundance per sample)
  otu_table_avg <- otu_table_df %>%
    group_by(Sample) %>%
    summarise(across(where(is.numeric), mean, na.rm = TRUE))
  
  # Convert back to matrix and ensure row names are correct
  otu_table_avg <- column_to_rownames(otu_table_avg, var = "Sample") %>%
    as.matrix()
  
  # Convert back to a phyloseq OTU table
  otu_table_avg_phylo <- otu_table(otu_table_avg, taxa_are_rows = FALSE)
  
  # Ensure tax_table matches the taxa in otu_table_avg_phylo
  common_taxa <- intersect(rownames(tax_table(ps.grouped_final)), colnames(otu_table_avg_phylo))
  tax_table_filtered <- tax_table(ps.grouped_final)[common_taxa, ]
  
  # Create a new sample metadata table
  meta_avg_phylo <- sample_data(data.frame(Sample = rownames(otu_table_avg_phylo),
                                           row.names = rownames(otu_table_avg_phylo)))
  
  # Set factor levels for Sample in the desired order
  meta_avg_phylo$Sample <- factor(meta_avg_phylo$Sample, 
                                  levels = c("beetroot powder", "BetaVia Pure WD",                     
                                             "carrot powder", "Cellulose",                           
                                             "dried apple pomace", "dried spinich powder",                
                                             "ginger powder", "Ground Beet Pulp","guar gum",                  
                                             "HiSmooth flax seed fiber", "Jerusalum articohoke powder (inulin)",
                                             "kappa carageenan ", "lemon powder",                       
                                             "locust bean gum", "Micronized cassava fiber",            
                                             "Mushroom chitosan", "Oat Bagasse hydrolysate-(B)",         
                                             "Oat Bagasse hydrolysate-(PBV)", "Oat Bagasse-Control",                 
                                             "Olive Flour", "pomegranate powder",                  
                                             "Psyllium Husk","xanthan gum","FBB16","FBB0"))  # Replace with actual category names
  
  # Create a new phyloseq object with matching OTU and tax tables
  ps.grouped.final.avg <- phyloseq(otu_table_avg_phylo, tax_table_filtered, meta_avg_phylo)
  
  OTU_Average<- ps.grouped.final.avg@otu_table |> as.data.frame()
  OTU_Average$Total<-rowSums(OTU_Average)
  OTU_Average$PPIa<-(OTU_Average$Beneficial- OTU_Average$Harmful)/OTU_Average$Total 
  
  return(OTU_Average)
}

############
OTU_Average <- SampleRead(ps.grouped_final)
OTU_Average_S1 <- SampleRead(ps.S1.grouped_final)
OTU_Average_S2 <- SampleRead(ps.S2.grouped_final)
OTU_Average_S3 <- SampleRead(ps.S3.grouped_final)
OTU_Average_S4 <- SampleRead(ps.S4.grouped_final)
OTU_Average_S5 <- SampleRead(ps.S5.grouped_final)
OTU_Average_S6 <- SampleRead(ps.S6.grouped_final)
OTU_Average_S7 <- SampleRead(ps.S7.grouped_final)
OTU_Average_S8 <- SampleRead(ps.S8.grouped_final)
OTU_Average_S9 <- SampleRead(ps.S9.grouped_final)
OTU_Average_S10 <- SampleRead(ps.S10.grouped_final)
OTU_Average_S11 <- SampleRead(ps.S11.grouped_final)
OTU_Average_S12 <- SampleRead(ps.S12.grouped_final)
OTU_Average_S13 <- SampleRead(ps.S13.grouped_final)
OTU_Average_S14 <- SampleRead(ps.S14.grouped_final)
OTU_Average_S15 <- SampleRead(ps.S15.grouped_final)
OTU_Average_S16 <- SampleRead(ps.S16.grouped_final)


#####Calculate  PPIa for each Category level


write.csv(OTU_Average, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_Sample.csv")
write.csv(OTU_Average_S1, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S1.csv")
write.csv(OTU_Average_S2, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S2.csv")
write.csv(OTU_Average_S3, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S3.csv")
write.csv(OTU_Average_S4, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S4.csv")
write.csv(OTU_Average_S5, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S5.csv")
write.csv(OTU_Average_S6, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S6.csv")
write.csv(OTU_Average_S7, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S7.csv")
write.csv(OTU_Average_S8, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S8.csv")
write.csv(OTU_Average_S9, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S9.csv")
write.csv(OTU_Average_S10, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S10.csv")
write.csv(OTU_Average_S11, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S11.csv")
write.csv(OTU_Average_S12, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S12.csv")
write.csv(OTU_Average_S13, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S13.csv")
write.csv(OTU_Average_S14, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S14.csv")
write.csv(OTU_Average_S15, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S15.csv")
write.csv(OTU_Average_S16, file = "/work/benson/bpeng4/Givaudan/Plot/Prebiotic_IDX_S16.csv")



#########GMHI Calculation
ps.species = phyloseq::tax_glom(ps.clean.re, taxrank = rank_names(ps.clean.re)[7])
#Output otu tables
OutOTU<-function(ps.run, filename = "otu_otuput"){
  OTU<-otu_table(ps.run) |>as.data.frame()
  Taxa<-tax_table(ps.run) |>as.data.frame()
  Taxa$Name<- apply(Taxa[,6:7], 1, function(x) paste(x, collapse = "_"))
  colnames(OTU)<-Taxa$Name
  OTU<-t(OTU) |> as.data.frame()
  return(OTU)
}

OTU<-OutOTU(ps.species)

write.csv(OTU, file = "/work/benson/bpeng4/Givaudan/species_relative_abundances.csv")
