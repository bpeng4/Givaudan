library(ggplot2)
library(ggpubr)
library(dplyr)
library(stringi)
library(stringr)
library(agricolae)



set1<-read.csv(file = "/Users/bopeng/Documents/Givaudan/Cache/FinalPlate1-12New.csv")
set2<-read.csv(file = "/Users/bopeng/Documents/Givaudan/Cache/Final_Givaudan_AiMS_Meta_Map-3.csv")
set3<- merge(set1, set2, by = "Sample_Num")

Givaudan_SCFA<- reshape(data = set3,
                      varying = c( "Acetate","Propionate","Isobutyrate", 
                                  "Butyrate","Isovalerate", "Valeric.Acid","sum.ace.pro.buty.", "sum.isobuty.isoval."), 
                      v.names = "SCFAs",
                      direction = "long")
Givaudan_SCFA$Type<-c(rep(x="Acetate",810),rep(x="Propionate",810),rep(x="Isobutyrate",810),rep(x="Butyrate",810),rep(x="Isovalerate",810),
                  rep(x="Valeric.Acid",810),rep(x="Total_SCFAs.",810),rep(x="Total_BCFAs.",810))

#Order the SCFAs levels
# Set desired order
desired_order <- c("Acetate", "Butyrate", "Propionate", "Valeric.Acid", "Isobutyrate", "Isovalerate",  "Total_SCFAs.",
                   "Total_BCFAs.")

# Make sure your SCFA type column is a factor with the desired order
Givaudan_SCFA$Type <- factor(Givaudan_SCFA$Type, levels = desired_order)

#Remove Spaces within sample names
Givaudan_SCFA$Sample_Name <- stri_trim_both(Givaudan_SCFA$Sample_Name)
Givaudan_SCFA$Sample_Name <- str_trim(str_replace_all(Givaudan_SCFA$Sample_Name, "\\p{Zs}+", ""))

#Order Samples by the sample names
sample_order <- c("Super B-glucan (SBG)","Oat B-glucans (OBG 70% (Low m.wt))- Garuda",
                  "lantamanen OBG-29% GF","OBG 28% (OatWell Bran)",
                  "Yeast B-glucans (YBG-Wellmune)","Gingest",
                  "Inulin","AXOS", "Agrifiber", "Acerola full spectrum", 
                  "Red Acerola 20% Vit C", "Green Acerola 34% Vit C", 
                  "Acerola red 20% vit C& acerola green vit C 34%", "Ascorbic acid (Vit c)", 
                  "Carrot juice + Green Acerola","Carrot juice pro vit A", 
                  "Acerola green + OBG 28%", "Resistant starch postbiotic candidate 1",
                  "Resistant starch postbiotic candidate 2")

# Set factor levels
Givaudan_SCFA$Sample <- factor(Givaudan_SCFA$Sample, levels = sample_order)

# Arrange data frame by the desired order
Givaudan_SCFA <- Givaudan_SCFA %>% arrange(Sample_Abbr)

#Remove bad data
Givaudan_SCFA_clean <- Givaudan_SCFA %>%
  filter(!is.na(SCFAs), !is.nan(SCFAs), !is.infinite(SCFAs))

###Anova Test
#anova_model <- aov(SCFAs ~ Type * Microbiome * Sample *Rep, data = Givaudan_SCFA_clean)
#summary(anova_model)

#Plot and with anova test for SCFAs comparison among Samples
# Step 1: Calculate 95th percentile for each Type
quantile_limits <- Givaudan_SCFA_clean %>%
  group_by(Type) %>%
  summarise(q95 = quantile(SCFAs, 0.95, na.rm = TRUE))

# Step 2: Join back to original data and filter
Givaudan_SCFA_trimmed <- Givaudan_SCFA_clean %>%
  left_join(quantile_limits, by = "Type") %>%
  filter(SCFAs <= q95)

# Step 3: Plot with trimmed data
ggplot(Givaudan_SCFA_trimmed, aes(x = Sample, y = SCFAs, fill = Sample)) +
  geom_boxplot() +
  facet_wrap(~Type, scales = "free_y") +
  stat_compare_means(method = "anova", size = 3, label.x = 2.5) +
  labs(
    x = "SCFAs Type", 
    y = "Concentration (mM)", 
    title = "SCFAs Concentration by Sample Type"
  ) +
  theme(
    axis.title = element_text(size = 16),
    axis.text.x = element_text(size = 5, angle = 90),
    axis.text.y = element_text(size = 13),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 17)
  )

#for each subset of microbiomes
 # Define plotting function
 Sample_SCFAs <- function(Givaudan_SCFA, Microbiome) {
   # Create plot
   ggplot(Givaudan_SCFA, aes(x = Sample, y = SCFAs, fill = Sample)) +
     geom_boxplot() +
     facet_wrap(~Type, scales = "free_y") +
     stat_compare_means(method = "anova", size = 3, label.x = 2.5) +
     labs(
       x = "SCFAs Type", 
       y = "Concentration (mM)", 
       title = paste0(Microbiome, " SCFAs Concentration by Sample Type")
     ) +
     theme(
       axis.title = element_text(size = 16),
       axis.text.x = element_text(size = 5, angle = 90),
       axis.text.y = element_text(size = 13),
       legend.title = element_text(size = 18),
       legend.text = element_text(size = 16),
       plot.title = element_text(size = 17)
     )
 }
 
 # Ensure `Givaudan_SCFA` is defined if using it for unique Microbiome extraction
 Microbiome_list <- unique(Givaudan_SCFA_trimmed$Microbiome)
 
 # Loop over each microbiome group
 for (i in Microbiome_list) {
   Sub_Givaudan_SCFA <- Givaudan_SCFA_trimmed[Givaudan_SCFA_trimmed$Microbiome == i, ]
   print(Sample_SCFAs(Sub_Givaudan_SCFA, i))
 }
 

## Do the Duncan Test for Acetate
# Mean and SE by treatment
 Givaudan_Acetate <- Givaudan_SCFA_trimmed[Givaudan_SCFA_trimmed$Type == "Acetate",]
 
 aggregate(SCFAs ~ Sample, Givaudan_Acetate, function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))
 
 # ANOVA + Duncan
 model <- aov(SCFAs ~ Sample, data = Givaudan_Acetate)
 summary(model)
 duncan_result <- duncan.test(model, "Sample", console = TRUE)
 
 write.csv(duncan_result[4], file = "/Users/bopeng/Documents/GitHub/Givaudan/Plot/Acetate_Duncan_Mean_SE.csv" )
 write.csv(duncan_result[6], file = "/Users/bopeng/Documents/GitHub/Givaudan/Plot/Acetate_Post_Hoc_Comparison.csv" )
 
 
## Write A function and do for loop for each SCFA Duncan Test
SCFAs_Duncan <- function(Givaudan_SCFA_sub, i) {
  # ANOVA + Duncan
  model <- aov(SCFAs ~ Sample, data = Givaudan_SCFA_sub)
  summary(model)
  duncan_result <- duncan.test(model, "Sample", console = TRUE)
  
  write.csv(duncan_result[4], file = paste("/Users/bopeng/Documents/GitHub/Givaudan/Plot/", i, "_Duncan_Mean_SE.csv") )
  write.csv(duncan_result[6], file = paste("/Users/bopeng/Documents/GitHub/Givaudan/Plot/", i, "_Post_Hoc_Comparison.csv") )
  
}

Givaudan_SCFAs <- unique(Givaudan_SCFA_trimmed$Type)

for (i in Givaudan_SCFAs) {
  Givaudan_SCFA_sub <- Givaudan_SCFA_trimmed[Givaudan_SCFA_trimmed$Type == i,]
  SCFAs_Duncan(Givaudan_SCFA_sub, i)
}












 
 