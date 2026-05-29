#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                         Honeybee gut microbiota analysis                     #
#                              Juan J. Ovalle-Barrera                          #
#                             Version 0.4 - 5/26/2026                          #
#                                                                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Contents
# 1. Data uploading, phyloseq and rarefaction
#   1.1 Subset creation
# 2. Alpha diversity analysis
#   2.1 Control vs. Pesticides
#   2.2 Control vs. Flavonoids
#   2.3 Pesticides vs. Pesticides+Flavonoids
# 3. Beta diversity analysis
#   3.1 Bray curtis distances
#   3.2 PERMANOVA
# 4. DeSeq, relative abundance analysis
# 5. Graphs 
#   5.1 Alpha diversity boxplot
#   5.2 NMDS
#   5.3 Relative abundance 

# Libraries
require("ggplot2")  
require("BiocManager")
require("dplyr")
require("vegan")
require("phyloseq")
require("tidyverse")
require("pheatmap")
require("ggsignif")

# 1. Data uploading, processing and rarefaction~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Data uploading
nochim <- read.csv("C:/Users/juanj/Documents/Abejas-Microbiota/Abejas-Microbiota/Analisis/seqtab_nochim_transp_1.csv", row.names = 1, header = TRUE)
taxa <- read.csv("C:/Users/juanj/Documents/Abejas-Microbiota/Abejas-Microbiota/Analisis/taxa_bee_16s.csv", row.names = 1)
metadata <- read.csv("C:/Users/juanj/Documents/Abejas-Microbiota/Abejas-Microbiota/Analisis/metadata_bee_16s.csv", row.names = 1)

# Matrixes
nochim_mat <- as.matrix(nochim)
taxa_mat <- as.matrix(taxa)

# Phyloseq database
ps <- phyloseq(otu_table(nochim_mat, taxa_are_rows = TRUE),
  tax_table(taxa_mat), sample_data(metadata)
)

ps

# Read counting
sample_sums(ps)

# Soting and read quantity checking
sort(sample_sums(ps))

# Phyla
table(tax_table(ps)[, "Phylum"])

# Rarefaction
ps_rare <- rarefy_even_depth(ps, rngseed = 123)

# 1.1 Subset creation from the main phyloseq

# Control and all individual molecules
ps_uni <- subset_samples(ps_rare, 
                         Treatment %in% c("Control_(N)", "Rutina", "Imidacloprid", "Kaempferol", 
                                          "p-Coumaric_acid", "Fipronil", "Sugar", "Deltametrina", "Glifosato"))

# Control and all individual pesticides
ps_pest <- subset_samples(ps_rare, 
                          Treatment %in% c("Control_(N)", "Imidaclorid", "Fipronil", "Sugar", "Deltametrina", "Glifosato"))

# Control and all individual flavonoids
ps_phyto <- subset_samples(ps_rare, 
                           Treatment %in% c("Control_(N)", "Rutina", "Kaempferol", 
                                            "p-Coumaric_acid", "Sugar"))

# to be continued...
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# 2. Alpha diversity analysis~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Alpha diversity for all treatments
alpha_div_all  <- estimate_richness(ps_rare, 
                               measures = c("Observed", 
                                            "Shannon", 
                                            "Simpson",
                                            "Chao1"))

alpha_div_all$Treatment <- as.character(metadata$Treatment)
alpha_div_all 

# Kruskall-wallis test for all treatments, non significant for any test
kruskal.test(Shannon ~ Treatment, data = alpha_div_all)
kruskal.test(Observed ~ Treatment, data = alpha_div_all)
kruskal.test(Chao1 ~ Treatment, data = alpha_div_all)

# Alpha diversity for each individual molecule
alpha_div_uni  <- estimate_richness(ps_uni, 
                                    measures = c("Observed", 
                                                 "Shannon", 
                                                 "Simpson",
                                                 "Chao1"))
alpha_div_uni$Treatment <- as.character(sample_data(ps_uni)$Treatment)
alpha_div_uni

# Kruskall-wallis test for each individual molecule, non significant ~
kruskal.test(Shannon ~ Treatment, data = alpha_div_uni)
kruskal.test(Observed ~ Treatment, data = alpha_div_uni)
kruskal.test(Chao1 ~ Treatment, data = alpha_div_uni)

# Alpha diversity for each individual pesticide
alpha_div_pest  <- estimate_richness(ps_pest, 
                                    measures = c("Observed", 
                                                 "Shannon", 
                                                 "Simpson",
                                                 "Chao1"))
alpha_div_pest$Treatment <- as.character(sample_data(ps_pest)$Treatment)
alpha_div_pest

# Kruskall-wallis test for each individual pesticide, not significant for any treatment
kruskal.test(Shannon ~ Treatment, data = alpha_div_pest)
kruskal.test(Observed ~ Treatment, data = alpha_div_pest)
kruskal.test(Chao1 ~ Treatment, data = alpha_div_pest)

# Alpha diversity for each individual phytochemical
alpha_div_phyto  <- estimate_richness(ps_phyto, 
                                     measures = c("Observed", 
                                                  "Shannon", 
                                                  "Simpson",
                                                  "Chao1"))
alpha_div_phyto$Treatment <- as.character(sample_data(ps_phyto)$Treatment)
alpha_div_phyto

# Kruskall-wallis test for each individual phytoicide, not significant for any treatment
kruskal.test(Shannon ~ Treatment, data = alpha_div_phyto)
kruskal.test(Observed ~ Treatment, data = alpha_div_phyto)
kruskal.test(Chao1 ~ Treatment, data = alpha_div_phyto)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# 3. Beta diversity analysis~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# 3.1 Bray curtis distances

# Calculate distances between individuals, Bray Curtis method
bray_all <- phyloseq::distance(ps_rare, method = "bray")

# Average distance per treatments
distances_all <- meandist(bray_all, sample_data(ps_rare)$Treatment)

View(distances_all)

# Calculate distances between individual molecules, Bray Curtis method
bray_uni <- phyloseq::distance(ps_uni, method = "bray")

# Average distance per treatments
distances_uni <- meandist(bray_uni, sample_data(ps_uni)$Treatment)

View(distances_uni)

# 3.2 PERMANOVA

# For all treatments
permanova_all <- adonis2(bray_all ~ Treatment,
                         data = as.data.frame(metadata_ordered),
                         permutations = 999)

# For all individual molecules
permanova_uni <- adonis2(bray_uni ~ sample_data(ps_uni)$Treatment,
                         data = as.data.frame(metadata),
                         permutations = 999)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Bray-Curtis NMDS
ord_nmds <- ordinate(ps_rare, method = "NMDS", distance = "bray")

# Stress for NMDS, less than 0.2
ord_nmds$stress

#Graph NMDS for all treatments simultaneously
nmds_points <- as.data.frame(ord_nmds$points)
nmds_points$Treatment <- sample_data(ps_rare)$Treatment

graph_nmds_all <- ggplot(nmds_points, aes(x = MDS1, y = MDS2, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +  
  stat_ellipse(aes(group = Treatment), 
        type = "t",level = 0.95)+
  labs(title = "NMDS - Bray Curtis",
       subtitle = paste("Stress =", round(ord_nmds$stress, 3))) +
  theme_bw()

ggsave("Figures/graph_nmds_all.jpg", plot = graph_nmds_all)

#PERMANOVA for all treatments

bray_dist <- phyloseq::distance(ps_rare, method = "bray")

permanova_all <- adonis2(bray_dist ~ Treatment,
                     data = as.data.frame(metadata),
                     permutations = 999)
permanova_all #Not significant when evaluating everything at once

sample_data(ps_rare)$Treated <- ifelse(
  sample_data(ps_rare)$Treatment == "Control_(N)", 
  "Control", 
  "Treated"
)
View(sample_data(ps_rare))

bray_dist <- phyloseq::distance(ps_rare, method = "bray")

adonis2(bray_dist ~ Treatment,
        data = as.data.frame(metadata),
        permutations = 999)

nrow(metadata)
nsamples(ps_rare)

#Trial subsets

#All individual molecules
ps_uni <- subset_samples(ps_rare, 
                         Treatment %in% c("Control_(N)", "Rutina", "Imidacloprid", "Kaempferol", 
                                          "p-Coumaric_acid", "Fipronil", "Sugar", "Deltametrina", "Glifosato"))

ord_nmds_uni <- ordinate(ps_uni, method = "NMDS", distance = "bray")
bray_uni <- phyloseq::distance(ps_uni, method = "bray")

nmds_points_uni <- as.data.frame(ord_nmds_uni$points)
nmds_points_uni$Treatment <- sample_data(ps_uni)$Treatment

permanova_uni <- adonis2(bray_uni ~ Treatment,
        data = data.frame(sample_data(ps_uni)),
        permutations = 999)

graph_nmds_uni <- ggplot(nmds_points_uni, aes(x = MDS1, y = MDS2, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +  
  stat_ellipse(aes(group = Treatment), 
               type = "t",level = 0.95)+
  labs(title = "NMDS - Bray Curtis",
       subtitle = paste("Stress =", round(ord_nmds_uni$stress, 3))) +
  theme_bw()

#Bray curtis distance represented
bray_matrix <- as.matrix(bray_uni)

# Agregar tratamiento a la matriz
treatments <- data.frame(sample_data(ps_uni))$Treatment

# Distancia promedio entre cada par de tratamientos
distances_all <- meandist(bray_uni, treatments)

View(distances_all)

dist_df <- as.data.frame(as.matrix(distances_all))

pheatmap::pheatmap(distances_all,
                   display_numbers = TRUE,
                   number_format = "%.2f",
                   main = "All treatments average Bray curtis distance")

#All pesticides

ps_pest <- subset_samples(ps_rare, 
                         Treatment %in% c("Control_(N)", "Imidaclorid", "Fipronil", "Sugar", "Deltametrina", "Glifosato"))

ord_nmds_pest <- ordinate(ps_pest, method = "NMDS", distance = "bray")
bray_pest <- phyloseq::distance(ps_pest, method = "bray")

nmds_points_pest <- as.data.frame(ord_nmds_pest$points)
nmds_points_pest$Treatment <- sample_data(ps_pest)$Treatment

adonis2(bray_pest ~ Treatment,
        data = data.frame(sample_data(ps_pest)),
        permutations = 999)

graph_nmds_pest <- ggplot(nmds_points_pest, aes(x = MDS1, y = MDS2, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +  
  stat_ellipse(aes(group = Treatment), 
               type = "t",level = 0.95)+
  labs(title = "NMDS - Bray Curtis",
       subtitle = paste("Stress =", round(ord_nmds_pest$stress, 3))) +
  theme_bw()

graph_nmds_pest

#Individual pesticides represented

bray_matrix_pest <- as.matrix(bray_pest)

# Agregar tratamiento a la matriz
treatments <- data.frame(sample_data(ps_pest))$Treatment

# Distancia promedio entre cada par de tratamientos
distances_pest <- meandist(bray_pest, treatments)

View(distances_pest)

dist_df_pest <- as.data.frame(as.matrix(distances_pest))

pheatmap::pheatmap(dist_df_pest,
                   display_numbers = TRUE,
                   number_format = "%.2f",
                   main = "Distancia Bray-Curtis entre pesticidas")

#All flavonoids

ps_phyto <- subset_samples(ps_rare, 
                         Treatment %in% c("Control_(N)", "Rutina", "Kaempferol", 
                                          "p-Coumaric_acid", "Sugar"))

ord_nmds_phyto <- ordinate(ps_phyto, method = "NMDS", distance = "bray")
bray_phyto <- phyloseq::distance(ps_phyto, method = "bray")

nmds_points_phyto <- as.data.frame(ord_nmds_phyto$points)
nmds_points_phyto$Treatment <- sample_data(ps_phyto)$Treatment

adonis2(bray_phyto ~ Treatment,
        data = data.frame(sample_data(ps_phyto)),
        permutations = 999)

graph_nmds_phyto <- ggplot(nmds_points_phyto, aes(x = MDS1, y = MDS2, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +  
  stat_ellipse(aes(group = Treatment), 
               type = "t",level = 0.95)+
  labs(title = "NMDS - Bray Curtis",
       subtitle = paste("Stress =", round(ord_nmds_phyto$stress, 3))) +
  theme_bw()

graph_nmds_phyto

#Pesticides vs pesticides + flavonoids

ps_pesticides <- subset_samples(ps_rare,
                           Treatment %in% c("Imidacloprid", "Fipronil", "Deltametrina", "Glifosato", "Kaempferol+Deltametrina", "Kaempferol+Fipronil", "Kaempferol+Glifosato", "Kaempferol+Imidaclorid", "	
p-Coumaric_acid+Deltametrina", "p-Coumaric_acid+Fipronil", "p-Coumaric_acid+Glifosato", "p-Coumaric_acid+Imidaclorid", "Rutina+Deltametrina", "Rutina+Fipronil", "	
Rutina+Glifosato", "Rutina+Imidaclorid")) 

bray_pesticides <- phyloseq::distance(ps_pesticides, method = "bray")

adonis2(bray_pesticides ~ Treatment,
        data = data.frame(sample_data(ps_pesticides)),
        permutations = 999)

#Per molecule analisis

ps_mol1 <- subset_samples(ps_rare,
                          Treatment %in% c("Control", "mol1"))

bray_mol1 <- phyloseq::distance(ps_mol1, method = "bray")

adonis2(bray_mol1 ~ Treatment,
        data = data.frame(sample_data(ps_mol1)),
        permutations = 999)

# Filtrar taxa con abundancia relativa promedio > 1%
ps_fil <- filter_taxa(ps_uni,
                      function(x) mean(x/sum(x)) > 0.01 & 
                        sum(x > 0) > 0.2 * nsamples(ps_uni),
                      TRUE)

ps_phylum <- tax_glom(ps_fil, taxrank = "Phylum")
plot_bar(ps_phylum, fill = "Phylum", x = "Treatment") +
  theme_bw()

# A nivel de Género
ps_genus <- tax_glom(ps_fil, taxrank = "Genus")
plot_bar(ps_genus, fill = "Genus", x = "Treatment") +
  theme_bw()


# Calcular métricas de diversidad alfa
alpha_div <- estimate_richness(ps_uni, 
                               measures = c("Observed", 
                                            "Shannon", 
                                            "Simpson",
                                            "Chao1"))

# Agregar metadata
alpha_div$Treatment <- sample_data(ps_uni)$Treatment
alpha_div$Sample <- rownames(alpha_div)

View(alpha_div)

alpha_div %>%
  pivot_longer(cols = c(Observed),
               names_to = "metric",
               values_to = "value") %>%
  ggplot(aes(x = Treatment, y = value, fill = Treatment)) +
  geom_boxplot(alpha = 0.7) +
  geom_signif(
    comparisons = list(
      c("Control_(N)", "Rutina"), c("Rutina", "Deltametrina")
    ),map_signif_level = TRUE) +
  geom_jitter(width = 0.2, size = 1.5) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

pairwise_adonis(
  ps_uni,
  Treatment,
  metadata = NULL,
  perm_design = NULL,
  p.adjust.method = "bonferroni",
  perm = 999
)
