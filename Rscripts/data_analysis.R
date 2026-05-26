#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                              #
#                                                                              #
#                                                                              #
#                                                                              #
# Libraries
require("ggplot2")  
require("BiocManager")
require("dplyr")
require("vegan")
BiocManager::install("phyloseq")
require("phyloseq")
require("tidyverse")


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

# Bray-Curtis NMDS
ord_nmds <- ordinate(ps_rare, method = "NMDS", distance = "bray")

# Stress for bray-curtis distance, less than 0.2
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

permanova <- adonis2(bray_dist ~ Treatment,
                     data = as.data.frame(metadata),
                     permutations = 999)
permanova #Not significant when evaluating everything at once

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

adonis2(bray_uni ~ Treatment,
        data = data.frame(sample_data(ps_uni)),
        permutations = 999)

graph_nmds_uni <- ggplot(nmds_points_uni, aes(x = MDS1, y = MDS2, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +  
  stat_ellipse(aes(group = Treatment), 
               type = "t",level = 0.95)+
  labs(title = "NMDS - Bray Curtis",
       subtitle = paste("Stress =", round(ord_nmds_uni$stress, 3))) +
  theme_bw()

#Only pesticides

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

#Only flavonoids

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
