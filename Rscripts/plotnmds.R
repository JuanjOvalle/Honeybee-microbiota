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
nochim1 <- read.csv("C:/Users/juanj/Documents/Abejas-Microbiota/Abejas-Microbiota/Analisis/seqtab_nochim_transp_1.csv", row.names = 1, header = TRUE)
taxa <- read.csv("C:/Users/juanj/Documents/Abejas-Microbiota/Abejas-Microbiota/Analisis/taxa_bee_16s.csv", row.names = 1)
metadata <- read.csv("C:/Users/juanj/Documents/Abejas-Microbiota/Abejas-Microbiota/Analisis/metadata_bee_16s.csv", row.names = 1)

# Matrixes
nochim_mat <- as.matrix(nochim1)
taxa_mat <- as.matrix(taxa)

# Phyloseq database
ps <- phyloseq(otu_table(nochim_mat, taxa_are_rows = TRUE),
  tax_table(taxa_mat)
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

# Bray-Curtis
ord_nmds <- ordinate(ps_rare, method = "NMDS", distance = "bray")#NMDS
ord_pcoa <- ordinate(ps_rare, method = "PCoA", distance = "bray")#PCoA

# Stress for bray-curtis distance, less than 0.2
ord_nmds$stress

# Graphs
plot_ordination(ps_rare, ord_nmds, 
                type = "samples", 
                color = "Treatment") +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "NMDS - Bray Curtis",
       subtitle = paste("Stress =", round(ord_nmds$stress, 3))) +
  theme_bw()

plot_ordination(ps_rare, ord_pcoa, 
                type = "samples", 
                color = "Treatment") +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = rainbow(20)) +
  theme_bw()

plot_ordination(ps_rare, ord_pcoa, type = "samples") +
  geom_point(size = 3, alpha = 0.7) +
  labs(title = "PCoA - Bray Curtis") +
  theme_bw()


# Agregar metadata al ps original
sample_data(ps) <- sample_data(metadata)

# Rarefy DESPUÉS de agregar metadata
ps_rare <- rarefy_even_depth(ps, rngseed = 123)

# Verificar que ps_rare tiene la metadata
sample_data(ps_rare)

# Correr ordenación de nuevo
ord_nmds <- ordinate(ps_rare, method = "NMDS", distance = "bray")

# Ver si la metadata está en ps_rare
head(sample_data(ps_rare))

# Ver los nombres de las columnas
colnames(sample_data(ps_rare))

nmds_points <- as.data.frame(ord_nmds$points)
nmds_points$Treatment <- sample_data(ps_rare)$Treatment

# Graficar
graph_nmds_all <- ggplot(nmds_points, aes(x = MDS1, y = MDS2, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +  
  stat_ellipse(aes(group = Treatment), 
        type = "t",level = 0.95)+
  labs(title = "NMDS - Bray Curtis",
       subtitle = paste("Stress =", round(ord_nmds$stress, 3))) +
  theme_bw()

ggsave("Figures/graph_nmds_all.jpg", plot = graph_nmds_all)

pcoa_points <- as.data.frame(ord_pcoa$vectors)
pcoa_points$Treatment <- sample_data(ps_rare)$Treatment

ggplot(pcoa_points, aes(x = Axis.1, y = Axis.2, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +
  stat_ellipse(aes(group = Treatment),
               type = "t",level = 0.95)+
  theme_bw()

bray_dist <- phyloseq::distance(ps_rare, method = "bray")

permanova <- adonis2(bray_dist ~ Treatment,
                     data = as.data.frame(metadata),
                     permutations = 999)
permanova

sample_data(ps_rare)$Treated <- ifelse(
  sample_data(ps_rare)$Treatment == "Control", 
  "Control", 
  "Treated"
)

bray_dist <- phyloseq::distance(ps_rare, method = "bray")

adonis2(bray_dist ~ Treatment,
        data = as.data.frame(metadata),
        permutations = 999)

nrow(metadata)
nsamples(ps_rare)

ps_uni <- subset_samples(ps_rare, 
                         Treatment %in% c("Control_(N)", "Rutina", "Imidacloprid", "Kaempferol", 
                                          "p-Coumaric_acid", "Fipronil", "Sugar", "Deltametrina", "Glifosato"))

ord_nmds_uni <- ordinate(ps_uni, method = "NMDS", distance = "bray")
bray_uni <- phyloseq::distance(ps_uni, method = "bray")

#Pesticides vs pesticides + flavonoids

adonis2(bray_uni ~ Treatment,
        data = data.frame(sample_data(ps_uni)),
        permutations = 999)

ps_pesticides <- subset_samples(ps_rare,
                           Treatment %in% c("Imidacloprid", "Fipronil", "Deltametrina", "Glifosato", "Kaempferol+Deltametrina", "Kaempferol+Fipronil", "Kaempferol+Glifosato", "Kaempferol+Imidaclorid", "	
p-Coumaric_acid+Deltametrina", "p-Coumaric_acid+Fipronil", "p-Coumaric_acid+Glifosato", "p-Coumaric_acid+Imidaclorid", "Rutina+Deltametrina", "Rutina+Fipronil", "	
Rutina+Glifosato", "Rutina+Imidaclorid")) 

bray_pesticides <- phyloseq::distance(ps_pesticides, method = "bray")

adonis2(bray_pesticides ~ Treatment,
        data = data.frame(sample_data(ps_pesticides)),
        permutations = 999)

#Per molecule analisis

ps_monomolecule <- function()

ps_mol1 <- subset_samples(ps_rare,
                          Treatment %in% c("Control", "mol1"))

bray_mol1 <- phyloseq::distance(ps_mol1, method = "bray")

adonis2(bray_mol1 ~ Treatment,
        data = data.frame(sample_data(ps_mol1)),
        permutations = 999)

#Graph monomolecules

nmds_points_mono <- as.data.frame(ord_nmds_uni$points)
nmds_points_mono$Treatment <- sample_data(ps_uni)$Treatment

# Graficar
graph_nmds_uni <- ggplot(nmds_points_mono, aes(x = MDS1, y = MDS2, color = Treatment)) +
  geom_point(size = 3, alpha = 0.7) +  
  stat_ellipse(aes(group = Treatment), 
               type = "t",level = 0.95)+
  labs(title = "NMDS - Bray Curtis",
       subtitle = paste("Stress =", round(ord_nmds_uni$stress, 3))) +
  theme_bw()
