# install.packages("ggrepel")  # si no lo tienes instalado
library(ggrepel)

# Necesitas que res_df tenga la taxonomía pegada, igual que hiciste para res_tax
res_df <- as.data.frame(res)
res_df <- cbind(res_df, as(tax_table(ps)[rownames(res_df), ], "matrix"))
res_df$Significativo <- ifelse(!is.na(res_df$padj) & res_df$padj < 0.05, "Sí", "No")

# Elige qué nivel taxonómico mostrar como etiqueta (ej. Genus)
res_df$Label <- ifelse(res_df$Significativo == "Sí", as.character(res_df$Genus), NA)

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Significativo)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Sí" = "firebrick", "No" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_text_repel(aes(label = Label), 
                  size = 3, 
                  max.overlaps = 20,     # aumenta si se omiten etiquetas
                  na.rm = TRUE,
                  show.legend = FALSE) +
  theme_minimal() +
  labs(title = "Volcano plot - Rutina vs Control_(N)",
       x = "log2 Fold Change", y = "-log10(padj)")
