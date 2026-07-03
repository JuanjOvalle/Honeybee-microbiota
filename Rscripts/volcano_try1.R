#ggrepel ya instalado en el coso base, importante

deseq_all <- phyloseq_to_deseq2(ps, ~Treatment)

cts <- counts(deseq_all)
geoMeans <- apply(cts, 1, function(row) {
  if (all(row == 0)) return(0)
  exp(sum(log(row[row > 0])) / length(row))
})

deseq_all <- estimateSizeFactors(deseq_all, geoMeans = geoMeans)

deseq_all <- DESeq(deseq_all, fitType = "local")

res <- results(deseq_all, contrast = c("Treatment", "Rutina", "Glifosato"))
res <- res[order(res$padj), ]
summary(res)
head(res)


res_df <- as.data.frame(res)
res_df$p.value <- ifelse(!is.na(res_df$padj) & res_df$padj < 0.05, "<0.05", ">0.05")

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = p.value)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("<0.05" = "firebrick", ">0.05" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = "Volcano plot - Rutina vs Control_(N)",
       x = "log2 Fold Change", y = "-log10(padj)")


# Lista de todos los tratamientos
tratamientos <- unique(data.frame(sample_data(ps))$Treatment)

# Todas las combinaciones posibles contra el control
# O si querés todos los pares entre sí:
pares <- combn(tratamientos, 2, simplify = FALSE)
res_df <- as.data.frame(res)
res_df <- cbind(res_df, as(tax_table(ps)[rownames(res_df), ], "matrix"))
res_df$Label <- ifelse(res_df$p.value == "<0.05", as.character(res_df$Genus), NA)


# Función para correr DESeq2 y generar volcano plot para cada par
volcano_par <- function(trat1, trat2, deseq_obj) {
  
  res <- results(deseq_obj, 
                 contrast = c("Treatment", trat1, trat2))
  
  res_df <- as.data.frame(res) %>%
    mutate(p.value = ifelse(!is.na(padj) & padj < 0.05, "<0.05", ">0.05"))
  res_df <- cbind(res_df, as(tax_table(ps)[rownames(res), ], "matrix"))
  res_df$Label <- ifelse(res_df$p.value == "<0.05", as.character(res_df$Genus), NA)
  
  p <- ggplot(res_df, aes(x = log2FoldChange, 
                          y = -log10(padj), 
                          color = p.value)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c("<0.05" = "firebrick", ">0.05" = "grey70")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_text_repel(aes(label = Label), 
                    size = 3, 
                    max.overlaps = 20,     
                    na.rm = TRUE,
                    show.legend = FALSE) +
    theme_minimal() +
    labs(title = paste("Volcano plot -", trat1, "vs", trat2),
         x = "log2 Fold Change", 
         y = "-log10(padj)")
  
  # Guardar figura
  ggsave(paste0("figures/volcano_", trat1, "_vs_", trat2, ".jpg"), 
         plot = p, width = 6, height = 4)
  
  return(p)
}

# Correr para todos los pares
plots <- lapply(pares, function(par) {
  volcano_par(par[1], par[2], deseq_all)
})