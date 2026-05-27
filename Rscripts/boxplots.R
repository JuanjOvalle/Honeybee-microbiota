# Phylum range
ps_phylum <- tax_glom(ps_uni, taxrank = "Phylum")
  
# Mutation
phylum_df <- psmelt(ps_phylum) %>%
  group_by(Sample) %>%
  mutate(rel_abund = Abundance / sum(Abundance) * 100) %>%
  ungroup()

# Filter by abundance, greater than 3%
phylum_pool <- phylum_df %>%
  group_by(Treatment, Phylum) %>%
  summarize(mean = mean(rel_abund), .groups = "drop") %>%
  group_by(Phylum) %>%
  summarize(pool = max(mean) < 3,
            mean = mean(mean),
            .groups = "drop")

# Plot
phylum_df %>%
  left_join(phylum_pool, by = "Phylum") %>%
  mutate(Phylum = if_else(pool, "Other", Phylum)) %>%
  group_by(Sample, Treatment, Phylum) %>%
  summarize(rel_abund = sum(rel_abund),
            mean = min(mean),
            .groups = "drop") %>%
  mutate(Phylum = fct_reorder(Phylum, mean, .desc = TRUE)) %>%
  ggplot(aes(x = Phylum, y = rel_abund, fill = Treatment)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.4,
                                              dodge.width = 0.8),
              pch = 21, stroke = 0, size = 1.8) +
  stat_summary(fun = median, geom = "crossbar",
               position = position_dodge(width = 0.8),
               width = 0.8, size = 0.25,
               show.legend = FALSE) +
  labs(x = NULL,
       y = "Relative Abundance (%)") +
  theme_classic()+
  facet_wrap(~ Treatment, nrow = 2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
