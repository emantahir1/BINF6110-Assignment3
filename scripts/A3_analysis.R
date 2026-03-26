# BINF6110 Assignment 3 - Shotgun Metagenomics Analysis
# Eman Tahir
# Vegan vs Omnivore gut microbiome comparison

library(phyloseq)
library(biomformat)
library(ggplot2)
library(vegan)
library(DESeq2)
library(reshape2)
library(patchwork)

# ── 1. LOAD DATA ──────────────────────────────────────────────────────────────

# Load the BIOM file
biom_data <- read_biom("~/Downloads/A3_metagenomics/results/table.biom")
physeq <- import_biom(biom_data)

# Fix taxonomy column names
colnames(tax_table(physeq)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Add metadata with correct sample names
metadata <- data.frame(
  Diet = c("Vegan", "Omnivore", "Omnivore", "Omnivore", "Vegan", "Vegan"),
  row.names = c("SRR8146944_bracken", "SRR8146970_bracken", "SRR8146972_bracken",
                "SRR8146976_bracken", "SRR8146983_bracken", "SRR8146985_bracken")
)
sample_data(physeq) <- metadata

# Rename samples
sample_names(physeq) <- c("Vegan1", "Omnivore1", "Omnivore2", "Omnivore3", "Vegan2", "Vegan3")

# Convert to relative abundance
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

# ── 2. FIGURE 1: TAXONOMIC ABUNDANCE ─────────────────────────────────────────

physeq_phylum <- tax_glom(physeq_rel, taxrank = "Phylum")
df <- psmelt(physeq_phylum)
df$Phylum <- gsub("p__", "", df$Phylum)

ggplot(df, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~Diet, scales = "free_x") +
  labs(title = "Taxonomic Abundance by Diet Group",
       y = "Relative Abundance", x = "Sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("~/Desktop/figure1_taxonomy.png", width = 10, height = 6, dpi = 300)

# ── 3. FIGURE 2: ALPHA DIVERSITY ─────────────────────────────────────────────

alpha_div <- estimate_richness(physeq, measures = c("Observed", "Shannon", "Simpson"))
alpha_div$Diet <- sample_data(physeq)$Diet
alpha_div$Sample <- rownames(alpha_div)

# Wilcoxon tests
wilcox.test(Observed ~ Diet, data = alpha_div)
wilcox.test(Shannon ~ Diet, data = alpha_div)
wilcox.test(Simpson ~ Diet, data = alpha_div)

alpha_long <- reshape2::melt(alpha_div, id.vars = c("Sample", "Diet"),
                             variable.name = "Measure", value.name = "Value")

ggplot(alpha_long, aes(x = Diet, y = Value, color = Diet)) +
  geom_point(size = 5, alpha = 0.8) +
  facet_wrap(~Measure, scales = "free_y") +
  scale_color_manual(values = c("Omnivore" = "#E41A1C", "Vegan" = "#377EB8")) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.3,
               color = "black", linewidth = 0.5) +
  labs(title = "Alpha Diversity Measures by Diet Group",
       x = "Diet", y = "Value") +
  theme_bw() +
  theme(legend.position = "none")

ggsave("~/Desktop/figure2_alpha.png", width = 10, height = 5, dpi = 300)

# ── 4. FIGURE 3: BETA DIVERSITY ──────────────────────────────────────────────

# PCoA
ord.pcoa.bray <- ordinate(physeq, method = "PCoA", distance = "bray")
pcoa_df <- data.frame(ord.pcoa.bray$vectors[, 1:2])
colnames(pcoa_df) <- c("Axis1", "Axis2")
pcoa_df$Diet <- sample_data(physeq)$Diet
pcoa_df$Sample <- rownames(pcoa_df)
eig <- ord.pcoa.bray$values$Eigenvalues
var_exp <- round(eig / sum(eig) * 100, 1)

# PERMANOVA
metadata_df <- as(sample_data(physeq), "data.frame")
adonis2(phyloseq::distance(physeq, method = "bray") ~ Diet, data = metadata_df)

# NMDS
ord.nmds.bray <- ordinate(physeq, method = "NMDS", distance = "bray")
nmds_df <- data.frame(ord.nmds.bray$points)
colnames(nmds_df) <- c("NMDS1", "NMDS2")
nmds_df$Diet <- sample_data(physeq)$Diet
nmds_df$Sample <- rownames(nmds_df)

p_pcoa <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, color = Diet, label = Sample)) +
  geom_point(size = 5) +
  geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("Omnivore" = "#E41A1C", "Vegan" = "#377EB8")) +
  labs(title = "PCoA (Bray-Curtis)",
       x = paste0("PCoA1 [", var_exp[1], "%]"),
       y = paste0("PCoA2 [", var_exp[2], "%]")) +
  theme_bw() +
  expand_limits(x = c(-0.6, 0.6))

p_nmds <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Diet, label = Sample)) +
  geom_point(size = 5) +
  geom_text(vjust = -0.8, size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("Omnivore" = "#E41A1C", "Vegan" = "#377EB8")) +
  labs(title = "NMDS (Bray-Curtis)",
       subtitle = paste0("Stress = ", round(ord.nmds.bray$stress, 3))) +
  theme_bw() +
  expand_limits(x = c(-0.6, 0.6))

p_pcoa + p_nmds + plot_layout(guides = "collect")
ggsave("~/Desktop/figure3_beta_combined.png", width = 14, height = 6, dpi = 300)

# ── 5. FIGURE 4: DIFFERENTIAL ABUNDANCE ──────────────────────────────────────

diagdds <- phyloseq_to_deseq2(physeq, ~ Diet)
diagdds <- DESeq(diagdds, test = "Wald", fitType = "local")
res <- results(diagdds, contrast = c("Diet", "Vegan", "Omnivore"))
res <- res[order(res$padj, na.last = NA), ]
res_df <- as.data.frame(res)
res_df$taxon <- rownames(res_df)
tax <- as.data.frame(tax_table(physeq))
res_df$Species <- tax[res_df$taxon, "Species"]
res_df$Genus <- tax[res_df$taxon, "Genus"]
res_df$label <- paste(gsub("g__", "", res_df$Genus), gsub("s__", "", res_df$Species))
res_df <- res_df[!is.na(res_df$pvalue), ]
top_taxa <- head(res_df[order(res_df$pvalue), ], 20)

ggplot(top_taxa, aes(x = log2FoldChange, y = reorder(label, log2FoldChange),
                     fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("TRUE" = "#377EB8", "FALSE" = "#E41A1C"),
                    labels = c("TRUE" = "Vegan", "FALSE" = "Omnivore"),
                    name = "Higher in") +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  labs(title = "Differential Abundance: Vegan vs Omnivore",
       subtitle = "Top 20 taxa by p-value (log2 fold change)",
       x = "Log2 Fold Change (Vegan / Omnivore)",
       y = "Taxon") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9))

ggsave("~/Desktop/figure4_diffabund.png", width = 10, height = 7, dpi = 300)