#######SCRIPT for generating figures2,3 manuscript###########
.libPaths("./rstudio/libs/4.4.1")
library(Seurat)
library(broom)
library(ggpubr)
library(ggplot2)
##load files##
seurat_14dpi <- readRDS("./filtered_postQC_postcb_postdoublet_postdowns_annotated_14dpi.rds")
Idents(seurat_14dpi) <- seurat_14dpi$celltypes
meta <- seurat_14dpi@meta.data
umap_df <- Embeddings(seurat_14dpi, "umap") %>%
  as.data.frame() %>%
  mutate(cluster = as.character(Idents(seurat_14dpi)))
#UMAP generation for 14 DPI#
umap_df <- Embeddings(seurat_14dpi, "umap") %>%
  as.data.frame() %>%
  mutate(cluster = as.character(Idents(seurat_14dpi)))

# Count cells per cluster
cluster_counts <- umap_df %>%
  count(cluster) %>%
  arrange(desc(n))

#Figure 2A
cluster_colors <- c(
  "CD2- GD T cells" = "plum",
  "B cells" = "orange",
  "Mixed CD8A+ AB T/NK cells" = "skyblue2",
  "CD4+CD8A- AB T cells" = "darkgreen",
  "Cytotoxic CD8A+ AB T cells" = "mediumseagreen",
  "NK cells" = "steelblue",
  "Monocytes" = "lightpink",
  "CD2+ GD T cells" = "darkblue",
  "ASCs" = "gold",
  "pDCs" = "red"
)
label_coords <- umap_df %>%
  group_by(cluster) %>%
  summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2), .groups = "drop")
umap_df$cluster <- factor(umap_df$cluster, levels = names(cluster_colors))


cluster_counts <- umap_df %>%
  count(cluster)

cluster_labels <- setNames(
  paste0(cluster_counts$cluster, " (", cluster_counts$n, " cells)"),
  cluster_counts$cluster
)

ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point(size = 0.7, alpha = 0.6) +
  scale_color_manual(
    values = cluster_colors,
    breaks = names(cluster_colors),
    labels = cluster_labels[names(cluster_colors)]
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 4),
    title.theme = element_text(face = "bold"),
    label.theme = element_text(face = "bold"),
    keywidth = 1.2, keyheight = 1.2
  )) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.box.background = element_rect(color = "black", linewidth = 0.5),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  ) +
  ggtitle(paste(length(unique(umap_df$cluster)), "celltypes")) +
  labs(color = "Celltype\n(Cell Count)", x = "UMAP 1", y = "UMAP 2")

#Figure 2B
#cell per sample
cell_counts <- meta %>%
  count(Sample, Treatment, celltypes, Sow, name = "n_cells")
#total pbmc per sample
sample_totals <- meta %>%
  count(Sample, Treatment,Sow, name = "total_cells")
#merge

cell_prop <- cell_counts %>%
  left_join(sample_totals, by = c("Sample", "Treatment","Sow")) %>%
  mutate(percent = (n_cells / total_cells) * 100)

cell_prop <- cell_prop %>%
  mutate(proportion = percent / 100)


##adding sow effect
lmm_emm <- cell_prop %>%
  group_by(celltypes) %>%
  do({
    m <- lmer(proportion ~ Treatment + (1 | Sow), data = .)
    emm <- emmeans(m, ~ Treatment)
    pairs(emm, adjust = "tukey") %>%
      as.data.frame()
  }) %>%
  ungroup()
#multiplicity across cell types - BH
lmm_emm <- lmm_emm %>%
  group_by(contrast) %>%
  mutate(p.adj.BH = p.adjust(p.value, method = "BH")) %>%
  ungroup()
lmm_emm$stars <- cut(
  lmm_emm$p.adj.BH,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("***", "**", "*", "")
)
write.csv(lmm_emm,"/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Project_Fang_10X/proportions_anova_tukey_sig_d14.csv")

plot_pbmc_simple <- ggbarplot(
  cell_prop,
  x = "celltypes",
  y = "proportion",
  fill = "Treatment",
  add = "mean_se",
  position = position_dodge(),
  palette = c("#56B4E9", "orange", "green")
  )  +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text = element_text(size = 10, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "top"
    ) +
    labs(
      x = "Cell Types",
      y = "Proportion"
    ) 
plot_pbmc_simple <- ggbarplot(
  cell_prop,
  x = "celltypes",
  y = "proportion",
  fill = "Treatment",
  add = "mean_se",
  position = position_dodge(),
  palette = c("#56B4E9", "orange", "green")
) +
  coord_flip() +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"),  # now y = celltypes
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "bottom"
  ) +
  labs(
    y = "Cell Types",
    x = "Proportion"
  )
plot_pbmc_simple

#Fig 3A
seurat_84dpi <- readRDS("./filtered_postQC_postcb_postdoublet_postdowns_postcellcycle_updated_annotation_84dpi.rds")
Idents(seurat_84dpi) <- seurat_84dpi$CellTypes
meta <- seurat_84dpi@meta.data
umap_df <- Embeddings(seurat_84dpi, "umap") %>%
  as.data.frame() %>%
  mutate(cluster = as.character(Idents(seurat_84dpi)))
# Count cells per cluster
cluster_counts <- umap_df %>%
  count(cluster) %>%
  arrange(desc(n))


cluster_colors <- c(
  "CD2- GD T cells" = "plum",
  "B cells" = "orange",
  "Transitional B-like" = "mediumseagreen",
  "AB T cells" = "gold",
  "cDCs" = "darkmagenta",
  "NK cells" = "skyblue2",
  "Monocytes" = "lightpink",
  "CD2+ GD T cells" = "steelblue",
  "ASC" = "darkgreen",
  "pDCs" = "red"
)
umap_df$cluster <- factor(umap_df$cluster, levels = names(cluster_colors))
label_coords <- umap_df %>%
  group_by(cluster) %>%
  summarize(umap_1 = mean(umap_1), umap_2 = mean(umap_2), .groups = "drop")
ggplot(umap_df, aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(size = 0.7, alpha = 0.6) +
  geom_text(
    data = label_coords,
    aes(label = cluster),
    color = "black", size = 4, fontface = "bold"
  ) +
  scale_color_manual(
    values = cluster_colors,
    labels = paste0(cluster_counts$cluster, " (", cluster_counts$n, " cells)")
  )+
  guides(color = guide_legend(
    override.aes = list(size = 4),
    title.theme = element_text(face = "bold"),
    label.theme = element_text(face = "bold"),
    keywidth = 1.2, keyheight = 1.2
  )) + 
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.box.background = element_rect(color = "black", size = 0.5),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  ) +
  ggtitle(paste(length(unique(umap_df$cluster)), "CellTypes")) +
  labs(color = "Celltype\n(Cell Count)", x = "UMAP 1", y = "UMAP 2")

#Figure 3B
#cell per sample
cell_counts <- meta %>%
  count(Sample, Treatment, CellTypes, Sow,Sex, name = "n_cells")
#total pbmc per sample
sample_totals <- meta %>%
  count(Sample, Treatment,Sow,Sex, name = "total_cells")
#merge

cell_prop <- cell_counts %>%
  left_join(sample_totals, by = c("Sample", "Treatment","Sow","Sex")) %>%
  mutate(proportion = (n_cells / total_cells) )


##adding sow effect
lmm_emm <- cell_prop %>%
  group_by(CellTypes) %>%
  do({
    m <- lmer(proportion ~ Treatment + Sex+ (1 | Sow), data = .)
    emm <- emmeans(m, ~ Treatment)
    pairs(emm, adjust = "tukey") %>%
      as.data.frame()
  }) %>%
  ungroup()
lmm_emm <- lmm_emm %>%
  group_by(contrast) %>%
  mutate(p.adj.BH = p.adjust(p.value, method = "BH")) %>%
  ungroup()
lmm_emm$stars <- cut(
  lmm_emm$p.value,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("***", "**", "*", "")
)
plot_pbmc_simple <- ggbarplot(
  cell_prop,
  x = "CellTypes",
  y = "proportion",
  fill = "Treatment",
  add = "mean_se",
  position = position_dodge(),
  palette = c("#56B4E9", "orange", "green")
)  + coord_flip() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "top"
  ) +
  labs(
    x = "Cell Types",
    y = "Proportion"
  )

plot_pbmc_simple  #add significant red stars on transitional b like(EC) and CD2+ GD (EC,PC), CD2- GD(EP)
p_zoom <- plot_pbmc_simple + coord_cartesian(ylim = c(0, 0.3))

plot_pbmc_simple <- ggboxplot(
  cell_prop,
  x = "CellTypes",
  y = "proportion",
  fill = "Treatment",
  #add = "mean_se",
  palette = c("#56B4E9", "orange", "green")
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "top"
  ) +
  labs(
    x = "Cell Types",
    y = "Proportion"
  )
plot_pbmc_simple

##AB proportions

seurat_84dpi <- readRDS("./filtered_postQC_postcb_postdoublet_postdowns_postcellcycle_ABT_subset_84dpi.rds")
Idents(seurat_84dpi) <- seurat_84dpi$Treatment
meta <- seurat_84dpi@meta.data
umap_df <- Embeddings(seurat_84dpi, "umap") %>%
  as.data.frame() %>%
  mutate(cluster = as.character(Idents(seurat_84dpi)))
# Count cells per cluster
cluster_counts <- umap_df %>%
  count(cluster) %>%
  arrange(desc(n))


cluster_colors <- c(
  "CD2- GD T cells" = "plum",
  "Cytotoxic CD8a+ T cell" = "mediumseagreen",
  "Naïve CD4+CD8a- ab T cells" = "darkmagenta",
  "Mixed CD4+CD8a+ ab T" = "steelblue",
  "Activated CD4+CD8a+ ab T cells" = "red"
)
cluster_colors <- c(
  "control" ="#56B4E9", 
  "extinct"="orange",
  "persistent"= "green"
)
DimPlot(seurat_84dpi, split.by = "Treatment", cols = cluster_colors)
umap_df$cluster <- factor(umap_df$cluster, levels = names(cluster_colors))
label_coords <- umap_df %>%
  group_by(cluster) %>%
  summarize(umap_1 = mean(umap_1), umap_2 = mean(umap_2), .groups = "drop")
ggplot(umap_df, aes(x = umap_1, y = umap_2, color = cluster)) +
  geom_point(size = 0.7, alpha = 0.6) +
  geom_text(
    data = label_coords,
    aes(label = cluster),
    color = "black", size = 4, fontface = "bold"
  ) +
  scale_color_manual(
    values = cluster_colors,
    labels = paste0(cluster_counts$cluster, " (", cluster_counts$n, " cells)")
  )+
  guides(color = guide_legend(
    override.aes = list(size = 4),
    title.theme = element_text(face = "bold"),
    label.theme = element_text(face = "bold"),
    keywidth = 1.2, keyheight = 1.2
  )) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.box.background = element_rect(color = "black", size = 0.5),
    axis.line = element_line(color = "black"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  ) +
  ggtitle(paste(length(unique(umap_df$cluster)), "CellTypes")) +
  labs(color = "Celltype\n(Cell Count)", x = "UMAP 1", y = "UMAP 2")

#Figure 3B
#cell per sample
cell_counts <- meta %>%
  count(Sample, Treatment, CellTypes, Sow,Sex, name = "n_cells")
#total pbmc per sample
sample_totals <- meta %>%
  count(Sample, Treatment,Sow,Sex, name = "total_cells")
#merge

cell_prop <- cell_counts %>%
  left_join(sample_totals, by = c("Sample", "Treatment","Sow","Sex")) %>%
  mutate(proportion = (n_cells / total_cells) )


##adding sow effect
lmm_emm <- cell_prop %>%
  group_by(CellTypes) %>%
  do({
    m <- lmer(proportion ~ Treatment + Sex+ (1 | Sow), data = .)
    emm <- emmeans(m, ~ Treatment)
    pairs(emm, adjust = "tukey") %>%
      as.data.frame()
  }) %>%
  ungroup()
lmm_emm <- lmm_emm %>%
  group_by(contrast) %>%
  mutate(p.adj.BH = p.adjust(p.value, method = "BH")) %>%
  ungroup()
lmm_emm$stars <- cut(
  lmm_emm$p.value,
  breaks = c(-Inf, 0.01, 0.05, 0.1, Inf),
  labels = c("***", "**", "*", "")
)
plot_pbmc_simple <- ggbarplot(
  cell_prop,
  x = "CellTypes",
  y = "proportion",
  fill = "Treatment",
  add = "mean_se",
  position = position_dodge(),
  palette = c("#56B4E9", "orange", "green")
)  +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "top"
  ) +
  labs(
    x = "Cell Types",
    y = "Proportion"
  )
