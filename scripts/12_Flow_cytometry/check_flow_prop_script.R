
## Flow cytometry for day 84 ###
.libPaths("/work/ABG/mkapoor/mkapoor/.ondemand-new/mkapoor/rstudio/libs/4.4.1")

seurat_84dpi <- readRDS("/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/filtered_postQC_postcb_postdoublet_postdowns_postcellcycle_updated_annotation_84dpi.rds")
#seurat_84dpi_subset_1 <- readRDS("/work/ABG/mkapoor/PRRSV/PRRSV_cellranger_v97/filtered_postQC_postcb_postdoublet_postdowns_postcellcycle_ABT_subset_84dpi.rds")
#number of cells ineach celltype in each treatment
table(seurat_84dpi$CellTypes, seurat_84dpi$Treatment)
#table(seurat_84dpi_subset_1$CellTypes, seurat_84dpi_subset_1$Treatment)
Idents(seurat_84dpi) <- "CellTypes" 
meta <- seurat_84dpi@meta.data
# Calculate total number of cells per sample
cell_type_counts <- meta %>%
  group_by(Sample) %>%
  mutate(total_cells = n()) %>% # Now count frequency of each celltype per sample
  group_by(CellTypes, Sample, Treatment) %>%
  mutate(frequency = n() / unique(total_cells)) %>%
  ungroup()
cell_type_summary <- cell_type_counts %>%
  dplyr::select(CellTypes, Sample, Treatment, frequency) %>%
  distinct()
#create heatmap
heatmap_df <- cell_type_summary %>%
  dplyr::select(Sample, CellTypes, frequency) %>%
  pivot_wider(names_from = CellTypes, values_from = frequency, values_fill = 0)

# Convert to long format for ggplot2
heatmap_long <- heatmap_df %>%
  pivot_longer(-Sample, names_to = "CellTypes", values_to = "frequency")

# Plot heatmap using ggplot
ggplot(heatmap_long, aes(x = CellTypes, y = Sample, fill = frequency)) +
  geom_tile(color = "gray80", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.5f", frequency)), size = 3, fontface = "bold") +
  scale_fill_gradientn(
    colors = c("blue", "white", "red"),
    name = "Frequency",
    limits = c(0, 0.5),
    oob = scales::squish
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  ) +
  labs(
    title = "Cell Type Frequency per Sample - High Resolution",
    x = "Cell Type", y = "Sample"
  )

#some distributions - Sex related#
sex_cluster_df <- as.data.frame(seurat_84dpi@meta.data)
sex_treat_df <- as.data.frame(seurat_84dpi@meta.data)
# Build contingency table: Cluster × Sex
sex_counts <- as.data.frame(table(Cluster = sex_cluster_df$louvain_res0_8,
                                  Sex = sex_cluster_df$Sex))

# Plot
ggplot(sex_counts, aes(x = as.factor(Cluster), y = Freq, fill = Sex)) +
  geom_bar(stat = "identity", position = "stack") +
  ylab("Count of Cells") + xlab("Louvain Cluster (res=0.25)") +
  ggtitle("Sex Distribution Across Louvain Clusters") +
  scale_fill_manual(values = c("Male" = "#1f77b4", "Female" = "#ff7f0e")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Build contingency table: Treatment × Sex
sex_counts_t <- as.data.frame(table(Treatment = sex_treat_df$Treatment,
                                    Sex = sex_treat_df$Sex))

# Plot
ggplot(sex_counts_t, aes(x = as.factor(Treatment), y = Freq, fill = Sex)) +
  geom_bar(stat = "identity", position = "stack") +
  ylab("Count of Cells") + xlab("Treatment") +
  ggtitle("Sex Distribution Across Treatment") +
  scale_fill_manual(values = c("Male" = "#1f77b4", "Female" = "#ff7f0e")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

###sow related ####
sow_cluster_df <- as.data.frame(seurat_84dpi@meta.data)
sow_treat_df <- as.data.frame(seurat_84dpi@meta.data)
# Build contingency table: Cluster × Sex
sow_counts <- as.data.frame(table(Cluster = sow_cluster_df$louvain_res0_8,
                                  Sow = sow_cluster_df$Sow))

# Plot
ggplot(sow_counts, aes(x = as.factor(Cluster), y = Freq, fill = Sow)) +
  geom_bar(stat = "identity", position = "fill") +
  ylab("Count of Cells") + xlab("Louvain Cluster (res=0.25)") +
  ggtitle("Sow Distribution Across Louvain Clusters") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Build contingency table: Treatment × Sow
sow_counts_t <- as.data.frame(table(Treatment = sow_treat_df$Treatment,
                                    Sow = sow_treat_df$Sow))

# Plot
ggplot(sow_counts_t, aes(x = as.factor(Treatment), y = Freq, fill = Sow)) +
  geom_bar(stat = "identity", position = "stack") +
  ylab("Count of Cells") + xlab("Treatment") +
  ggtitle("Sow Distribution Across Treatment")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

##get umap with cell numbers
Idents(seurat_84dpi) <- seurat_84dpi$CellTypes

umap_df <- Embeddings(seurat_84dpi, "umap") %>%
  as.data.frame() %>%
  mutate(cluster = as.character(Idents(seurat_84dpi)))
# Count cells per cluster
cluster_counts <- umap_df %>%
  count(cluster) %>%
  arrange(desc(n))


#cluster_colors <- setNames(
#  hue_pal()(length(unique(umap_df$cluster))),
#  sort(unique(umap_df$cluster))
#)
#cols <- c('lightpink', 'red', 'orange', 'gold', 'darkgreen', 'mediumseagreen', 'skyblue2', 'steelblue', 'navy', 'plum3', 'darkmagenta', 'black', 'grey')
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

# Get total counts across all treatments
Idents(seurat_84dpi) <- seurat_84dpi$Treatment
celltype_time_prop<-prop.table(table(Idents(seurat_84dpi), seurat_84dpi$"CellTypes"), margin = 2)
celltype_time_prop<- as.data.frame(celltype_time_prop)
colnames(celltype_time_prop)[1] <- "Treatment"
colnames(celltype_time_prop)[2] <- "CellTypes"
colnames(celltype_time_prop)[3] <- "Proportion"
colors <- c("#56B4E9", "orange", "green")
Idents(seurat_84dpi) <- seurat_84dpi$CellTypes
cell_counts <- table(Idents(seurat_84dpi)) %>% sort(decreasing = TRUE)
celltype_order <- names(cell_counts)
celltype_time_prop$CellTypes <- factor(celltype_time_prop$CellTypes, levels = celltype_order)
plot3 <- ggplot(celltype_time_prop, aes(x = CellTypes, y = Proportion, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Proportion of cells within CellTypes by Treatment",
       x = "CellTypes",
       y = "Proportion") +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 12, face = "bold"),
        legend.position = "top") +
  scale_fill_manual(values = colors)

plot3

# Flow cytometry for day 14 ###

seurat_14dpi <- readRDS("/work/abg/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Project_Fang_10X/filtered_postQC_postcb_postdoublet_postdowns_annotated_14dpi.rds")

Idents(seurat_14dpi) <- seurat_14dpi$celltypes
meta <- seurat_14dpi@meta.data
# Calculate total number of cells per sample
cell_type_counts <- meta %>%
  group_by(Sample) %>%
  mutate(total_cells = n()) %>% # Now count frequency of each celltype per sample
  group_by(celltypes, Sample, Treatment) %>%
  mutate(frequency = n() / unique(total_cells)) %>%
  ungroup()
cell_type_summary <- cell_type_counts %>%
  dplyr::select(celltypes, Sample, Treatment, frequency) %>%
  distinct()
#create heatmap
heatmap_df <- cell_type_summary %>%
  dplyr::select(Sample, celltypes, frequency) %>%
  pivot_wider(names_from = celltypes, values_from = frequency, values_fill = 0)

# Convert to long format for ggplot2
heatmap_long <- heatmap_df %>%
  pivot_longer(-Sample, names_to = "celltypes", values_to = "frequency")

# Plot heatmap using ggplot
ggplot(heatmap_long, aes(x = celltypes, y = Sample, fill = frequency)) +
  geom_tile(color = "gray80", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.5f", frequency)), size = 3, fontface = "bold") +
  scale_fill_gradientn(
    colors = c("blue", "white", "red"),
    name = "Frequency",
    limits = c(0, 0.5),
    oob = scales::squish
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    legend.position = "right"
  ) +
  labs(
    title = "Cell Type Frequency per Sample - High Resolution",
    x = "Cell Type", y = "Sample"
  )

###sow related ####
sow_cluster_df <-data.frame(Sow = meta$Sow, louvain_res1_5 = meta$louvain_res1_5)
sow_treat_df <- data.frame(Sow = meta$Sow, Treatment = meta$Treatment)
# Build contingency table: Cluster × Sow
sow_counts <- as.data.frame(table(louvain_res1_5 = sow_cluster_df$louvain_res1_5,
                                  Sow = sow_cluster_df$Sow))

##get umap with cell numbers
Idents(seurat_14dpi) <- seurat_14dpi$celltypes

umap_df <- Embeddings(seurat_14dpi, "umap") %>%
  as.data.frame() %>%
  mutate(cluster = as.character(Idents(seurat_14dpi)))

# Count cells per cluster
cluster_counts <- umap_df %>%
  count(cluster) %>%
  arrange(desc(n))


#cluster_colors <- setNames(
#  hue_pal()(length(unique(umap_df$cluster))),
#  sort(unique(umap_df$cluster))
#)
#cols <- c('lightpink', 'red', 'orange', 'gold', 'darkgreen', 'mediumseagreen', 'skyblue2', 'steelblue', 'navy', 'plum3', 'darkmagenta', 'black', 'grey')
cluster_colors <- c(
  "CD2- GD T cells" = "plum",
  "B cells" = "orange",
  "Innate CD8A+ ab T/NK cells" = "skyblue2",
  "Naïve CD4+CD8A- ab T cells" = "darkgreen",
  "Cytotoxic CD8A+ ab T cells" = "mediumseagreen",
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
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point(size = 0.7, alpha = 0.6) +
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
  ggtitle(paste(length(unique(umap_df$cluster)), "celltypes")) +
  labs(color = "Celltype\n(Cell Count)", x = "UMAP 1", y = "UMAP 2")

DimPlot(seurat_14dpi, group.by = "louvain_res1_5", label =TRUE)
#### cluster
###############
#get umap with cell numbers
Idents(seurat_14dpi) <- seurat_14dpi$louvain_res1_5

umap_df <- Embeddings(seurat_14dpi, "umap") %>%
  as.data.frame() %>%
  mutate(cluster = as.character(Idents(seurat_14dpi)))

# Count cells per cluster
cluster_counts <- umap_df %>%
  count(cluster) %>%
  arrange(desc(n))


cluster_colors <- setNames(
 hue_pal()(length(unique(umap_df$cluster))),
 sort(unique(umap_df$cluster))
)

label_coords <- umap_df %>%
  group_by(cluster) %>%
  summarize(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2), .groups = "drop")

ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
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
  ggtitle(paste(length(unique(umap_df$cluster)), "celltypes")) +
  labs(color = "Cluster\n(Cell Count)", x = "UMAP 1", y = "UMAP 2")

#first lets plot some stacked barplots to see what cluster is driven by treatments
Idents(seurat_14dpi) <- seurat_14dpi$louvain_res1_5
Idents(seurat_14dpi) <- seurat_14dpi$Treatment
celltype_time_prop<-prop.table(table(Idents(seurat_14dpi), seurat_14dpi$"louvain_res1_5"), margin = 2)
celltype_time_prop<- as.data.frame(celltype_time_prop)
colnames(celltype_time_prop)[1] <- "Treatment"
colnames(celltype_time_prop)[2] <- "cluster"
colnames(celltype_time_prop)[3] <- "Proportion"
colors <- c("dodgerblue", "orange", "green")
plot1 <- ggplot(celltype_time_prop, aes(x = cluster, y = Proportion, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +  # Added position = "dodge" for unstacked bars
  theme_minimal() +
  labs(title = "Proportion of cells within cluster by treatement",
       x = "Cluster",
       y = "Proportion") +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 12, face = "bold"),
        legend.position = "top") +
  scale_fill_manual(values = colors)
plot1

#now lets compare on celltypoe level
Idents(seurat_14dpi) <- seurat_14dpi$Treatment
celltype_time_prop<-prop.table(table(Idents(seurat_14dpi), seurat_14dpi$"celltypes"), margin = 2)
celltype_time_prop<- as.data.frame(celltype_time_prop)
colnames(celltype_time_prop)[1] <- "Treatment"
colnames(celltype_time_prop)[2] <- "celltypes"
colnames(celltype_time_prop)[3] <- "Proportion"
colors <- c("#56B4E9", "orange", "green")
plot2 <- ggplot(celltype_time_prop, aes(x = celltypes, y = Proportion, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +  # Added position = "dodge" for unstacked bars
  theme_minimal() +
  labs(title = "Proportion of cells within CellTypes by treatement",
       x = "CellTypes",
       y = "Proportion") +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 12, face = "bold"),
        legend.position = "top") +
  scale_fill_manual(values = colors)
plot2

# Get total counts across all treatments
Idents(seurat_14dpi) <- seurat_14dpi$celltypes
cell_counts <- table(Idents(seurat_14dpi)) %>% sort(decreasing = TRUE)
celltype_order <- names(cell_counts)
celltype_time_prop$celltypes <- factor(celltype_time_prop$celltypes, levels = celltype_order)
plot3 <- ggplot(celltype_time_prop, aes(x = celltypes, y = Proportion, fill = Treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Proportion of cells within CellTypes by Treatment",
       x = "CellTypes",
       y = "Proportion") +
  theme(axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 12, face = "bold"),
        legend.position = "top") +
  scale_fill_manual(values = colors)

plot3

# Plot
ggplot(sow_counts, aes(x = as.factor(louvain_res1_5), y = Freq, fill = Sow)) +
  geom_bar(stat = "identity", position = "fill") +
  ylab("Count of Cells") + xlab("Louvain Cluster (res=1.5)") +
  ggtitle("Sow Distribution Across Louvain Clusters") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Build contingency table: Treatment × Sow
sow_counts_t <- as.data.frame(table(Treatment = sow_treat_df$Treatment,
                                    Sow = sow_treat_df$Sow))

# Plot
ggplot(sow_counts_t, aes(x = as.factor(Treatment), y = Freq, fill = Sow)) +
  geom_bar(stat = "identity", position = "stack") +
  ylab("Count of Cells") + xlab("Treatment") +
  ggtitle("Sow Distribution Across Treatment")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#number of cells ineach celltype in each treatment
table(seurat_14dpi$celltypes, seurat_14dpi$Treatment)






###FLOw- 2 comparsion to sc proprtions##
seurat_14dpi <- readRDS("/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Project_Fang_10X/filtered_postQC_postcb_postdoublet_postdowns_annotated_14dpi.rds")

Idents(seurat_14dpi) <- seurat_14dpi$celltypes
meta <- seurat_14dpi@meta.data %>%
  dplyr::select(Sample, Treatment, celltypes)
#get total cells per sample and persamplepercelltypes
scrna_props <- meta %>%
  add_count(Sample, name = "n_cells_total") %>%                    # total cells per sample
  add_count(Sample, celltypes, name = "n_cells") %>%               # cells per (sample, celltype)
  distinct(Sample, Treatment, celltypes, n_cells, n_cells_total) %>%
  mutate(prop_scrna = n_cells / n_cells_total)

all_samples  <- unique(scrna_props$Sample)
all_celltypes <- unique(scrna_props$celltypes)
#if no celltypes_sample ccombintaion then fill with 0
scrna_props <- scrna_props %>%
  complete(Sample = all_samples,
           celltypes = all_celltypes,
           fill = list(n_cells = 0, prop_scrna = 0)) %>%
  #add NA for incomplete ones
  group_by(Sample) %>%
  fill(Treatment, .direction = "downup") %>%
  ungroup()

write.csv(scrna_props,"/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/scrna_props_14dpi.csv", row.names = FALSE)

ggplot(scrna_props, aes(x = celltypes, y = Sample, fill = prop_scrna)) +
  geom_tile(color = "gray90", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.3f", prop_scrna)), size = 2.7, fontface = "bold") +
  scale_fill_gradientn(colors = c("blue", "white", "red"), name = "Proportion") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank()
  ) +
  labs(title = "scRNA cell-type proportions per sample- 14 dpi",
       x = "Cell type", y = "Sample")

####84 dpi
seurat_84dpi <- readRDS("/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/filtered_postQC_postcb_postdoublet_postdowns_postcellcycle_updated_annotation_84dpi.rds")

Idents(seurat_84dpi) <- seurat_84dpi$CellTypes
meta <- seurat_84dpi@meta.data %>%
  dplyr::select(Sample, Treatment, CellTypes)
#get total cells per sample and persamplepercelltypes
scrna_props <- meta %>%
  add_count(Sample, name = "n_cells_total") %>%                    # total cells per sample
  add_count(Sample, CellTypes, name = "n_cells") %>%               # cells per (sample, celltype)
  distinct(Sample, Treatment, CellTypes, n_cells, n_cells_total) %>%
  mutate(prop_scrna = n_cells / n_cells_total)

all_samples  <- unique(scrna_props$Sample)
all_celltypes <- unique(scrna_props$CellTypes)
#if no celltypes_sample ccombintaion then fill with 0
scrna_props <- scrna_props %>%
  complete(Sample = all_samples,
           CellTypes = all_celltypes,
           fill = list(n_cells = 0, prop_scrna = 0)) %>%
  #add NA for incomplete ones
  group_by(Sample) %>%
  fill(Treatment, .direction = "downup") %>%
  ungroup()

write.csv(scrna_props,"/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/scrna_props_84dpi.csv", row.names = FALSE)

ggplot(scrna_props, aes(x = CellTypes, y = Sample, fill = prop_scrna)) +
  geom_tile(color = "gray90", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.3f", prop_scrna)), size = 2.7, fontface = "bold") +
  scale_fill_gradientn(colors = c("blue", "white", "red"), name = "Proportion") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank()
  ) +
  labs(title = "scRNA cell-type proportions per sample- 84 dpi",
       x = "Cell type", y = "Sample")
###subset
seurat_84dpi_AB <- readRDS("/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/filtered_postQC_postcb_postdoublet_postdowns_postcellcycle_ABT_subset_84dpi.rds")
meta_AB   <- seurat_84dpi_AB@meta.data%>%
  dplyr::select(Sample, Treatment, CellTypes)

meta <- seurat_84dpi@meta.data %>%
  dplyr::select(Sample, Treatment, CellTypes)

# Count totals (sample_celltype) for the subset
subset_counts <- meta_AB %>%
  dplyr::count(Sample, CellTypes, name = "n_cells_subset")

# Count totals for main 84 dpi all cells
total_counts <- meta %>%
  dplyr::count(Sample, name = "total_cells_all")

# Merge to compute proportions relative to all cells
ab_props_vs_all <- subset_counts %>%
  left_join(total_counts, by = "Sample") %>%
  mutate(prop_of_all = n_cells_subset / total_cells_all)

head(ab_props_vs_all)
#zero fill 
all_samples<- unique(meta$Sample)
all_AB <- unique(meta_AB$CellTypes)

ab_props_vs_all <- ab_props_vs_all %>%
  tidyr::complete(Sample = all_samples,
                  CellTypes = all_AB,
                  fill = list(n_cells_subset = 0, prop_of_all = 0))
#plot

ggplot(ab_props_vs_all, aes(x = CellTypes, y = Sample, fill = prop_of_all)) +
  geom_tile(color = "gray90", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.3f", prop_of_all)), size = 2.7, fontface = "bold") +
  scale_fill_gradientn(colors = c("blue", "white", "red"), name = "Proportion") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid = element_blank()
  ) +
  labs(title = "AB T cell subset proportions relative to all cells",
       x = "Cell type", y = "Sample")


# 5 populations: GD, B, NK, Myeloid, AB

flow <- read.csv("/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/Flow_proportions.csv", sep= '')
sc_14 <- read.csv("/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/scrna_props_14dpi.csv")
sc_84 <- read.csv("/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/scrna_props_84dpi.csv")


flow14_5 <- flow %>%
  filter(dpi == 14, DataType == "PrcLv") %>%
  mutate(
    pop5 = case_when(
      CellType == "gdTCR"    ~ "GD",
      CellType == "AllBcells"~ "B",
      CellType == "NKcells"  ~ "NK",
      CellType == "CD172sp"  ~ "Myeloid", 
      CellType == "abTCR"    ~ "AB",
      TRUE ~ NA_character_
    ),
    prop_flow = Value / 100   
  ) %>%
  filter(!is.na(pop5)) %>%
  group_by(Sample, Treatment, pop5) %>%
  summarise(prop_flow = sum(prop_flow), .groups = "drop")
flow14_5 <- flow14_5 %>%
  mutate(
    Treatment = case_when(
      Treatment == "neg"     ~ "control",
      Treatment == "persist" ~ "persistent",
      Treatment == "extinct" ~ "extinct",
      TRUE ~ Treatment
    )
  )
valid_samples <- unique(sc_14$Sample)
flow14_5_filtered <- flow14_5 %>%
  filter(Sample %in% valid_samples)

flow84_5 <- flow %>%
  filter(dpi == 84, DataType == "PrcLv") %>%
  mutate(
    pop5 = case_when(
      CellType == "gdTCR"    ~ "GD",
      CellType == "AllBcells"~ "B",
      CellType == "NKcells"  ~ "NK",
      CellType == "CD172sp"  ~ "Myeloid", 
      CellType == "abTCR"    ~ "AB",
      TRUE ~ NA_character_
    ),
    prop_flow = Value / 100   
  ) %>%
  filter(!is.na(pop5)) %>%
  group_by(Sample, Treatment, pop5) %>%
  summarise(prop_flow = sum(prop_flow), .groups = "drop")
flow84_5 <- flow84_5 %>%
  mutate(
    Treatment = case_when(
      Treatment == "neg"     ~ "control",
      Treatment == "persist" ~ "persistent",
      Treatment == "extinct" ~ "extinct",
      TRUE ~ Treatment
    )
  )
sc14_5 <- sc_14 %>%
  mutate(
    pop5 = case_when(
      celltypes %in% c("CD2+ GD T cells" ,"CD2- GD T cells") ~ "GD",
      celltypes %in% c("B cells","ASCs") ~ "B",
      celltypes %in% c("NK cells") ~ "NK",
      celltypes %in% c("Monocytes","pDCs") ~ "Myeloid",
      TRUE ~ "AB"   
    )
  ) %>%
  filter(!is.na(pop5)) %>%  
  group_by(Sample, Treatment, pop5) %>%
  summarise(
    n = sum(n_cells, na.rm = TRUE),
    n_total = max(n_cells_total, na.rm = TRUE),
    prop_sc = n / n_total,
    .groups = "drop"
  )

sc84_5 <- sc_84 %>%
  mutate(
    pop5 = case_when(
      CellTypes %in% c("CD2+ GD T cells" ,"CD2- GD T cells") ~ "GD",
      CellTypes %in% c("B cells","ASC", "Transitional B-like") ~ "B",
      CellTypes %in% c("NK cells") ~ "NK",
      CellTypes %in% c("Monocytes","pDCs","cDCs") ~ "Myeloid",
      CellTypes %in% c("AB T cells") ~ "AB",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(pop5)) %>%          # <<< drop unmapped cell types
  group_by(Sample, Treatment, pop5) %>%
  summarise(
    n       = sum(n_cells, na.rm = TRUE),
    n_total = max(n_cells_total, na.rm = TRUE),
    prop_sc = n / n_total,
    .groups = "drop"
  )


#compare
comp14 <- inner_join(flow14_5, sc14_5,
                     by = c("Sample","Treatment","pop5"))
cor.test(comp14$prop_flow, comp14$prop_sc, method = 'spearman') 
cor_tests <- comp14 %>%
  group_by(pop5) %>%
  summarise(
    rho  = cor(prop_flow, prop_sc, method = "spearman", use = "complete.obs"),
    pval = cor.test(prop_flow, prop_sc, method = "spearman")$p.value,
    n    = n(),
    .groups = "drop"
  )

cor_tests
cor_by_trt <- comp14 %>%
  group_by(Treatment) %>%
  summarise(
    rho  = cor(prop_flow, prop_sc, method = "spearman", use = "complete.obs"),
    pval = cor.test(prop_flow, prop_sc, method = "spearman")$p.value,
    n    = n(),
    .groups = "drop"
  )

cor_by_trt

comp84 <- inner_join(flow84_5, sc84_5,
                     by = c("Sample","Treatment","pop5"))
cor.test(comp84$prop_flow, comp84$prop_sc, method = 'spearman') 
cor_tests_84 <- comp84 %>%
  group_by(pop5) %>%
  summarise(
    rho  = cor(prop_flow, prop_sc, method = "spearman", use = "complete.obs"),
    pval = cor.test(prop_flow, prop_sc, method = "spearman")$p.value,
    n    = n(),
    .groups = "drop"
  )

cor_tests_84
cor_by_trt_84 <- comp84 %>%
  group_by(Treatment) %>%
  summarise(
    rho  = cor(prop_flow, prop_sc, method = "spearman", use = "complete.obs"),
    pval = cor.test(prop_flow, prop_sc, method = "spearman")$p.value,
    n    = n(),
    .groups = "drop"
  )

cor_by_trt_84




###plotting ###
cor_tests_plot <- cor_tests%>%
  mutate(
    pop5 = factor(pop5, levels = c("AB", "B", "GD", "Myeloid", "NK")),
    sig = case_when(
      pval < 0.05  ~ "*"
    )
  )

ggplot(cor_tests_plot, aes(x = pop5, y = rho)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey70") +
  geom_bar(stat = "identity", aes(fill = pop5)) +
  geom_text(aes(label = sprintf("%.3f", rho)), vjust = -0.2, size = 3.5) +
  geom_text(aes(y = rho + 0.05, label = sig), vjust = 0, size = 4) +
  coord_cartesian(ylim = c(0, 1.1)) +
  theme_bw() +
  labs(
    title = "Correlation by celltype - 14 DPI",
    x = "Cellular Phenotype",
    y = "Spearman correlation"
  )+
  theme(
    axis.text.x = element_text( hjust = 1, face = 'bold'),
    axis.text.y = element_text(face = 'bold')
  )
##trt 
cor_trts_plot <- cor_by_trt%>%
  mutate(
    pop5 = factor(Treatment, levels = c("control", "extinct", "persistent")),
    sig = case_when(
      pval < 0.05  ~ "*"
    )
  )

ggplot(cor_trts_plot, aes(x = pop5, y = rho)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey70") +
  geom_bar(stat = "identity", aes(fill = Treatment)) +
  geom_text(aes(label = sprintf("%.3f", rho)), vjust = -0.2, size = 3.5) +
  geom_text(aes(y = rho + 0.05, label = sig), vjust = 0, size = 4) +
  coord_cartesian(ylim = c(0, 1.1)) +
  theme_bw() +
  labs(
    title = "Correlation by Treatment - 14 DPI",
    x = "Treatment",
    y = "Spearman correlation"
  )+
  theme(
    axis.text.x = element_text( hjust = 1, face = 'bold'),
    axis.text.y = element_text(face = 'bold')
  )



summ <- comp14 %>%
  group_by(Treatment, pop5) %>%
  summarise(
    flow_mean = mean(prop_flow, na.rm = TRUE),
    flow_sd   = sd(prop_flow,   na.rm = TRUE),
    sc_mean   = mean(prop_sc,   na.rm = TRUE),
    sc_sd     = sd(prop_sc,     na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(flow_mean, sc_mean, flow_sd, sc_sd),
    names_to = c("Method", ".value"),
    names_pattern = "(flow|sc)_(mean|sd)"
  )

ggplot(summ, aes(x = pop5, y = mean, color = Method)) +
  geom_point(
    position = position_dodge(width = 0.4),
    size = 3
  ) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd, group = Method),
    position = position_dodge(width = 0.4),
    width = 0.1
  ) +
  facet_wrap(~ Treatment) +
  theme_bw() +
  labs(
    title = "Trt level: 14 DPI",
    x = "Cellular Phenotypes",
    y = "Proportion",
    color = "Method"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = 'bold'),
    axis.text.y = element_text(face = 'bold')
  )


cor_txct <- comp84 %>%
  group_by(Treatment, pop5) %>%
  summarise(
    rho   = cor(prop_flow, prop_sc, method = "spearman", use = "complete.obs"),
    pval  = cor.test(prop_flow, prop_sc, method = "spearman")$p.value,
    MAE   = mean(abs(prop_sc - prop_flow), na.rm = TRUE),
    Ratio = mean(prop_sc / prop_flow,    na.rm = TRUE),
    n     = n(),
    .groups = "drop"
  )
cor_plot <- cor_txct %>%
  mutate(
    pop5 = factor(pop5, levels = c("AB", "B", "GD", "Myeloid", "NK")),
    MAE_cat = case_when(
      MAE < 0.05 ~ "Close (<0.05)",
      MAE < 0.10 ~ "Okay (0.05-0.1)",
      TRUE       ~ "poor (>0.1)"
    ),
    MAE_cat = factor(MAE_cat, levels = c("Close (<0.05)","Okay (0.05-0.1)","poor (>0.1)")),
    p_star = case_when(
      pval < 0.05  ~ "*",
      TRUE         ~ ""
    ),
    ratio_lab = sprintf("R=%.2f", Ratio)
  )

ggplot(cor_plot, aes(x = Treatment, y = pop5)) +
  geom_tile(fill = "grey95", color = "white") +
  geom_point(
    aes(fill = rho, shape = MAE_cat),
    size   = 15,   
    color  = "black",
    stroke = 0.5
  ) +

  geom_text(
    aes(label = p_star),
    color = "white",
    fontface = "bold",
    size = 5,
    vjust = -0.5
  ) +

  geom_text(
    aes(label = ratio_lab),
    color = "white",
    size = 3,
    fontface = "bold",
    vjust = 1.8   
  ) +
  
  scale_fill_gradient2(
    low  = "steelblue",
    mid  = "white",
    high = "firebrick",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman ρ"
  ) +

  scale_shape_manual(
    values = c(
      "Close (<0.05)"    = 21,   # big circle
      "Okay (0.05-0.1)" = 22,   # square
      "poor (>0.1)"   = 24    # triangle
    ),
    name = "MAE level"
  ) +
  
  theme_bw() +
  labs(
    title    = "scRNA vs Flow - Day 84",
    subtitle = "Fill = correlation, Shape = MAE, Star = p-value, Ration = sc/flow ",
    x        = "Treatment",
    y        = "Cell type"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid  = element_blank()
  )

write.csv(cor_plot, file="/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/day84_sc_flow_plot_metrics.csv")


comp14$n_total <- NULL; comp14$n <- NULL
write.csv(comp14, file="/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/day14_sc_flow.csv")

comp84$n_total <- NULL; comp84$n <- NULL
write.csv(comp84, file="/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/day84_sc_flow.csv")
