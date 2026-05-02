#######SCRIPT for GO analysis using cluter-profiler###########
.libPaths("/work/ABG/mkapoor/mkapoor/.ondemand-new/mkapoor/rstudio/libs/4.4.1")

library(enrichplot)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(msigdbr)
library(org.Hs.eg.db)
library(DOSE)
library(clusterProfiler)


###compare if ExC mono have more enriched terms than all genes exp in the data#####
#day 14#
seurat_14dpi <- readRDS("/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Project_Fang_10X/filtered_postQC_postcb_postdoublet_postdowns_annotated_14dpi.rds")
#cconvert to human ortholog gene names
orthoGenes <- read.delim("/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/PigToHuman_GeneOrthos_v97.txt") # read in gene ortholog file
orthoGenes <- subset(orthoGenes, Human.homology.type == 'ortholog_one2one') # subset to only one to one orthologs
genes <- as.data.frame(rownames(seurat_14dpi[['RNA']]@data)) # extract pig gene names from dataset
colnames(genes) <- 'gene'
genes <- data.frame(gene = genes)

pigGenes <- read_delim('/work/ABG/mkapoor/mkapoor/PIPseq_Salmonella/2023_oct_seq/gtf/Sus_scrofa.Sscrofa11.1.97_modified06302021_JEW_SKS.csv' ,) # read in file with an updated gene symbol annotation for Sus scrofa v97 annotation build
pigGenes$FinalList <-gsub("_", "-", pigGenes$gene_name) # replace all underscores with dashes since this occurred when processing data in a previous step , total length:1292513
pigGenes <- pigGenes[pigGenes$FinalList %in% genes$gene, ] # slim down to only genes in our dataset
orthos <- intersect(pigGenes$gene_id, orthoGenes$Gene.stable.ID) # find which genes are one-to-one orthologs
#length(orthos) # how many genes are orthologs?- 16147
pigGenes <- pigGenes[pigGenes$gene_id %in% orthos, ]
#dim(pigGenes) - 1086339      25
orthoGenes <- orthoGenes[orthoGenes$Gene.stable.ID %in% pigGenes$gene_id, ] # slim down to only ortho genes in our dataset
orthoGenes <- orthoGenes %>% distinct(orthoGenes$Gene.stable.ID, orthoGenes$Human.gene.stable.ID, .keep_all = TRUE)  
# Merge pigGenes with orthoGenes to get human gene symbols
pig2human <- merge(pigGenes, orthoGenes, by.x = "gene_id", by.y = "Gene.stable.ID")
colnames(pig2human)[which(colnames(pig2human) == "Human.gene.symbol")] <- "human_gene_symbol"
############################
###GO analysis #####
############################
expressed_genes <- rownames(seurat_14dpi)[rowSums(seurat_14dpi@assays$RNA@counts > 0) > 0]
expressed_genes_mapped <- pig2human %>%
  filter(FinalList %in% expressed_genes) %>%
  pull(Human.gene.name) %>%
  unique()

#day 84#
seurat_84dpi <- readRDS("/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/filtered_postQC_postcb_postdoublet_postdowns_postcellcycle_updated_annotation_84dpi.rds")
#convert to human ortholog gene names
orthoGenes <- read.delim("/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/PigToHuman_GeneOrthos_v97.txt") # read in gene ortholog file
orthoGenes <- subset(orthoGenes, Human.homology.type == 'ortholog_one2one') # subset to only one to one orthologs
genes <- as.data.frame(rownames(seurat_84dpi[['RNA']]@data)) # extract pig gene names from dataset
colnames(genes) <- 'gene'
genes <- data.frame(gene = genes)

pigGenes <- read_delim('/work/ABG/mkapoor/mkapoor/PIPseq_Salmonella/2023_oct_seq/gtf/Sus_scrofa.Sscrofa11.1.97_modified06302021_JEW_SKS.csv' ,) # read in file with an updated gene symbol annotation for Sus scrofa v97 annotation build
pigGenes$FinalList <-gsub("_", "-", pigGenes$gene_name) # replace all underscores with dashes since this occurred when processing data in a previous step , total length:1292513
pigGenes <- pigGenes[pigGenes$FinalList %in% genes$gene, ] # slim down to only genes in our dataset
orthos <- intersect(pigGenes$gene_id, orthoGenes$Gene.stable.ID) # find which genes are one-to-one orthologs
#length(orthos) # how many genes are orthologs?- 16147
pigGenes <- pigGenes[pigGenes$gene_id %in% orthos, ]
#dim(pigGenes) - 1086339      25
orthoGenes <- orthoGenes[orthoGenes$Gene.stable.ID %in% pigGenes$gene_id, ] # slim down to only ortho genes in our dataset
orthoGenes <- orthoGenes %>% distinct(orthoGenes$Gene.stable.ID, orthoGenes$Human.gene.stable.ID, .keep_all = TRUE)  
# Merge pigGenes with orthoGenes to get human gene symbols
pig2human <- merge(pigGenes, orthoGenes, by.x = "gene_id", by.y = "Gene.stable.ID")
colnames(pig2human)[which(colnames(pig2human) == "Human.gene.symbol")] <- "human_gene_symbol"
expressed_genes <- rownames(seurat_84dpi)[rowSums(seurat_84dpi@assays$RNA@counts > 0) > 0]
expressed_genes_mapped_84 <- pig2human %>%
  filter(FinalList %in% expressed_genes) %>%
  pull(Human.gene.name) %>%
  unique()
############################
######Extinct vs Peristent####
###########################

#MONOCYTES #
de_results_mono<- read.csv("/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Project_Fang_10X/MAST_DE_14dpi/Monocytes/Mast_Mono_ext_per.csv")
de_results_mono<- read.csv("/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/MAST_DE_84dpi/Monocytes/Mast_Monocytes_ext_per_sigFDR.csv")
colnames(de_results_mono)[2] <- "gene"
de_results_mono$gene <- as.character(de_results_mono$gene)
de_results_mono <- na.omit(de_results_mono)
# check to see if format is correct
head(de_results_mono)
dim(de_results_mono)
if ("direction" %in% colnames(de_results_mono)) {
  de_results_mono$direction <- NULL
}
de_results_mono$direction <- ifelse(de_results_mono$fdr < 0.05 & de_results_mono$logFC > 0.05, "Upregulated",
                                    ifelse(de_results_mono$fdr < 0.05 & de_results_mono$logFC < -0.05, "Downregulated", "Non-significant"))
head(de_results_mono)
#Join DEG genes with pig-to-human map
shared_up <- de_results_mono %>%
  filter(direction == "Upregulated")
shared_up_human <- shared_up %>%
  left_join(pig2human, by = c("gene" = "FinalList")) %>%
  filter(!is.na(Human.gene.name)) %>%
  pull(Human.gene.name) %>%
  unique()

# Shared downregulated genes
shared_down <- de_results_mono %>%
  filter(direction == "Downregulated")
shared_down_human <- shared_down %>%
  left_join(pig2human, by = c("gene" = "FinalList")) %>%
  filter(!is.na(Human.gene.name)) %>%
  pull(Human.gene.name) %>%
  unique()
#enrichment
comparelist <- list(shared_up = shared_up_human,
                    shared_down = shared_down_human)
names(comparelist)<-c("14_up_E/P","14_down_E/P") 
cclust<-compareCluster(geneCluster = comparelist, 
                       fun = enrichGO,
                       OrgDb= org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont= "BP",pAdjustMethod = "BH",
                       qvalueCutoff  = 0.05,
                       universe=expressed_genes_mapped_84,  minGSSize = 10)
clust_simple <- clusterProfiler::simplify(cclust, cutoff=0.35, by="p.adjust", select_fun=min) 
GO_results_minGsize5 <- data.frame(clust_simple)
GO_results_minGsize5_no_simplify <- data.frame(cclust)
go_genes <- clust_simple@compareClusterResult %>%
  dplyr::select(ID, Description, geneID, Cluster,p.adjust,Count) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(GO_ID = ID, Gene = geneID)
go_magnitude <- go_genes %>%
  left_join(de_results_mono, by = c("Gene" = "gene")) %>%
  group_by(GO_ID, Description, Cluster,p.adjust,Count) %>%
  summarise(
    N_genes = n(),
    median_abs_logFC = median(abs(logFC), na.rm = TRUE)
  ) %>%
  arrange(desc(median_abs_logFC))
top_terms <- go_magnitude %>% top_n(35, median_abs_logFC)


p<-ggplot(top_terms, aes(x = median_abs_logFC,
                         y = reorder(Description, median_abs_logFC),
                         color = p.adjust,
                         size = Count)) +
  geom_point() +
  facet_wrap(~Cluster) +
  scale_color_gradient(low = "red", high = "yellow", name = "FDR (p.adjust)") +
  scale_size_continuous(name = "Gene Count") +
  labs(
    x = "Median |log2FC|",
    y = "GO term"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10, face = "bold"))
p
write.csv(GO_results_minGsize5, "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Project_Fang_10X/GSEA_GO/GO_Mono/GO_terms/mono_GO_EP.csv", row.names = FALSE)
ggsave(file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Project_Fang_10X/GSEA_GO/GO_Mono/GO_terms/mono_EP_GO.png", p, width = 11, height = 5, units = "in")

#run pathway
library(msigdbr)
m_t2g_symbol <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
head(m_t2g_symbol)
cclust <- compareCluster(
  geneCluster   = comparelist,
  fun           = "enricher",
  TERM2GENE     = m_t2g_symbol,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)
pathway <- data.frame(cclust)
path_genes <- cclust@compareClusterResult %>%
  dplyr::select(ID, Description, geneID, Cluster,p.adjust,Count) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(GO_ID = ID, Gene = geneID)
path_magnitude <- path_genes %>%
  left_join(de_results_mono, by = c("Gene" = "gene")) %>%
  group_by(GO_ID, Description, Cluster,p.adjust,Count) %>%
  summarise(
    N_genes = n(),
    mean_abs_logFC = mean(abs(logFC), na.rm = TRUE),
    median_abs_logFC = median(abs(logFC), na.rm = TRUE)
  ) %>%
  arrange(desc(median_abs_logFC))
#write.csv(go_magnitude, "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/compare_14_84/GO_unique_genes/GO_term_mag_d14_minGsize20_sim0.35", row.names = FALSE)
top_terms <- path_magnitude %>% top_n(35, median_abs_logFC)


p<-ggplot(top_terms, aes(x = median_abs_logFC,
                         y = reorder(Description, median_abs_logFC),
                         color = p.adjust,
                         size = Count)) +
  geom_point() +
  facet_wrap(~Cluster) +
  scale_color_gradient(low = "red", high = "yellow", name = "FDR (p.adjust)") +
  scale_size_continuous(name = "Gene Count") +
  labs(
    x = "Median |log2FC|",
    y = "MsigDB term"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10, face = "bold"))
write.csv(pathway, "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Project_Fang_10X/GSEA_GO/GO_Mono/GO_terms/mono_EP_GO_path.csv", row.names = FALSE)
ggsave(file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Project_Fang_10X/GSEA_GO/GO_Mono/GO_terms/mono_EP_GO_path.png", p, width = 11, height = 5, units = "in")

go_df <- as.data.frame(clust_simple)
go_df$FirstGene <- sapply(strsplit(as.character(go_df$geneID), "/"), function(x) {paste(head(x,5), collapse =",")})
top_terms <- go_df[order(go_df$p.adjust), ][1:min(30, nrow(go_df)), ]
up_df   <- go_df %>% filter(grepl("up",   Cluster, ignore.case = TRUE))
down_df <- go_df %>% filter(grepl("down", Cluster, ignore.case = TRUE))
up_df$Description   <- with(up_df, reorder(Description,  -log10(p.adjust)))
down_df$Description <- with(down_df, reorder(Description, -log10(p.adjust)))

p <- ggplot() +
geom_point(
  data = up_df,
  aes(x = -log10(p.adjust),
      y = Description,
      size = Count,
      color = -log10(p.adjust))
) +
  scale_color_gradient(low = "orange", high = "darkred", name = "-log10(FDR)") +
  ggnewscale::new_scale_color() + 
geom_point(
  data = down_df,
  aes(x = -log10(p.adjust),
      y = Description,
      size = Count,
      color = -log10(p.adjust))
) +
  scale_color_gradient(low = "green", high = "darkgreen", name = "-log10(FDR)") +
  ggrepel::geom_text_repel(
    data = bind_rows(up_df, down_df),
    aes(x = -log10(p.adjust),
        y = Description,
        label = FirstGene),
    size = 4,force = 2,
    max.overlaps = 50,
    box.padding = 0.6,point.padding = 0.5,
    segment.color = "grey50"
  ) +
  
  facet_wrap(~ Cluster) +
  scale_size_continuous(name = "Gene Count") +
  labs(x = "-log10(FDR)", y = "GO Term") +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.y = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold")
  )

p
ggsave(file = "/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/GSEA_GO/GO_mono/GO_terms/mono_PC_GO_path_v2.png", plot, width = 12, height = 5, units = "in")


############################
######Extinct vs Control####
###########################

#MONOCYTES #
de_results_mono<- read.csv("/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Project_Fang_10X/MAST_DE_14dpi/Monocytes/Mast_Mono_ext_con.csv")
colnames(de_results_mono)[2] <- "gene"
de_results_mono$gene <- as.character(de_results_mono$gene)
de_results_mono <- na.omit(de_results_mono)
# check to see if format is correct
head(de_results_mono)
dim(de_results_mono)
if ("direction" %in% colnames(de_results_mono)) {
  de_results_mono$direction <- NULL
}
de_results_mono$direction <- ifelse(de_results_mono$fdr < 0.05 & de_results_mono$logFC > 0.05, "Upregulated",
                                    ifelse(de_results_mono$fdr < 0.05 & de_results_mono$logFC < -0.05, "Downregulated", "Non-significant"))
head(de_results_mono)
#Join DEG genes with pig-to-human map
shared_up <- de_results_mono %>%
  filter(direction == "Upregulated")
shared_up_human <- shared_up %>%
  left_join(pig2human, by = c("gene" = "FinalList")) %>%
  filter(!is.na(Human.gene.name)) %>%
  pull(Human.gene.name) %>%
  unique()

# Shared downregulated genes
shared_down <- de_results_mono %>%
  filter(direction == "Downregulated")
shared_down_human <- shared_down %>%
  left_join(pig2human, by = c("gene" = "FinalList")) %>%
  filter(!is.na(Human.gene.name)) %>%
  pull(Human.gene.name) %>%
  unique()
go_results_down_mono <- enrichGO(
  gene          = unique(shared_down_human),
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",         # Use gene symbols or other appropriate key
  ont           = "ALL",             # Biological Process ontology
  universe      = expressed_genes_mapped,  # Background
  pAdjustMethod = "BH",             # Adjust p-values using Benjamini-Hochberg
  qvalueCutoff  = 0.05, ,  minGSSize = 10)              # Significance threshold
go_results_down_mono_s <- simplify(go_results_down_mono, cutoff=0.35, by="p.adjust", select_fun=min) 
go_results_up_mono <- enrichGO(
  gene          = unique(shared_up_human),
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",         # Use gene symbols or other appropriate key
  ont           = "ALL",             # Biological Process ontology
  universe      = expressed_genes_mapped,  # Background
  pAdjustMethod = "BH",             # Adjust p-values using Benjamini-Hochberg
  qvalueCutoff  = 0.05, ,  minGSSize = 10)              # Significance threshold
go_results_up_mono_s <- simplify(go_results_up_mono, cutoff=0.35, by="p.adjust", select_fun=min) 

go_genes <- go_results_up_mono_s@result %>%
  dplyr::select(ID, Description, geneID,p.adjust,Count) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(GO_ID = ID, Gene = geneID)
go_magnitude <- go_genes %>%
  left_join(de_results_mono, by = c("Gene" = "gene")) %>%
  group_by(GO_ID, Description,p.adjust,Count) %>%
  summarise(
    N_genes = n(),
    mean_abs_logFC = mean(abs(logFC), na.rm = TRUE),
    median_abs_logFC = median(abs(logFC), na.rm = TRUE)
  ) %>%
  arrange(desc(median_abs_logFC))
top_terms <- go_magnitude %>% top_n(35, median_abs_logFC)


ggplot(top_terms, aes(x = median_abs_logFC,
                      y = reorder(Description, median_abs_logFC),
                      color = p.adjust,
                      size = Count)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "yellow", name = "FDR (p.adjust)") +
  scale_size_continuous(name = "Gene Count") +
  labs(
    title = "Magnitude of transcriptional change per GO term (D14)",
    x = "Median |log2FC|",
    y = "GO term"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10, face = "bold"))
#enrichment
comparelist <- list(shared_up = shared_up_human,
                    shared_down = shared_down_human)
names(comparelist)<-c("14_up_E/C","14_down_E/C") 
cclust<-compareCluster(geneCluster = comparelist, 
                       fun = enrichGO,
                       OrgDb= org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont= "BP",pAdjustMethod = "BH",
                       qvalueCutoff  = 0.05,
                       universe=expressed_genes_mapped,  minGSSize = 10)
clust_simple <- simplify(cclust, cutoff=0.35, by="p.adjust", select_fun=min) 
GO_results_minGsize5 <- data.frame(cclust_simple)
GO_results_minGsize5_no_simplify <- data.frame(cclust)
go_genes <- clust_simple@compareClusterResult %>%
  dplyr::select(ID, Description, geneID, Cluster,p.adjust,Count) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(GO_ID = ID, Gene = geneID)
go_magnitude <- go_genes %>%
  left_join(de_results_mono, by = c("Gene" = "gene")) %>%
  group_by(GO_ID, Description, Cluster,p.adjust,Count) %>%
  summarise(
    N_genes = n(),
    mean_abs_logFC = mean(abs(logFC), na.rm = TRUE),
    median_abs_logFC = median(abs(logFC), na.rm = TRUE)
  ) %>%
  arrange(desc(median_abs_logFC))
top_terms <- go_magnitude %>% top_n(35, median_abs_logFC)


ggplot(top_terms, aes(x = median_abs_logFC,
                      y = reorder(Description, median_abs_logFC),
                      color = p.adjust,
                      size = Count)) +
  geom_point() +
  facet_wrap(~Cluster) +
  scale_color_gradient(low = "red", high = "yellow", name = "FDR (p.adjust)") +
  scale_size_continuous(name = "Gene Count") +
  labs(
    title = "Magnitude of transcriptional change per GO term (D14)",
    x = "Median |log2FC|",
    y = "GO term"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10, face = "bold"))
#run pathway
library(msigdbr)
m_t2g_symbol <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
head(m_t2g_symbol)
cclust <- compareCluster(
  geneCluster   = comparelist,
  fun           = "enricher",
  TERM2GENE     = m_t2g_symbol,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

path_genes <- cclust@compareClusterResult %>%
  dplyr::select(ID, Description, geneID, Cluster,p.adjust,Count) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(GO_ID = ID, Gene = geneID)
path_magnitude <- path_genes %>%
  left_join(de_results_mono, by = c("Gene" = "gene")) %>%
  group_by(GO_ID, Description, Cluster,p.adjust,Count) %>%
  summarise(
    N_genes = n(),
    mean_abs_logFC = mean(abs(logFC), na.rm = TRUE),
    median_abs_logFC = median(abs(logFC), na.rm = TRUE)
  ) %>%
  arrange(desc(median_abs_logFC))
#write.csv(go_magnitude, "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/compare_14_84/GO_unique_genes/GO_term_mag_d14_minGsize20_sim0.35", row.names = FALSE)
top_terms <- path_magnitude %>% top_n(35, median_abs_logFC)


ggplot(top_terms, aes(x = median_abs_logFC,
                      y = reorder(Description, median_abs_logFC),
                      color = p.adjust,
                      size = Count)) +
  geom_point() +
  facet_wrap(~Cluster) +
  scale_color_gradient(low = "red", high = "yellow", name = "FDR (p.adjust)") +
  scale_size_continuous(name = "Gene Count") +
  labs(
    title = "Magnitude of transcriptional change per Hallmark term (D14)",
    x = "Median |log2FC|",
    y = "MsigDB term"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10, face = "bold"))



############################
######Persistent vs Control####
###########################

#MONOCYTES #
de_results_mono<- read.csv("/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Project_Fang_10X/MAST_DE_14dpi/Monocytes/Mast_Mono_per_con.csv")
colnames(de_results_mono)[2] <- "gene"
de_results_mono$gene <- as.character(de_results_mono$gene)
de_results_mono <- na.omit(de_results_mono)
# check to see if format is correct
head(de_results_mono)
dim(de_results_mono)
if ("direction" %in% colnames(de_results_mono)) {
  de_results_mono$direction <- NULL
}
de_results_mono$direction <- ifelse(de_results_mono$fdr < 0.05 & de_results_mono$logFC > 0.05, "Upregulated",
                                    ifelse(de_results_mono$fdr < 0.05 & de_results_mono$logFC < -0.05, "Downregulated", "Non-significant"))
head(de_results_mono)
#Join DEG genes with pig-to-human map
shared_up <- de_results_mono %>%
  filter(direction == "Upregulated")
shared_up_human <- shared_up %>%
  left_join(pig2human, by = c("gene" = "FinalList")) %>%
  filter(!is.na(Human.gene.name)) %>%
  pull(Human.gene.name) %>%
  unique()

# Shared downregulated genes
shared_down <- de_results_mono %>%
  filter(direction == "Downregulated")
shared_down_human <- shared_down %>%
  left_join(pig2human, by = c("gene" = "FinalList")) %>%
  filter(!is.na(Human.gene.name)) %>%
  pull(Human.gene.name) %>%
  unique()
#enrichment
comparelist <- list(shared_up = shared_up_human,
                    shared_down = shared_down_human)
names(comparelist)<-c("14_up_P/C","14_down_P/C") 
cclust<-compareCluster(geneCluster = comparelist, 
                       fun = enrichGO,
                       OrgDb= org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont= "BP",pAdjustMethod = "BH",
                       qvalueCutoff  = 0.05,
                       universe=expressed_genes_mapped,  minGSSize = 10)
clust_simple <- simplify(cclust, cutoff=0.35, by="p.adjust", select_fun=min) 
GO_results_minGsize5 <- data.frame(clust_simple)
GO_results_minGsize5_no_simplify <- data.frame(cclust)
go_genes <- clust_simple@compareClusterResult %>%
  dplyr::select(ID, Description, geneID, Cluster,p.adjust,Count) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(GO_ID = ID, Gene = geneID)
go_magnitude <- go_genes %>%
  left_join(de_results_mono, by = c("Gene" = "gene")) %>%
  group_by(GO_ID, Description, Cluster,p.adjust,Count) %>%
  summarise(
    N_genes = n(),
    mean_abs_logFC = mean(abs(logFC), na.rm = TRUE),
    median_abs_logFC = median(abs(logFC), na.rm = TRUE)
  ) %>%
  arrange(desc(median_abs_logFC))
top_terms <- go_magnitude %>% top_n(35, median_abs_logFC)


ggplot(top_terms, aes(x = median_abs_logFC,
                      y = reorder(Description, median_abs_logFC),
                      color = p.adjust,
                      size = Count)) +
  geom_point() +
  facet_wrap(~Cluster) +
  scale_color_gradient(low = "red", high = "yellow", name = "FDR (p.adjust)") +
  scale_size_continuous(name = "Gene Count") +
  labs(
    title = "Magnitude of transcriptional change per GO term (D14)",
    x = "Median |log2FC|",
    y = "GO term"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10, face = "bold"))
#run pathway
library(msigdbr)
m_t2g_symbol <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
head(m_t2g_symbol)
cclust <- compareCluster(
  geneCluster   = comparelist,
  fun           = "enricher",
  TERM2GENE     = m_t2g_symbol,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

path_genes <- cclust@compareClusterResult %>%
  dplyr::select(ID, Description, geneID, Cluster,p.adjust,Count) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(GO_ID = ID, Gene = geneID)
path_magnitude <- path_genes %>%
  left_join(de_results_mono, by = c("Gene" = "gene")) %>%
  group_by(GO_ID, Description, Cluster,p.adjust,Count) %>%
  summarise(
    N_genes = n(),
    mean_abs_logFC = mean(abs(logFC), na.rm = TRUE),
    median_abs_logFC = median(abs(logFC), na.rm = TRUE)
  ) %>%
  arrange(desc(median_abs_logFC))
#write.csv(go_magnitude, "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/compare_14_84/GO_unique_genes/GO_term_mag_d14_minGsize20_sim0.35", row.names = FALSE)
top_terms <- path_magnitude %>% top_n(35, median_abs_logFC)


ggplot(top_terms, aes(x = median_abs_logFC,
                      y = reorder(Description, median_abs_logFC),
                      color = p.adjust,
                      size = Count)) +
  geom_point() +
  facet_wrap(~Cluster) +
  scale_color_gradient(low = "red", high = "yellow", name = "FDR (p.adjust)") +
  scale_size_continuous(name = "Gene Count") +
  labs(
    title = "Magnitude of transcriptional change per Hallmark term (D14)",
    x = "Median |log2FC|",
    y = "MsigDB term"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10, face = "bold"))

############################
######Extinct vs Persistent - 84####
###########################

#MONOCYTES #
de_results_mono<- read.csv("/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/MAST_DE_84dpi/Monocytes/Mast_Monocytes_ext_per_sigFDR.csv")
colnames(de_results_mono)[2] <- "gene"
de_results_mono$gene <- as.character(de_results_mono$gene)
de_results_mono <- na.omit(de_results_mono)
# check to see if format is correct
head(de_results_mono)
dim(de_results_mono)
if ("direction" %in% colnames(de_results_mono)) {
  de_results_mono$direction <- NULL
}
de_results_mono$direction <- ifelse(de_results_mono$fdr < 0.05 & de_results_mono$logFC > 0.05, "Upregulated",
                                    ifelse(de_results_mono$fdr < 0.05 & de_results_mono$logFC < -0.05, "Downregulated", "Non-significant"))
head(de_results_mono)
dim(de_results_mono)
#Join DEG genes with pig-to-human map
shared_up <- de_results_mono %>%
  filter(direction == "Upregulated")
shared_up_human <- shared_up %>%
  left_join(pig2human, by = c("gene" = "FinalList")) %>%
  filter(!is.na(Human.gene.name)) %>%
  pull(Human.gene.name) %>%
  unique()

# Shared downregulated genes
shared_down <- de_results_mono %>%
  filter(direction == "Downregulated")
shared_down_human <- shared_down %>%
  left_join(pig2human, by = c("gene" = "FinalList")) %>%
  filter(!is.na(Human.gene.name)) %>%
  pull(Human.gene.name) %>%
  unique()
#enrichment
comparelist <- list(shared_up = shared_up_human,
                    shared_down = shared_down_human)
names(comparelist)<-c("84_up_E/P","84_down_E/P") 
cclust<-compareCluster(geneCluster = comparelist, 
                       fun = enrichGO,
                       OrgDb= org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont= "BP",pAdjustMethod = "BH",
                       qvalueCutoff  = 0.05,
                       universe=expressed_genes_mapped_84,  minGSSize = 10)
clust_simple <- clusterProfiler::simplify(cclust, cutoff=0.35, by="p.adjust", select_fun=min) 
GO_results_minGsize5 <- data.frame(clust_simple)
GO_results_minGsize5_no_simplify <- data.frame(cclust)
go_genes <- clust_simple@compareClusterResult %>%
  dplyr::select(ID, Description, geneID, Cluster,p.adjust,Count) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(GO_ID = ID, Gene = geneID)
go_magnitude <- go_genes %>%
  left_join(de_results_mono, by = c("Gene" = "gene")) %>%
  group_by(GO_ID, Description, Cluster,p.adjust,Count) %>%
  summarise(
    N_genes = n(),
    median_abs_logFC = median(abs(logFC), na.rm = TRUE)
  ) %>%
  arrange(desc(median_abs_logFC))
top_terms <- go_magnitude %>% top_n(35, median_abs_logFC)


p<-ggplot(top_terms, aes(x = median_abs_logFC,
                         y = reorder(Description, median_abs_logFC),
                         color = p.adjust,
                         size = Count)) +
  geom_point() +
  facet_wrap(~Cluster) +
  scale_color_gradient(low = "red", high = "yellow", name = "FDR (p.adjust)") +
  scale_size_continuous(name = "Gene Count") +
  labs(
    x = "Median |log2FC|",
    y = "GO term"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10, face = "bold"))
write.csv(GO_results_minGsize5, "/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/GSEA_GO/GO_mono/GO_terms/mono_EP_GO.csv", row.names = FALSE)
ggsave(file = "/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/GSEA_GO/GO_mono/GO_terms/mono_EP_GO.png", p, width = 11, height = 5, units = "in")

#run pathway
library(msigdbr)
m_t2g_symbol <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
head(m_t2g_symbol)
cclust <- compareCluster(
  geneCluster   = comparelist,
  fun           = "enricher",
  TERM2GENE     = m_t2g_symbol,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

path_genes <- cclust@compareClusterResult %>%
  dplyr::select(ID, Description, geneID, Cluster,p.adjust,Count) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(GO_ID = ID, Gene = geneID)
path_magnitude <- path_genes %>%
  left_join(de_results_mono, by = c("Gene" = "gene")) %>%
  group_by(GO_ID, Description, Cluster,p.adjust,Count) %>%
  summarise(
    N_genes = n(),
    mean_abs_logFC = mean(abs(logFC), na.rm = TRUE),
    median_abs_logFC = median(abs(logFC), na.rm = TRUE)
  ) %>%
  arrange(desc(median_abs_logFC))
#write.csv(go_magnitude, "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/compare_14_84/GO_unique_genes/GO_term_mag_d14_minGsize20_sim0.35", row.names = FALSE)
top_terms <- path_magnitude %>% top_n(35, median_abs_logFC)


p<-ggplot(top_terms, aes(x = median_abs_logFC,
                         y = reorder(Description, median_abs_logFC),
                         color = p.adjust,
                         size = Count)) +
  geom_point() +
  facet_wrap(~Cluster) +
  scale_color_gradient(low = "red", high = "yellow", name = "FDR (p.adjust)") +
  scale_size_continuous(name = "Gene Count") +
  labs(
    x = "Median |log2FC|",
    y = "MsigDB term"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10, face = "bold"))
write.csv(top_terms, "/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/GSEA_GO/GO_mono/GO_terms/mono_PC_GO_path.csv", row.names = FALSE)
ggsave(file = "/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/GSEA_GO/GO_mono/GO_terms/mono_PC_GO_path.png", p, width = 11, height = 5, units = "in")

############################
######Persistent vs Control- 84####
###########################

#MONOCYTES #
de_results_mono<- read.csv("/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/MAST_DE_84dpi/Monocytes/Mast_Monocytes_per_con_sigFDR.csv")
colnames(de_results_mono)[2] <- "gene"
de_results_mono$gene <- as.character(de_results_mono$gene)
de_results_mono <- na.omit(de_results_mono)
# check to see if format is correct
head(de_results_mono)
dim(de_results_mono)
if ("direction" %in% colnames(de_results_mono)) {
  de_results_mono$direction <- NULL
}
de_results_mono$direction <- ifelse(de_results_mono$fdr < 0.05 & de_results_mono$logFC > 0.05, "Upregulated",
                                    ifelse(de_results_mono$fdr < 0.05 & de_results_mono$logFC < -0.05, "Downregulated", "Non-significant"))
head(de_results_mono)
dim(de_results_mono)
#Join DEG genes with pig-to-human map
shared_up <- de_results_mono %>%
  filter(direction == "Upregulated")
shared_up_human <- shared_up %>%
  left_join(pig2human, by = c("gene" = "FinalList")) %>%
  filter(!is.na(Human.gene.name)) %>%
  pull(Human.gene.name) %>%
  unique()

# Shared downregulated genes
shared_down <- de_results_mono %>%
  filter(direction == "Downregulated")
shared_down_human <- shared_down %>%
  left_join(pig2human, by = c("gene" = "FinalList")) %>%
  filter(!is.na(Human.gene.name)) %>%
  pull(Human.gene.name) %>%
  unique()
#enrichment
comparelist <- list(shared_up = shared_up_human,
                    shared_down = shared_down_human)
names(comparelist)<-c("84_up_P/C","84_down_P/C") 
cclust<-compareCluster(geneCluster = comparelist, 
                       fun = enrichGO,
                       OrgDb= org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont= "BP",pAdjustMethod = "BH",
                       qvalueCutoff  = 0.05,
                       universe=expressed_genes_mapped_84,  minGSSize = 10)
clust_simple <- simplify(cclust, cutoff=0.35, by="p.adjust", select_fun=min) 
GO_results_minGsize5 <- data.frame(clust_simple)
GO_results_minGsize5_no_simplify <- data.frame(cclust)
go_genes <- clust_simple@compareClusterResult %>%
  dplyr::select(ID, Description, geneID, Cluster,p.adjust,Count) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(GO_ID = ID, Gene = geneID)
go_magnitude <- go_genes %>%
  left_join(de_results_mono, by = c("Gene" = "gene")) %>%
  group_by(GO_ID, Description, Cluster,p.adjust,Count) %>%
  summarise(
    N_genes = n(),
    median_abs_logFC = median(abs(logFC), na.rm = TRUE)
  ) %>%
  arrange(desc(median_abs_logFC))
top_terms <- go_magnitude %>% top_n(35, median_abs_logFC)


p<-ggplot(top_terms, aes(x = median_abs_logFC,
                      y = reorder(Description, median_abs_logFC),
                      color = p.adjust,
                      size = Count)) +
  geom_point() +
  facet_wrap(~Cluster) +
  scale_color_gradient(low = "red", high = "yellow", name = "FDR (p.adjust)") +
  scale_size_continuous(name = "Gene Count") +
  labs(
    x = "Median |log2FC|",
    y = "GO term"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10, face = "bold"))
write.csv(top_terms, "/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/GSEA_GO/GO_mono/GO_terms/mono_PC_GO.csv", row.names = FALSE)
ggsave(file = "/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/GSEA_GO/GO_mono/GO_terms/mono_PC_GO.png", p, width = 11, height = 5, units = "in")

#run pathway
library(msigdbr)
m_t2g_symbol <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
head(m_t2g_symbol)
cclust <- compareCluster(
  geneCluster   = comparelist,
  fun           = "enricher",
  TERM2GENE     = m_t2g_symbol,
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

path_genes <- cclust@compareClusterResult %>%
  dplyr::select(ID, Description, geneID, Cluster,p.adjust,Count) %>%
  tidyr::separate_rows(geneID, sep = "/") %>%
  dplyr::rename(GO_ID = ID, Gene = geneID)
path_magnitude <- path_genes %>%
  left_join(de_results_mono, by = c("Gene" = "gene")) %>%
  group_by(GO_ID, Description, Cluster,p.adjust,Count) %>%
  summarise(
    N_genes = n(),
    mean_abs_logFC = mean(abs(logFC), na.rm = TRUE),
    median_abs_logFC = median(abs(logFC), na.rm = TRUE)
  ) %>%
  arrange(desc(median_abs_logFC))
#write.csv(go_magnitude, "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/compare_14_84/GO_unique_genes/GO_term_mag_d14_minGsize20_sim0.35", row.names = FALSE)
top_terms <- path_magnitude %>% top_n(35, median_abs_logFC)


p<-ggplot(top_terms, aes(x = median_abs_logFC,
                      y = reorder(Description, median_abs_logFC),
                      color = p.adjust,
                      size = Count)) +
  geom_point() +
  facet_wrap(~Cluster) +
  scale_color_gradient(low = "red", high = "yellow", name = "FDR (p.adjust)") +
  scale_size_continuous(name = "Gene Count") +
  labs(
    x = "Median |log2FC|",
    y = "MsigDB term"
  ) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10, face = "bold"))
write.csv(top_terms, "/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/GSEA_GO/GO_mono/GO_terms/mono_PC_GO_path.csv", row.names = FALSE)
ggsave(file = "/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/GSEA_GO/GO_mono/GO_terms/mono_PC_GO_path.png", p, width = 11, height = 5, units = "in")

go_df <- as.data.frame(clust_simple)
go_df$FirstGene <- sapply(strsplit(as.character(go_df$geneID), "/"), function(x) {paste(head(x,1), collapse =",")})
top_terms <- go_df[order(go_df$p.adjust), ][1:min(20, nrow(go_df)), ]

plot <- ggplot(top_terms, aes(
  x = -log10(p.adjust),
  y = reorder(Description, -log10(p.adjust)),
  size = Count,
  color = -log10(p.adjust)
)) +
  geom_point() +
  geom_text_repel(aes(label = FirstGene),
                  size = 3,
                  max.overlaps = 20,
                  box.padding = 0.5,
                  segment.color = "grey50") +
  scale_color_gradient(low = "green", high = "darkgreen") +
  theme_minimal(base_size = 15) +theme(axis.text.y = element_text(face = "bold"),axis.text.x = element_text(face = "bold")) +
  labs(
    x = "-log10(FDR)",
    y = "GO Term",
    size = "Gene Count",
    color = "-log10(FDR)"
  )
plot
ggsave(file = "/work/ABG/mkapoor/mkapoor/PRRSV/PRRSV_cellranger_v97/GSEA_GO/GO_mono/GO_terms/mono_PC_GO_path_v2.png", plot, width = 12, height = 5, units = "in")




##creative 
library(plotly)

term_df <- as.data.frame(cclust) %>%
  dplyr::mutate(hover = paste0(Description, "<br>Genes: ", geneID)) %>%
  dplyr::arrange(p.adjust) %>% dplyr::slice_head(n = 20)

p <- ggplot(term_df, aes(-log10(p.adjust),
                         reorder(Description, -log10(p.adjust)),
                         size = Count, color = -log10(p.adjust),
                         text = hover)) +
  geom_point() + theme_minimal()
p

ggplotly(p, tooltip = "text")
