## Psuedobulk by  celltype to see if we can find biomarkers for prediction - script to run ###
#Using edgeR GLM
#using limma voom
#using monocyte signature residuals 
###############################################################################

.libPaths("/work/ABG/mkapoor/.ondemand-new/mkapoor/rstudio/libs/4.4.1")
library(Seurat)
library(tidyverse)
library(cowplot)
library(edgeR)
library(dplyr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(png)
library(RColorBrewer)
library(SingleCellExperiment)
library(SummarizedExperiment)
#library(MAST)
library(edgeR)
library(data.table)
library(limma)
library(variancePartition)
library(lmerTest)
library(EnhancedVolcano)
#load day 14
seurat_14dpi <- readRDS("/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Project_Fang_10X/filtered_postQC_postcb_postdoublet_postdowns_annotated_14dpi.rds")
Idents(seurat_14dpi) <- seurat_14dpi$celltypes

#We will separarte each pseudobulk profiles by its celltype##

#create sce + metadata #
counts <- seurat_14dpi@assays$RNA@counts 
metadata <- seurat_14dpi@meta.data
# Set up metadata as desired for aggregation and DE analysis
metadata$celltypes <- factor(seurat_14dpi@active.ident)

# Create single cell experiment object
sce <- SingleCellExperiment(assays = list(counts = counts), 
                            colData = metadata)

#remove lowly expressed with less than 1 counts
sce_2 <- sce[rowSums(counts(sce) > 1) >= 1, ]
dim(sce_2)

counts_14 <- counts(sce_2)
meta_14   <- as.data.frame(colData(sce_2))
stopifnot(ncol(counts_14) == nrow(meta_14))
stopifnot(identical(colnames(counts_14), rownames(meta_14)))
pbmc_counts14 <- t(rowsum(t(as.matrix(counts_14)), group = meta_14$Sample))
meta_pbmc14 <- meta_14 %>%
  dplyr::distinct(Sample, .keep_all = TRUE) %>%
  dplyr::select(Sample, Treatment, Sow)

##EP##
meta_pe <- meta_pbmc14 %>%
  dplyr::filter(Treatment %in% c("persistent","extinct")) %>%
  droplevels()
rownames(meta_pe) <- meta_pe$Sample
keep <- intersect(colnames(pbmc_counts14), rownames(meta_pe))
pbmc_counts14_pe <- pbmc_counts14[, keep, drop = FALSE]
meta_pe<- meta_pe[keep, , drop = FALSE]
stopifnot(identical(colnames(pbmc_counts14_pe), rownames(meta_pe)))

#cell type fraction
tab  <- with(metadata, table(Sample, celltypes))           # cells per Sample×celltype
frac <- prop.table(tab, 1) |> as.data.frame.matrix()       # row-wise fractions
frac <- frac[rownames(meta_pe), , drop = FALSE]            # align to meta_pe rows
pcs <- prcomp(frac, scale. = TRUE)$x
meta_pe$compPC1 <- pcs[, 1]
meta_pe$compPC2 <- if (ncol(pcs) >= 2) pcs[, 2] else 0
##try limma voom for duplicate correlation
y   <- DGEList(pbmc_counts14_pe)
keep <- filterByExpr(y, model.matrix(~ 0 + Treatment + compPC1 + compPC2, meta_pe))
y   <- y[keep,, keep.lib.sizes=FALSE]
y   <- calcNormFactors(y)
X   <- model.matrix(~ 0 + Treatment + compPC1 + compPC2, data = meta_pe)
colnames(X) <- sub("Treatment", "", colnames(X))
v   <- voom(y, X, plot = FALSE)
corfit <- duplicateCorrelation(v, design = X, block = meta_pe$Sow)
v   <- voom(y, X, plot = FALSE)  # re-voom
fit <- lmFit(v, X, block = meta_pe$Sow, correlation = corfit$consensus)
ct  <- makeContrasts(extinct - persistent, levels = colnames(X))
fit2 <- contrasts.fit(fit, ct)
fit2 <- eBayes(fit2, robust = TRUE)
res  <- topTable(fit2, number = Inf)
head(res)
sc_de <-read.csv("/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/biomarker_from_sc/candidate_DEG_PI_D14.csv")
#sc_de <-read.csv("/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Project_Fang_10X/MAST_DE_14dpi/CD8_NK/Mast_CD8_NK_ext_per_sigFDR.csv")
#sc_de <-na.omit(sc_de)
genes_sc_sig <- sc_de$gene[sc_de$fdr < 0.05 & sc_de$absLFC > 0.25]
pb <- as.data.frame(res)
pb$gene <- rownames(pb)

keep <- with(pb, (adj.P.Val < 0.25) | (abs(logFC) >= 0.05))
res_sub <- pb[keep, , drop = FALSE]
res_sub <- as.data.frame(res_sub)
res_sub$highlight <- ifelse(rownames(res_sub) %in% genes_sc_sig, "sc_DEG", "pb_DEG")
View(res_sub)
genes_use <- intersect(genes_sc_sig, rownames(res_sub))
sc_sub  <- sc_de[match(genes_use, sc_de$gene), ]
bulk_sub <- res[match(genes_use, rownames(res)), ]
compare_df <- data.frame(
  gene      = genes_use,
  logFC_sc  = sc_sub$logFC,
  logFC_pb  = bulk_sub$logFC,
  FDR_sc    = sc_sub$fdr,
  FDR_pb    = bulk_sub$adj.P.Val
)
write.csv(compare_df,
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X//biomarker_from_sc/discr_pb_sc_VE_D14.csv",
          row.names = FALSE)
cor_test <- cor.test(compare_df$logFC_sc, compare_df$logFC_pb, method = "spearman") 
cor_test
write.csv(res,
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/limma_voom/cd8_nk/limma_voom_results_pb_EP.csv",
          row.names = FALSE)
write.csv(res_sub,
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/limma_voom/cd8_nk/limma_voom_results_pb_FC_0.05_EP.csv",
          row.names = FALSE)
write.csv(compare_df,
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/limma_voom/cd8_nk/limma_voom_results_compare_pb_sc_EP_0.05.csv",
          row.names = FALSE)
topN <- 25

plot_df <- compare_df %>%
  mutate(abs_sc = abs(logFC_sc)) %>%
  arrange(desc(abs_sc)) %>%
  slice_head(n = topN) %>%
  mutate(gene = factor(gene, levels = rev(gene)))
plot_long <- plot_df %>%
  select(gene, logFC_sc, FDR_sc, logFC_pb, FDR_pb) %>%
  pivot_longer(
    cols = c(logFC_sc, logFC_pb),
    names_to = "source",
    values_to = "logFC"
  ) %>%
  mutate(
    source = recode(source,
                    logFC_sc = "scRNA-seq (Monocytes)",
                    logFC_pb = "Pseudobulk"),
    FDR = ifelse(source == "scRNA-seq (Monocytes)", FDR_sc, FDR_pb),
    neglog10FDR = -log10(pmax(FDR, 1e-300))
  )

p <- ggplot(plot_long, aes(x = logFC, y = gene)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(aes(color = neglog10FDR), size = 3) +
  facet_wrap(~ source, ncol = 2, scales = "free_x") +
  scale_color_gradient(
    low = "grey70",
    high = "#B2182B",   # deep red
    name = "-log10(FDR)"
  ) +
  labs(
    x = "log2FC",
    y = "Genes",
    color = "-log10(FDR)",
    title = "Discriminative Monocyte PI"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  theme(
    axis.title.x = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 12),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title.y = element_text(face = "bold", size = 12)
  )
p
ggsave("/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/biomarker_from_sc/cand_PI.png",p)
#do regression based 
y_pe   <- edgeR::DGEList(pbmc_counts14_pe)
y_pe   <- edgeR::calcNormFactors(y_pe)
lcpm_pe <- edgeR::cpm(y_pe, log = TRUE, prior.count = 1)
meta_pe$Treatment <- factor(meta_pe$Treatment, levels = c("persistent","extinct"))
dim(pbmc_counts14_pe); table(meta_pe$Treatment)
genes_sc_sig <- sc_de$primerid[sc_de$fdr < 0.05 & abs(sc_de$logFC) > 0.05]
genes_use <- intersect(genes_sc_sig, rownames(lcpm_pe))
pb_mat <- as.matrix(lcpm_pe[genes_use, , drop = FALSE])
##Per-sample composition from metdata
tab  <- with(metadata, table(Sample, celltypes))             # cells per Sample×celltype
frac <- prop.table(tab, 1) |> as.data.frame.matrix()
frac <- frac[colnames(lcpm_pe), , drop = FALSE]
stopifnot(identical(rownames(frac), colnames(lcpm_pe)))
pcs <- prcomp(frac, scale. = TRUE)$x
compPC1 <- pcs[, 1]; compPC2 <- pcs[, 2]
mono_sig <- colMeans(pb_mat, na.rm = TRUE)
df <- data.frame(
  Sample    = colnames(lcpm_pe),
  mono_sig  = as.numeric(mono_sig),
  Treatment = meta_pe[colnames(lcpm_pe), "Treatment", drop = TRUE],
  Sow       = meta_pe[colnames(lcpm_pe), "Sow", drop = TRUE],
  compPC1   = compPC1,
  compPC2   = compPC2,
  stringsAsFactors = FALSE
); df$Treatment <- factor(df$Treatment, levels = c("persistent","extinct")) ; df$Sow <- factor(df$Sow)

#adjuts 
adj_fit <- lm(mono_sig ~ compPC1 + compPC2, data = df)
df$mono_sig_adj <- resid(adj_fit)
fit_lmm <- lm(mono_sig_adj ~ Treatment + Sow, data = df)
print(summary(fit_lmm))
genes_detected <- genes_use
gene_expr_summary <- data.frame(
  Gene = genes_detected,
  Mean_logCPM = rowMeans(lcpm_pe[genes_detected, , drop = FALSE], na.rm = TRUE),
  SD_logCPM   = apply(lcpm_pe[genes_detected, , drop = FALSE], 1, sd, na.rm = TRUE)
)
write.csv(gene_expr_summary,
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/regression_based/cd8_nk/GE_pb_summary_EP_logFC0.05.csv",
          row.names = FALSE)
write.csv(lcpm_pe[genes_use, ],
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/regression_based/cd8_nk/GE_pb_sample_EP_logFC0.05.csv")

plot <- ggplot(df, aes(Treatment, mono_sig_adj, fill = Treatment)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 3, aes(shape = Sow)) +
  labs(y = "cd8_nk residuals")+
  theme_bw(base_size = 14)
ggsave("/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/regression_based/cd8_nk/plot_EP.png", plot)
range(gene_expr_summary$Mean_logCPM)
range(gene_expr_summary$SD_logCPM)

###EC####
###########
meta_pe <- meta_pbmc14 %>%
  dplyr::filter(Treatment %in% c("control","extinct")) %>%
  droplevels()
rownames(meta_pe) <- meta_pe$Sample
keep <- intersect(colnames(pbmc_counts14), rownames(meta_pe))
pbmc_counts14_pe <- pbmc_counts14[, keep, drop = FALSE]
meta_pe          <- meta_pe[keep, , drop = FALSE]
stopifnot(identical(colnames(pbmc_counts14_pe), rownames(meta_pe)))#

tab  <- with(metadata, table(Sample, celltypes))           # cells per Sample×celltype
frac <- prop.table(tab, 1) |> as.data.frame.matrix()       # row-wise fractions
frac <- frac[rownames(meta_pe), , drop = FALSE]            # align to meta_pe rows
pcs <- prcomp(frac, scale. = TRUE)$x
meta_pe$compPC1 <- pcs[, 1]
meta_pe$compPC2 <- if (ncol(pcs) >= 2) pcs[, 2] else 0

#regression
y_pe   <- edgeR::DGEList(pbmc_counts14_pe)
y_pe   <- edgeR::calcNormFactors(y_pe)
lcpm_pe <- edgeR::cpm(y_pe, log = TRUE, prior.count = 1)
meta_pe$Treatment <- factor(meta_pe$Treatment, levels = c("control","extinct"))
dim(pbmc_counts14_pe); table(meta_pe$Treatment)

sc_de <-read.csv("/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Project_Fang_10X/MAST_DE_14dpi/CD8_NK/Mast_CD8_NK_ext_con_sigFDR.csv")
sc_de <- na.omit(sc_de)
genes_sc_sig <- sc_de$primerid[sc_de$fdr < 0.05 & abs(sc_de$logFC) > 0.05]
genes_use <- intersect(genes_sc_sig, rownames(lcpm_pe))
pb_mat <- as.matrix(lcpm_pe[genes_use, , drop = FALSE])
length(lcpm_pe)
##Per-sample composition from metdata
tab  <- with(metadata, table(Sample, celltypes))             # cells per Sample×celltype
frac <- prop.table(tab, 1) |> as.data.frame.matrix()
frac <- frac[colnames(lcpm_pe), , drop = FALSE]
stopifnot(identical(rownames(frac), colnames(lcpm_pe)))
pcs <- prcomp(frac, scale. = TRUE)$x
compPC1 <- pcs[, 1]; compPC2 <- pcs[, 2]
mono_sig <- colMeans(pb_mat, na.rm = TRUE)

df <- data.frame(
  Sample    = colnames(lcpm_pe),
  mono_sig  = as.numeric(mono_sig),
  Treatment = meta_pe[colnames(lcpm_pe), "Treatment", drop = TRUE],
  Sow       = meta_pe[colnames(lcpm_pe), "Sow", drop = TRUE],
  compPC1   = compPC1,
  compPC2   = compPC2,
  stringsAsFactors = FALSE
); df$Treatment <- factor(df$Treatment, levels = c("control","extinct")) ; df$Sow <- factor(df$Sow)

#adjuts 
adj_fit <- lm(mono_sig ~ compPC1 + compPC2, data = df)
df$mono_sig_adj <- resid(adj_fit)
fit_lmm <- lm(mono_sig_adj ~ Treatment , data = df)
print(summary(fit_lmm))
plot <- ggplot(df, aes(Treatment, mono_sig_adj, fill = Treatment)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 3, aes(shape = Sow)) +
  labs(y = "cd8_nk residuals") +
  theme_bw(base_size = 14)

genes_detected <- genes_use
gene_expr_summary <- data.frame(
  Gene = genes_detected,
  Mean_logCPM = rowMeans(lcpm_pe[genes_detected, , drop = FALSE], na.rm = TRUE),
  SD_logCPM   = apply(lcpm_pe[genes_detected, , drop = FALSE], 1, sd, na.rm = TRUE)
)
write.csv(gene_expr_summary ,
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/regression_based/cd8_nk/GE_pb_summary_EC_logFC0.05.csv",
          row.names = FALSE)
write.csv(lcpm_pe[genes_use, ],
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/regression_based/cd8_nk/GE_pb_sample_EC_logFC0.05.csv")
ggsave("/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/regression_based/cd8_nk/plot_EC.png", plot)
range(gene_expr_summary$Mean_logCPM)
range(gene_expr_summary$SD_logCPM)

###PC####
###########
meta_pe <- meta_pbmc14 %>%
  dplyr::filter(Treatment %in% c("control","persistent")) %>%
  droplevels()
rownames(meta_pe) <- meta_pe$Sample
keep <- intersect(colnames(pbmc_counts14), rownames(meta_pe))
pbmc_counts14_pe <- pbmc_counts14[, keep, drop = FALSE]
meta_pe          <- meta_pe[keep, , drop = FALSE]
stopifnot(identical(colnames(pbmc_counts14_pe), rownames(meta_pe)))

tab  <- with(metadata, table(Sample, celltypes))           # cells per Sample×celltype
frac <- prop.table(tab, 1) |> as.data.frame.matrix()       # row-wise fractions
frac <- frac[rownames(meta_pe), , drop = FALSE]            # align to meta_pe rows
pcs <- prcomp(frac, scale. = TRUE)$x
meta_pe$compPC1 <- pcs[, 1]
meta_pe$compPC2 <- if (ncol(pcs) >= 2) pcs[, 2] else 0
stopifnot(identical(rownames(frac), rownames(meta_pe)))

##try limma voom for duplicate correlation
y   <- DGEList(pbmc_counts14_pe)
keep <- filterByExpr(y, model.matrix(~ 0 + Treatment + compPC1 + compPC2, meta_pe))
y   <- y[keep,, keep.lib.sizes=FALSE]
y   <- calcNormFactors(y)
X   <- model.matrix(~ 0 + Treatment + compPC1 + compPC2, data = meta_pe)
colnames(X) <- sub("Treatment", "", colnames(X))
v   <- voom(y, X, plot = FALSE)
corfit <- duplicateCorrelation(v, design = X, block = meta_pe$Sow)
v   <- voom(y, X, plot = FALSE)  # re-voom
fit <- lmFit(v, X, block = meta_pe$Sow, correlation = corfit$consensus)
ct  <- makeContrasts( persistent- control, levels = colnames(X))
fit2 <- contrasts.fit(fit, ct)
fit2 <- eBayes(fit2, robust = TRUE)
res  <- topTable(fit2, number = Inf)
head(res)
sc_de <-read.csv("/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Project_Fang_10X/MAST_DE_14dpi/CD8_NK/Mast_CD8_NK_per_con_sigFDR.csv")
sc_de <- na.omit(sc_de)
genes_sc_sig <- sc_de$primerid[sc_de$fdr < 0.05 & abs(sc_de$logFC) > 0.05]

keep <- with(res, (adj.P.Val < 0.25) | (abs(logFC) >= 0.05))
res_sub <- res[keep, , drop = FALSE]
res_sub <- as.data.frame(res_sub)
res_sub$highlight <- ifelse(rownames(res_sub) %in% genes_sc_sig, "sc_DEG", "pb_DEG")

genes_use <- intersect(genes_sc_sig, rownames(res_sub))
sc_sub  <- sc_de[match(genes_use, sc_de$primerid), ]
bulk_sub <- res[match(genes_use, rownames(res)), ]
compare_df <- data.frame(
  gene      = genes_use,
  logFC_sc  = sc_sub$logFC,
  logFC_pb  = bulk_sub$logFC,
  FDR_sc    = sc_sub$fdr,
  FDR_pb    = bulk_sub$adj.P.Val
)
cor_test <- cor.test(compare_df$logFC_sc, compare_df$logFC_pb, method = "spearman") 
cor_test
write.csv(res,     file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/limma_voom/cd8_nk/limma_voom_results_pb_PC.csv",
          row.names = FALSE)
write.csv(res_sub,
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/limma_voom/cd8_nk/limma_voom_results_pb_FC_0.05_PC.csv",
          row.names = FALSE)
write.csv(compare_df,
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/limma_voom/cd8_nk/limma_voom_results_compare_pb_sc_PC_0.05.csv",
          row.names = FALSE)


#regression
y_pe   <- edgeR::DGEList(pbmc_counts14_pe)
y_pe   <- edgeR::calcNormFactors(y_pe)
lcpm_pe <- edgeR::cpm(y_pe, log = TRUE, prior.count = 1)
meta_pe$Treatment <- factor(meta_pe$Treatment, levels = c("control","persistent"))
dim(pbmc_counts14_pe); table(meta_pe$Treatment)

#check if monocyte DE show same as pseudobulk DE on MDS - just EP
genes_sc_sig <- sc_de$primerid[sc_de$fdr < 0.05 & abs(sc_de$logFC) > 0.05]
genes_use <- intersect(genes_sc_sig, rownames(lcpm_pe))
pb_mat <- as.matrix(lcpm_pe[genes_use, , drop = FALSE])
length(lcpm_pe)
##Per-sample composition from metdata
tab  <- with(metadata, table(Sample, celltypes))             # cells per Sample×celltype
frac <- prop.table(tab, 1) |> as.data.frame.matrix()
frac <- frac[colnames(lcpm_pe), , drop = FALSE]
stopifnot(identical(rownames(frac), colnames(lcpm_pe)))
pcs <- prcomp(frac, scale. = TRUE)$x
compPC1 <- pcs[, 1]; compPC2 <- pcs[, 2]
mono_sig <- colMeans(pb_mat, na.rm = TRUE)

df <- data.frame(
  Sample    = colnames(lcpm_pe),
  mono_sig  = as.numeric(mono_sig),
  Treatment = meta_pe[colnames(lcpm_pe), "Treatment", drop = TRUE],
  Sow       = meta_pe[colnames(lcpm_pe), "Sow", drop = TRUE],
  compPC1   = compPC1,
  compPC2   = compPC2,
  stringsAsFactors = FALSE
); df$Treatment <- factor(df$Treatment, levels = c("control","persistent")) ; df$Sow <- factor(df$Sow)

#adjuts 
adj_fit <- lm(mono_sig ~ compPC1 + compPC2, data = df)
df$mono_sig_adj <- resid(adj_fit)
fit_lmm <- lm(mono_sig_adj ~ Treatment+Sow , data = df)
print(summary(fit_lmm))
plot <- ggplot(df, aes(Treatment, mono_sig_adj, fill = Treatment)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 3, aes(shape = Sow)) +
  labs(y = "cd8_nk residuals")+
  theme_bw(base_size = 14)

ggsave("/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/regression_based/cd8_nk/plot_PC.png", plot)

genes_detected <- genes_use
gene_expr_summary <- data.frame(
  Gene = genes_detected,
  Mean_logCPM = rowMeans(lcpm_pe[genes_detected, , drop = FALSE], na.rm = TRUE),
  SD_logCPM   = apply(lcpm_pe[genes_detected, , drop = FALSE], 1, sd, na.rm = TRUE)
)
write.csv(gene_expr_summary,
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/regression_based/cd8_nk/GE_pb_summary_PC_logFC0.05.csv",
          row.names = FALSE)
write.csv(lcpm_pe[genes_use, ],
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/regression_based/cd8_nk/GE_pb_sample_PC_logFC0.05.csv")
range(gene_expr_summary$Mean_logCPM)
range(gene_expr_summary$SD_logCPM)
