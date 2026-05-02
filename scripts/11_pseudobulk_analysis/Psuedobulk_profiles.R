## Psuedobulk by  celltype to see if we can find biomarkers for prediction ###
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
meta_pe <- meta_pbmc14 %>%
  dplyr::filter(Treatment %in% c("persistent","extinct")) %>%
  droplevels()
rownames(meta_pe) <- meta_pe$Sample
keep <- intersect(colnames(pbmc_counts14), rownames(meta_pe))
pbmc_counts14_pe <- pbmc_counts14[, keep, drop = FALSE]
meta_pe<- meta_pe[keep, , drop = FALSE]
stopifnot(identical(colnames(pbmc_counts14_pe), rownames(meta_pe)))

#Run EdgeR pipleine 
tab  <- with(metadata, table(Sample, celltypes))           # cells per Sample×celltype
frac <- prop.table(tab, 1) |> as.data.frame.matrix()       # row-wise fractions
frac <- frac[rownames(meta_pe), , drop = FALSE]            # align to meta_pe rows
pcs <- prcomp(frac, scale. = TRUE)$x
meta_pe$compPC1 <- pcs[, 1]
meta_pe$compPC2 <- if (ncol(pcs) >= 2) pcs[, 2] else 0

stopifnot(identical(rownames(frac), rownames(meta_pe)))
y_pe   <- edgeR::DGEList(pbmc_counts14_pe)
y_pe   <- edgeR::calcNormFactors(y_pe)
design <- model.matrix(~ 0 + Treatment +compPC1 +compPC2, data = meta_pe) #dropped sow because no df of to estimate doispersion. 
colnames(design) <- sub("Treatment", "", colnames(design))
design
keep <- filterByExpr(y_pe, design) #filter low expressed 
y_pe <- y_pe[keep, , keep.lib.sizes = FALSE]
y_pe <- estimateDisp(y_pe, design)
fit <- glmQLFit(y_pe, design)
qlf <- glmQLFTest(fit, contrast = makeContrasts(persistent - extinct, levels = design))
res <- topTags(qlf, n = Inf)$table
head(res)
summary(decideTests(qlf))
plotMD(qlf, column=1, status=decideTests(qlf)[,1],
       main="edgeR pseudobulk (extinct vs persistent)", xlab="Average logCPM")

EnhancedVolcano(res,
                lab = rownames(res),
                x   = 'logFC',
                y   = 'FDR',
                pCutoff = 0.05,
                FCcutoff = 1,
                title = 'Pseudobulk DE (extinct vs persistent)\nwith composition adjusted')

# Highlight overlap with scDEGs
genes_sc_sig <- sc_de$primerid[sc_de$fdr < 0.05 & abs(sc_de$logFC) > 0.25]
res$highlight <- ifelse(rownames(res) %in% genes_sc_sig, "sc_DEG", "other")
colCustom <- rep("grey60", nrow(res)); names(colCustom) <- rownames(res)
colCustom[rownames(res) %in% genes_sc_sig] <- "royalblue"

# EnhancedVolcano(res,
#                 lab = ifelse(res$highlight == "sc_DEG", rownames(res), ""),
#                 x   = 'logFC',
#                 y   = 'FDR',
#                 selectLab = rownames(res)[res$highlight == "sc_DEG"],
#                 pCutoff = 0.05, FCcutoff = 1,
#                 title = 'Overlap of sc and pseudobulk DE genes\n(composition adjusted)',
#                 subtitle = 'Red = sig in pseudobulk; Blue = sc DEGs',
#                 colCustom = colCustom, labSize = 3)

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
plotSA(fit, main="Residual stdev vs abundance") # variance trend 
corfit$consensus  #corr of RE
sc_de <-read.csv("/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Project_Fang_10X/MAST_DE_14dpi/B_cells/Mast_Bcell_ext_per_sigFDR.csv")
genes_sc_sig <- sc_de$primerid[sc_de$fdr < 0.05 & abs(sc_de$logFC) > 0.05]

keep <- with(res, (adj.P.Val < 0.25) | (abs(logFC) >= 0.05))
res_sub <- res[keep, , drop = FALSE]
res_sub <- as.data.frame(res_sub)
res_sub$highlight <- ifelse(rownames(res_sub) %in% genes_sc_sig, "sc_DEG", "pb_DEG")
ggplot(res_sub, aes(x = logFC, y = adj.P.Val, color = highlight)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("other" = "grey70", "sc_DEG" = "royalblue")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", linewidth = 0.3) +
  labs(
    x = "log2 Fold Change",
    y = "FDR",
    color = "Highlight" # 816 our of 1093 DEG from sc detected in pb when 0.05 logfc , 15 when 0.25 logfc but no signifcance until 0.4 FDR
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
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
cor_test <- cor.test(compare_df$logFC_sc, compare_df$logFC_pb, method = "spearman") #0.65 p-value < 2.2e-16 correlation of 816 DEGs , 0.95 of 15 DEGs at 0.25 logfc
write.csv(res,
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/limma_voom/B/limma_voom_results_pb_EP.csv",
          row.names = FALSE)
write.csv(res_sub,
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/limma_voom/B/limma_voom_results_pb_FC_0.05_EP.csv",
          row.names = FALSE)
write.csv(compare_df,
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/limma_voom/B/limma_voom_results_compare_pb_sc_EP_0.05.csv",
          row.names = FALSE)

#do regression based 
lcpm_pe <- edgeR::cpm(y_pe, log = TRUE, prior.count = 1)
meta_pe$Treatment <- factor(meta_pe$Treatment, levels = c("persistent","extinct"))
dim(pbmc_counts14_pe); table(meta_pe$Treatment)

#check if monocyte DE show same as pseudobulk DE on MDS - just EP
sc_de <-read.csv("/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Project_Fang_10X/MAST_DE_14dpi/Monocytes/Mast_Mono_ext_per_sigFDR.csv")
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
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/regression_based/B/GE_pb_summary_EP_logFC0.05.csv",
          row.names = FALSE)
write.csv(lcpm_pe[genes_use, ],
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/regression_based/B/GE_pb_sample_EP_logFC0.05.csv")

plot <- ggplot(df, aes(Treatment, mono_sig_adj, fill = Treatment)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 3, aes(shape = Sow)) +
  labs(y = "B residuals")+
  theme_bw(base_size = 14)
ggsave("/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/regression_based/B/plot.png", plot)
#doesn not yeild significnace, no detectuin after adjusting for cell composition
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
#Run EdgeR pipleine 
tab  <- with(metadata, table(Sample, celltypes))           # cells per Sample×celltype
frac <- prop.table(tab, 1) |> as.data.frame.matrix()       # row-wise fractions
frac <- frac[rownames(meta_pe), , drop = FALSE]            # align to meta_pe rows
pcs <- prcomp(frac, scale. = TRUE)$x
meta_pe$compPC1 <- pcs[, 1]
meta_pe$compPC2 <- if (ncol(pcs) >= 2) pcs[, 2] else 0

stopifnot(identical(rownames(frac), rownames(meta_pe)))
y_pe   <- edgeR::DGEList(pbmc_counts14_pe)
y_pe   <- edgeR::calcNormFactors(y_pe)
design <- model.matrix(~ 0 + Treatment +compPC1 +compPC2, data = meta_pe) #dropped sow because no df of to estimate doispersion. 
colnames(design) <- sub("Treatment", "", colnames(design))
design
keep <- filterByExpr(y_pe, design) #filter low expressed 
y_pe <- y_pe[keep, , keep.lib.sizes = FALSE]
y_pe <- estimateDisp(y_pe, design)
fit <- glmQLFit(y_pe, design)
qlf <- glmQLFTest(fit, contrast = makeContrasts(extinct - control, levels = design))
res <- topTags(qlf, n = Inf)$table
head(res)
summary(decideTests(qlf))

##try limma voom for duplicate correlation - no df for EC
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
ct  <- makeContrasts(extinct - control, levels = colnames(X))
fit2 <- contrasts.fit(fit, ct)
fit2 <- eBayes(fit2, robust = TRUE)
res  <- topTable(fit2, number = Inf)
head(res)
plotSA(fit, main="Residual stdev vs abundance") # variance trend 
corfit$consensus  #corr of RE

y_pe   <- edgeR::DGEList(pbmc_counts14_pe)
y_pe   <- edgeR::calcNormFactors(y_pe)
lcpm_pe <- edgeR::cpm(y_pe, log = TRUE, prior.count = 1)
meta_pe$Treatment <- factor(meta_pe$Treatment, levels = c("control","extinct"))
dim(pbmc_counts14_pe); table(meta_pe$Treatment)

#check if monocyte DE show same as pseudobulk DE on MDS - just EC
sc_de <-read.csv("/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Project_Fang_10X/MAST_DE_14dpi/B_cells/Mast_B_cells_ext_con_sigFDR.csv")
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
  labs(y = "B residuals") +
  theme_bw(base_size = 14)

genes_detected <- genes_use
gene_expr_summary <- data.frame(
  Gene = genes_detected,
  Mean_logCPM = rowMeans(lcpm_pe[genes_detected, , drop = FALSE], na.rm = TRUE),
  SD_logCPM   = apply(lcpm_pe[genes_detected, , drop = FALSE], 1, sd, na.rm = TRUE)
)
write.csv(gene_expr_summary ,
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/regression_based/B/GE_pb_summary_EC_logFC0.05.csv",
          row.names = FALSE)
write.csv(lcpm_pe[genes_use, ],
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/regression_based/B/GE_pb_sample_EC_logFC0.05.csv")
ggsave("/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/regression_based/B/plot_EC.png", plot)
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
##Run EdgeR pipleine 
tab  <- with(metadata, table(Sample, celltypes))           # cells per Sample×celltype
frac <- prop.table(tab, 1) |> as.data.frame.matrix()       # row-wise fractions
frac <- frac[rownames(meta_pe), , drop = FALSE]            # align to meta_pe rows
pcs <- prcomp(frac, scale. = TRUE)$x
meta_pe$compPC1 <- pcs[, 1]
meta_pe$compPC2 <- if (ncol(pcs) >= 2) pcs[, 2] else 0

stopifnot(identical(rownames(frac), rownames(meta_pe)))
y_pe   <- edgeR::DGEList(pbmc_counts14_pe)
y_pe   <- edgeR::calcNormFactors(y_pe)
design <- model.matrix(~ 0 + Treatment+compPC1 +compPC2, data = meta_pe) #dropped sow because no df of to estimate doispersion. 
colnames(design) <- sub("Treatment", "", colnames(design))
design
keep <- filterByExpr(y_pe, design) #filter low expressed 
y_pe <- y_pe[keep, , keep.lib.sizes = FALSE]
y_pe <- estimateDisp(y_pe, design)
fit <- glmQLFit(y_pe, design)
qlf <- glmQLFTest(fit, contrast = makeContrasts(persistent - control, levels = design))
res <- topTags(qlf, n = Inf)$table
head(res)
summary(decideTests(qlf))

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
plotSA(fit, main="Residual stdev vs abundance") # variance trend 
corfit$consensus  #corr of RE
genes_sc_sig <- sc_de$primerid[sc_de$fdr < 0.05 & abs(sc_de$logFC) > 0.25]

keep <- with(res, (adj.P.Val < 0.25) | (abs(logFC) >= 0.25))
res_sub <- res[keep, , drop = FALSE]
res_sub <- as.data.frame(res_sub)
res_sub$highlight <- ifelse(rownames(res_sub) %in% genes_sc_sig, "sc_DEG", "pb_DEG")
ggplot(res_sub, aes(x = logFC, y = adj.P.Val, color = highlight)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("other" = "grey70", "sc_DEG" = "royalblue")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.3) +
  labs(
    x = "log2 Fold Change",
    y = "FDR",
    color = "Highlight" # 273 our of 1542 DEG from sc detected in pb when 0.25 logfc , 2704 out of 3911 when 0.25 logfc but no signifcance until 0.4 FDR
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
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
cor_test <- cor.test(compare_df$logFC_sc, compare_df$logFC_pb, method = "spearman") #0.29 p-value < .88e-06 correlation of 273 DEGs , 0.29 of 2704 DEGs at 0.25 logfc
write.csv(res,
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/limma_voom/B/limma_voom_results_pb_PC.csv",
          row.names = FALSE)
write.csv(res_sub,
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/limma_voom/B/limma_voom_results_pb_FC_0.25_PC.csv",
          row.names = FALSE)
write.csv(compare_df,
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/limma_voom/B/limma_voom_results_compare_pb_sc_PC_0.25.csv",
          row.names = FALSE)
ggsave("/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/regression_based/B/plot_PC.png", plot)
range(gene_expr_summary$Mean_logCPM)
range(gene_expr_summary$SD_logCPM)

#regression
y_pe   <- edgeR::DGEList(pbmc_counts14_pe)
y_pe   <- edgeR::calcNormFactors(y_pe)
lcpm_pe <- edgeR::cpm(y_pe, log = TRUE, prior.count = 1)
meta_pe$Treatment <- factor(meta_pe$Treatment, levels = c("control","persistent"))
dim(pbmc_counts14_pe); table(meta_pe$Treatment)

#check if monocyte DE show same as pseudobulk DE on MDS - just EP
sc_de <-read.csv("/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Project_Fang_10X/MAST_DE_14dpi/B_cells/Mast_Bcell_per_con_sigFDR.csv")
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
  labs(y = "Monocyte residuals")+
  theme_bw(base_size = 14)


genes_detected <- genes_use
gene_expr_summary <- data.frame(
  Gene = genes_detected,
  Mean_logCPM = rowMeans(lcpm_pe[genes_detected, , drop = FALSE], na.rm = TRUE),
  SD_logCPM   = apply(lcpm_pe[genes_detected, , drop = FALSE], 1, sd, na.rm = TRUE)
)
write.csv(gene_expr_summary,
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/regression_based/B/GE_pb_summary_PC_logFC0.05.csv",
          row.names = FALSE)
write.csv(lcpm_pe[genes_use, ],
          file = "/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/regression_based/B/GE_pb_sample_PC_logFC0.05.csv")

#plot - meanlogCPM= abg abundance of gene measured; sdlogcpm = noise/variability, cv = sd/mean, auc = predcitove power for evry gene across all samples
#test roc if it sepeartes phenotypes; if fail NA
# auc to separate ability of predicting groups
#ci_low, ci_hihg - show uncertainity

genes <-read.csv("/work/ABG/mkapoor/mkapoor/Project_Fang_10X/14dpi_PRRSV/Project_Fang_10X/Psuedobulk/GE_pb_summary_PC_logFC0.05.csv")
gene_expr_summary$CV_log <- with(gene_expr_summary, SD_logCPM / (Mean_logCPM + 1e-8))
meta_pe$Treatment <- factor(meta_pe$Treatment, levels = c("persistent", "extinct"))
auc_data <- lapply(gene_expr_summary$Gene, function(g) {
  expr <- as.numeric(lcpm_pe[g, ])
  ro <- tryCatch(
    roc(meta_pe$Treatment, expr, levels = c("persistent", "extinct")),
    error = function(e) NULL
  )
  if (!is.null(ro)) {
    ci_vals <- ci.auc(ro, conf.level = 0.95)
    data.frame(
      Gene = g,
      AUC  = as.numeric(auc(ro)),
      CI_low = ci_vals[1],
      CI_high = ci_vals[3]
    )
  } else {
    data.frame(Gene = g, AUC = NA, CI_low = NA, CI_high = NA)
  }
})
#combine by gene
auc_df <- bind_rows(auc_data)
gene_expr_auc <- gene_expr_summary %>%
  left_join(auc_df, by = "Gene")
gene_expr_auc <- gene_expr_auc %>%
  mutate(CI_width = CI_high - CI_low) %>%
  arrange(desc(AUC)) #width to check narrowness of uncertainity

best_genes <- gene_expr_auc %>%
  filter(!is.na(AUC),
         Mean_logCPM > 6,      # detectable
         CV_log      < 0.30,   # stable across sample
         AUC         > 0.80,   # predictive string separation
         CI_low      >0.6)  %>% #95%CI ; if rnadom >0.5
  arrange(desc(AUC), CI_width, CV_log)

print(head(best_genes, 20), row.names = FALSE)


boxplot(as.numeric(lcpm_pe["PLAC8", ]) ~ meta_pe$Treatment,
        main = paste("PLAC8 AUC =", round(auc(roc(meta_pe$Treatment, lcpm_pe["PLAC8", ])), 3)))

boxplot(as.numeric(lcpm_pe["LUC7L3", ]) ~ meta_pe$Treatment,
        main = paste("LUC7L3 AUC =", round(auc(roc(meta_pe$Treatment, lcpm_pe["LUC7L3", ])), 3)))

boxplot(as.numeric(lcpm_pe["SRGN", ]) ~ meta_pe$Treatment,
        main = paste("SRGN AUC =", round(auc(roc(meta_pe$Treatment, lcpm_pe["SRGN", ])), 3)))

boxplot(as.numeric(lcpm_pe["ATP6", ]) ~ meta_pe$Treatment,
        main = paste("ATP6 AUC =", round(auc(roc(meta_pe$Treatment, lcpm_pe["ATP6", ])), 3)))
range(lcpm_pe["ATP6", ])
