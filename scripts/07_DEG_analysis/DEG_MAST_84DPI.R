library(Seurat)
library(MAST)
library(data.table)
library(dplyr)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(stringr)
library(ggplot2)
library(lme4)
## DE for day 84 ###
seurat_84dpi <- readRDS("./PRRSV/filtered_postQC_postcb_postdoublet_postdowns_postcellcycle_updated_annotation_84dpi.rds")
# get cell number per celltype per Sow per Sex per Treatment
tab <- table(seurat_84dpi$Sex, 
            seurat_84dpi$Sow, 
            seurat_84dpi$Treatment, 
            seurat_84dpi$CellTypes)
tab_df <- as.data.frame(tab)
colnames(tab_df) <- c("Sex", "Sow", "Treatment", "CellType", "Count")
print(tab_df)
ggplot(tab_df, aes(x = Sow, y = Count, fill = Sex)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Treatment ~ CellType) +
  theme_minimal() + RotatedAxis()
  labs(title = "Cell counts by Sow, Treatment, and Cell Type")

exprs_data <- GetAssayData(seurat_84dpi, layer = "data")
metadata <- seurat_84dpi@meta.data
# Create a SingleCellAssay object for MAST
sca <- FromMatrix(exprsArray = as.matrix(exprs_data), 
                  cData = metadata)
sca <- SceToSingleCellAssay(sca, class = "SingleCellAssay")
print(dim(sca))
sca <- sca[freq(sca) > 0.05, ]
print(dim(sca))
sca_per_ext <- sca[, sca$Treatment %in% c("persistent", "extinct")]
sca_per_ext$Treatment <- droplevels(sca_per_ext$Treatment)
#Split by cell type
cell_types <- metadata$CellTypes
# Split the SCE object by cCellTypes# Split the SCE object by cell sca_per_ext
split_sca_per_ext <- split(seq_len(ncol(sca_per_ext)), cell_types)
# Create standalone SCE objects for each cell type
sca_per_ext_list <- lapply(split_sca_per_ext, function(indices) {
  sca_per_ext[, indices]
})
celltype_list <- list(
  "B_cells"         = sca_per_ext_list[["B cells"]],
  "ASCs"            = sca_per_ext_list[["ASC"]],
  "CD2_neg_T"       = sca_per_ext_list[["CD2- GD T cells"]],
  "CD2_pos_T"       = sca_per_ext_list[["CD2+ GD T cells"]],
  "Monocytes"       = sca_per_ext_list[["Monocytes"]],
  "pDCs"            = sca_per_ext_list[["pDCs"]],
  "cDCs"            = sca_per_ext_list[["cDCs"]],
  "NK_cells"        = sca_per_ext_list[["NK cells"]],
  "Transitional_B"  = sca_per_ext_list[["Transitional B-like"]]
)

for (ct_name in names(celltype_list)) {
  sca_obj <- celltype_list[[ct_name]]
  print(paste0("Processing: ", ct_name))
  
  cdr <- colSums(assay(sca_obj) > 0)
  colData(sca_obj)$cdr <- scale(cdr)
  
  Treatment <- factor(colData(sca_obj)$Treatment)
  Treatment <- relevel(Treatment, "persistent")
  colData(sca_obj)$Treatment <- Treatment
  
  zlmCond <- zlm(formula = ~ Treatment + cdr + Sex + (1 | Sow),
                 sca = sca_obj,
                 method = 'glmer',
                 ebayes = FALSE,
                 strictConvergence = FALSE,
                 fitArgsD = list(nAGQ = 0))
  
  summaryCond <- summary(zlmCond, doLRT = 'Treatmentextinct')
  summaryDt <- summaryCond$datatable
  
  fcHurdle <- merge(
    summaryDt[contrast == 'Treatmentextinct' & component == 'H', .(primerid, `Pr(>Chisq)`)],
    summaryDt[contrast == 'Treatmentextinct' & component == 'logFC', .(primerid, logFC = coef, ci.hi, ci.lo)],
    by = 'primerid'
  )
  fcHurdle[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdleSig <- merge(fcHurdle, as.data.table(mcols(sca_obj)), by = 'primerid')
  setorder(fcHurdleSig, fdr)
  fcHurdleSig_DE <- fcHurdleSig[fdr < 0.05]
  cat(paste0("# of genes with FDR < 0.05 in ", ct_name, " Extinct vs Persistent: ", nrow(fcHurdleSig_DE), "\n"))
  
  out_dir <- "./PRRSV/MAST_DE_84dpi"
  out_path <- file.path(out_dir, ct_name, paste0("Mast_", ct_name, "_ext_per.csv"))
  dir.create(file.path(out_dir, ct_name), showWarnings = FALSE, recursive = TRUE)
  write.csv(fcHurdleSig, file = out_path)
}

#cannot find DE genes for cDCs
