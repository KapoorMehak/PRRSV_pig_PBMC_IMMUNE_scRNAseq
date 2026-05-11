library(Seurat)
library(MAST)
library(data.table)
library(dplyr)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(stringr)
library(ggplot2)
library(lme4)
###DEfor AB subset analysis at res 0.5
seurat_84dpi_subset_1 <- readRDS("./PRRSV/filtered_postQC_postcb_postdoublet_postdowns_postcellcycle_ABT_subset_84dpi.rds")

exprs_data <- GetAssayData(seurat_84dpi_subset_1, layer = "data")
metadata <- seurat_84dpi_subset_1@meta.data
# Create a SingleCellAssay object for MAST
sca <- FromMatrix(exprsArray = as.matrix(exprs_data), 
                  cData = metadata)
sca <- SceToSingleCellAssay(sca, class = "SingleCellAssay")
print(dim(sca))
#expressed in atleast 5% of cells
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
  "cytotoxic_cd8a"       = sca_per_ext_list[["Cytotoxic CD8a+ T cells"]],
  "Activated_CD4CD8a_pos"            = sca_per_ext_list[["Activated CD4+CD8a+ ab T cells"]],
  "Naïve_CD4_CD8a_neg"            = sca_per_ext_list[["Naïve CD4+CD8a- ab T cells"]],
#  "Mixed_CD4_CD8a_pos"        = sca_per_ext_list[["Mixed CD4+CD8a+ ab T"]],
  "CD2_pos_GD"  = sca_per_ext_list[["CD2+ GD T cells"]]
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
  





