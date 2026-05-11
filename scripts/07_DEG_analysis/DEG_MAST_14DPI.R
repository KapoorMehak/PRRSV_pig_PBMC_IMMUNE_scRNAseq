### DE for all celltypes in treatments ####
# load libraries #
library(Seurat)
library(miloR)
library(kableExtra)
#library(SeuratDisk)
library(zellkonverter)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
#library(scRNAseq)
library(scuttle)
library(irlba)
library(BiocParallel)
library(ggplot2)
library(tidyr)
library(hdf5r)
library(HDF5Array)
library(scCustomize)
library(tidyr)
library(dplyr)
library(ggplot2)
library(miloDE)
library(edgeR)
library(qvalue)
# analysis libraries
library(scuttle)
suppressMessages(library(miloR))
suppressMessages(library(uwot))
library(scran)
suppressMessages(library(dplyr))
library(reshape2)
library(tidyverse)
library(cowplot)
library(Matrix)
library(ComplexHeatmap)
library(ggplot2)
library(GGally)
library(limma)
library(reshape2)
library(data.table)
library(knitr)
library(stringr)
library(NMF)
library(rsvd)
library(RColorBrewer)
library(MAST)



## DE for day 14 ###
seurat_14dpi <- readRDS("./PRRSV/filtered_postQC_postcb_postdoublet_postdowns_14dpi.rds")
Idents(seurat_14dpi) <-  seurat_14dpi$louvain_res1_5
levels(seurat_14dpi) <- c('10','14', #monocytes
                          '29', #pDCs
                          '1','2','5','21','23', #B
                          '25', #antibody-secreting cells
                          '6','7','16','28', #CD4 T
                          '3','13', #CD8AB T
                          '8','9','12','19', #innate CD8/NK
                          '11','18','20', #NK
                          '22', #CD2+ GD T
                          '0','4','15','17','24','26','27') #CD2- GD T) 
seurat_14dpi$neworder <- Idents(seurat_14dpi) # Reorder the clusters based on putative cell type IDs we came up with from looking at the data
Idents(seurat_14dpi) <- seurat_14dpi$neworder
Mono <- rep('Monocytes', 2)
pDC <- 'pDCs'
B <- rep('B cells', 5)
ASC <- 'ASC'
CD4T <- rep('CD4+ ab T cells', 4)
CD8T <- rep('CD8ab+ ab T cells', 2)
CD8TNK <- rep('CD8a+ ab T/NK cells', 4)
NK <- rep('NK cells',3)
CD2posGD <- 'CD2+ GD T cells'
CD2negGD <- rep('CD2- GD T cells', 7)

CellTypes <- c(Mono, pDC, B, ASC, CD4T, CD8T, CD8TNK, NK, CD2posGD, CD2negGD)

seurat_14dpi$celltypes <-seurat_14dpi$neworder
Idents(seurat_14dpi) <- seurat_14dpi$celltypes
names(CellTypes) <- levels(seurat_14dpi) # assign CellTypes to cluster numbers
seurat_14dpi <- RenameIdents(seurat_14dpi, CellTypes) # change dataset identity to cell types in Seurat object
seurat_14dpi$celltypes <- Idents(seurat_14dpi)
Idents(seurat_14dpi) <- seurat_14dpi$neworder
Idents(seurat_14dpi) <- seurat_14dpi$celltypes
cols <- c('lightpink', 'red', 'orange', 'gold', 'darkgreen', 'mediumseagreen', 'skyblue2', 'steelblue', 'navy', 'plum3', 'darkmagenta', 'black', 'grey')
DimPlot(seurat_14dpi, cols = cols)
Idents(seurat_14dpi) <- seurat_14dpi$louvain_res1_5
DimPlot(seurat_14dpi,label = TRUE)
Idents(seurat_14dpi) %>% head()



exprs_data <- GetAssayData(seurat_14dpi, layer = "data")
metadata <- seurat_14dpi@meta.data
library(MAST)
freq_expressed <- 0.2
FCTHRESHOLD <- log2(1.5)
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
cell_types <- sca_per_ext$celltypes
# Split the SCE object by cell sca_per_ext
split_sca_per_ext <- split(seq_len(ncol(sca_per_ext)), cell_types)
# Create standalone SCE objects for each cell type
sca_per_ext_list <- lapply(split_sca_per_ext, function(indices) {
  sca_per_ext[, indices]
})
per_ext_Bcells<-sca_per_ext_list$'B cells'
per_ext_ASCs<-sca_per_ext_list$'ASC'
per_ext_GD_Tcells<-sca_per_ext_list$'CD2- GD T cells'
per_ext_CD2pos_GD<-sca_per_ext_list$'CD2+ GD T cells'
per_ext_CD4_Tcells<-sca_per_ext_list$'CD4+ ab T cells'
per_ext_Monocytes<-sca_per_ext_list$'Monocytes'
per_ext_pDCs<-sca_per_ext_list$'pDCs'
per_ext_CD8_NK<- sca_per_ext_list$'CD8a+ ab T/NK cells'
per_ext_NK<-sca_per_ext_list$'NK cells'
per_ext_CD8<-sca_per_ext_list$'CD8ab+ ab T cells'


### B cells in persistent vs extinct ##

cdr <-colSums(assay(per_ext_Bcells)>0)
colData(per_ext_Bcells)$cdr<- scale(cdr)
#In order to have more interpretable coefficients, we'll set the reference level of the factor to be the "unstimulated" cells.

Treatment <- factor(colData(per_ext_Bcells)$Treatment)
Treatment<-relevel(Treatment,"persistent")
colData(per_ext_Bcells)$Treatment<-Treatment
#Run MAST
zlmCond <- zlm(formula = ~ Treatment + cdr + (1 | Sow), 
               sca = per_ext_Bcells, 
               method = 'glmer', 
               ebayes = FALSE, 
               strictConvergence = FALSE,
               fitArgsD = list(nAGQ = 0))
#only test the condition coefficient via contrast matrix
summaryCond <- summary(zlmCond, doLRT='Treatmentextinct') 
##print the top 4 genes by contrast using the logFC
print(summaryCond, n=4)
## by discrete Z-score
print(summaryCond, n=4, by='D')
## by continuous Z-score
print(summaryCond, n=4, by='C')
##Make data table of results
summaryDt <- summaryCond$datatable
fcHurdle <- merge(
  summaryDt[contrast == 'Treatmentextinct' & component == 'H', .(primerid, `Pr(>Chisq)`)],
  summaryDt[contrast == 'Treatmentextinct' & component == 'logFC', .(primerid, logFC = coef, ci.hi, ci.lo)],
  by = 'primerid'
)
#add fdr
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
#add gene names via primerid
fcHurdleSig <- merge(fcHurdle, as.data.table(mcols(per_ext_Bcells)), by='primerid')
#order by fdr
setorder(fcHurdleSig, fdr)
fcHurdleSig_DE<-fcHurdleSig[fdr < 0.05]
## number of significant genes
cat("# of genes with FDR < 0.05 in B cells Extinct vs Persistent : ", nrow(fcHurdleSig_DE), "\n")

## plot with inverse logit transformed x-axis
ggplot(predicted_sig)+aes(x=invlogit(etaD),y=muC,xse=seD,yse=seC,col=sample)+
  facet_wrap(~primerid,scales="free_y")+theme_linedraw()+
  geom_point(size=0.5)+scale_x_continuous("Proportion expression")+
  scale_y_continuous("Estimated Mean")+
  stat_ell(aes(x=etaD,y=muC),level=0.95, invert='x')

mat_to_plot <- assay(ext_con_Bcells[entrez_to_plot,])
rownames(mat_to_plot) <- symbols_to_plot
heatmap(mat_to_plot,annCol=colData(ext_con_Bcells)[,"Treatment"],main="DE genes",col=rev(colorRampPalette(colors = brewer.pal(name="PiYG",n=10))(20)))

### ASC cells in persistent vs extinct ##

cdr <-colSums(assay(per_ext_ASCs)>0)
colData(per_ext_ASCs)$cdr<- scale(cdr)

Treatment <- factor(colData(per_ext_ASCs)$Treatment)
Treatment<-relevel(Treatment,"persistent")
colData(per_ext_ASCs)$Treatment<-Treatment
#Run per_ext_ASCs
zlmCond <- zlm(formula = ~ Treatment + cdr + (1 | Sow), 
               sca = per_ext_ASCs, 
               method = 'glmer', 
               ebayes = FALSE, 
               strictConvergence = FALSE,
               fitArgsD = list(nAGQ = 0))

#only test the condition coefficient via contrast matrix
summaryCond <- summary(zlmCond, doLRT='Treatmentextinct') 
##print the top 4 genes by contrast using the logFC
print(summaryCond, n=4)
## by discrete Z-score
print(summaryCond, n=4, by='D')
## by continuous Z-score
print(summaryCond, n=4, by='C')
##Make data table of results
summaryDt <- summaryCond$datatable
fcHurdle <- merge(
  summaryDt[contrast == 'Treatmentextinct' & component == 'H', .(primerid, `Pr(>Chisq)`)],
  summaryDt[contrast == 'Treatmentextinct' & component == 'logFC', .(primerid, logFC = coef, ci.hi, ci.lo)],
  by = 'primerid'
)
#add fdr
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
#add gene names via primerid
fcHurdleSig <- merge(fcHurdle, as.data.table(mcols(per_ext_ASCs)), by='primerid')
#order by fdr
setorder(fcHurdleSig, fdr)
fcHurdleSig_DE<-fcHurdleSig[fdr < 0.05]
## number of significant genes
cat("# of genes with FDR < 0.05 in ASC cells Extinct vs Persistent: ", nrow(fcHurdleSig_DE), "\n")


### CD2- GD T cells in persistent vs extinct ##

cdr <-colSums(assay(per_ext_GD_Tcells)>0)
colData(per_ext_GD_Tcells)$cdr<- scale(cdr)

Treatment <- factor(colData(per_ext_GD_Tcells)$Treatment)
Treatment<-relevel(Treatment,"persistent")
colData(per_ext_GD_Tcells)$Treatment<-Treatment
#Run per_ext_GD_Tcells
zlmCond <- zlm(formula = ~ Treatment + cdr + (1 | Sow), 
               sca =per_ext_GD_Tcells, 
               method = 'glmer', 
               ebayes = FALSE, 
               strictConvergence = FALSE,
               fitArgsD = list(nAGQ = 0))
#only test the condition coefficient via contrast matrix
summaryCond <- summary(zlmCond, doLRT='Treatmentextinct') 
##print the top 4 genes by contrast using the logFC
print(summaryCond, n=4)
## by discrete Z-score
print(summaryCond, n=4, by='D')
## by continuous Z-score
print(summaryCond, n=4, by='C')
##Make data table of results
summaryDt <- summaryCond$datatable
fcHurdle <- merge(
  summaryDt[contrast == 'Treatmentextinct' & component == 'H', .(primerid, `Pr(>Chisq)`)],
  summaryDt[contrast == 'Treatmentextinct' & component == 'logFC', .(primerid, logFC = coef, ci.hi, ci.lo)],
  by = 'primerid'
)
#add fdr
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
#add gene names via primerid
fcHurdleSig <- merge(fcHurdle, as.data.table(mcols(per_ext_GD_Tcells)), by='primerid')
#order by fdr
setorder(fcHurdleSig, fdr)
fcHurdleSig_DE<-fcHurdleSig[fdr < 0.05]
## number of significant genes
cat("# of genes with FDR < 0.05 in CD2-GD T cells Extinct vs Persistent: ", nrow(fcHurdleSig_DE), "\n")



### CD2+ GD T cells in persistent vs extinct ##

cdr <-colSums(assay(per_ext_CD2pos_GD)>0)
colData(per_ext_CD2pos_GD)$cdr<- scale(cdr)

Treatment <- factor(colData(per_ext_CD2pos_GD)$Treatment)
Treatment<-relevel(Treatment,"persistent")
colData(per_ext_CD2pos_GD)$Treatment<-Treatment
#Run per_ext_CD2pos_GD
zlmCond <- zlm(formula = ~ Treatment + cdr + (1 | Sow), 
               sca = per_ext_CD2pos_GD, 
               method = 'glmer', 
               ebayes = FALSE, 
               strictConvergence = FALSE,
               fitArgsD = list(nAGQ = 0))
#summary(zlm(~ Treatment + (1|Sow), sca = per_ext_CD2pos_GD, method = "glmer", ebayes = FALSE, fitArgsD = list(nAGQ = 0)))$datatable
#only test the condition coefficient via contrast matrix
summaryCond <- summary(zlmCond, doLRT='Treatmentextinct') 
##print the top 4 genes by contrast using the logFC
print(summaryCond, n=4)
## by discrete Z-score
print(summaryCond, n=4, by='D')
## by continuous Z-score
print(summaryCond, n=4, by='C')
##Make data table of results
summaryDt <- summaryCond$datatable
fcHurdle <- merge(
  summaryDt[contrast == 'Treatmentextinct' & component == 'H', .(primerid, `Pr(>Chisq)`)],
  summaryDt[contrast == 'Treatmentextinct' & component == 'logFC', .(primerid, logFC = coef, ci.hi, ci.lo)],
  by = 'primerid'
)
#add fdr
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
#add gene names via primerid
fcHurdleSig <- merge(fcHurdle, as.data.table(mcols(per_ext_CD2pos_GD)), by='primerid')
#order by fdr
setorder(fcHurdleSig, fdr)
fcHurdleSig_DE<-fcHurdleSig[fdr < 0.05]
## number of significant genes
cat("# of genes with FDR < 0.05 in CD2+GD T cells Extinct vs Persistent: ", nrow(fcHurdleSig_DE), "\n")


### CD4 AB T cells in persistent vs extinct ##

cdr <-colSums(assay(per_ext_CD4_Tcells)>0)
colData(per_ext_CD4_Tcells)$cdr<- scale(cdr)
Treatment <- factor(colData(per_ext_CD4_Tcells)$Treatment)
Treatment<-relevel(Treatment,"persistent")
colData(per_ext_CD4_Tcells)$Treatment<-Treatment
#Run per_ext_CD4_Tcells
zlmCond <- zlm(formula = ~ Treatment  +cdr+ (1 | Sow), 
               sca = per_ext_CD4_Tcells, 
               method = 'glmer', 
               ebayes = FALSE, 
               strictConvergence = FALSE,
               fitArgsD = list(nAGQ = 0))
#summary(zlm(~ Treatment + (1|Sow), sca = per_ext_CD4_Tcells, method = "glmer", ebayes = FALSE, fitArgsD = list(nAGQ = 0)))$datatable
#only test the condition coefficient via contrast matrix
summaryCond <- summary(zlmCond, doLRT='Treatmentextinct') 
a##print the top 4 genes by contrast using the logFC
print(summaryCond, n=4)
## by discrete Z-score
print(summaryCond, n=4, by='D')
## by continuous Z-score
print(summaryCond, n=4, by='C')
##Make data table of results
summaryDt <- summaryCond$datatable
fcHurdle <- merge(
  summaryDt[contrast == 'Treatmentextinct' & component == 'H', .(primerid, `Pr(>Chisq)`)],
  summaryDt[contrast == 'Treatmentextinct' & component == 'logFC', .(primerid, logFC = coef, ci.hi, ci.lo)],
  by = 'primerid'
)
#add fdr
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
#add gene names via primerid
fcHurdleSig <- merge(fcHurdle, as.data.table(mcols(per_ext_CD4_Tcells)), by='primerid')
#order by fdr
setorder(fcHurdleSig, fdr)
## number of significant genes
fcHurdleSig_DE<-fcHurdleSig[fdr < 0.05]
## number of significant genes
cat("# of genes with FDR < 0.05 in CD4 T cells Extinct vs Persistent: ", nrow(fcHurdleSig_DE), "\n")

### Monocytes in persistent vs extinct ##

cdr <-colSums(assay(per_ext_Monocytes)>0)
colData(per_ext_Monocytes)$cdr<- scale(cdr)

Treatment <- factor(colData(per_ext_Monocytes)$Treatment)
Treatment<-relevel(Treatment,"persistent")
colData(per_ext_Monocytes)$Treatment<-Treatment
#Run per_con_Monocytes
zlmCond <- zlm(formula = ~ Treatment  +cdr+ (1 | Sow), 
               sca = per_ext_Monocytes, 
               method = 'glmer', 
               ebayes = FALSE, 
               strictConvergence = FALSE,
               fitArgsD = list(nAGQ = 0))
#summary(zlm(~ Treatment + (1|Sow), sca = per_ext_Monocytes, method = "glmer", ebayes = FALSE, fitArgsD = list(nAGQ = 0)))$datatable
#only test the condition coefficient via contrast matrix
summaryCond <- summary(zlmCond, doLRT='Treatmentextinct') 
##print the top 4 genes by contrast using the logFC
print(summaryCond, n=4)
## by discrete Z-score
print(summaryCond, n=4, by='D')
## by continuous Z-score
print(summaryCond, n=4, by='C')
##Make data table of results
summaryDt <- summaryCond$datatable
fcHurdle <- merge(
  summaryDt[contrast == 'Treatmentextinct' & component == 'H', .(primerid, `Pr(>Chisq)`)],
  summaryDt[contrast == 'Treatmentextinct' & component == 'logFC', .(primerid, logFC = coef, ci.hi, ci.lo)],
  by = 'primerid'
)
#add fdr
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
#add gene names via primerid
fcHurdleSig <- merge(fcHurdle, as.data.table(mcols(per_ext_Monocytes)), by='primerid')
#order by fdr
setorder(fcHurdleSig, fdr)
fcHurdleSig_DE<-fcHurdleSig[fdr < 0.05]
## number of significant genes
cat("# of genes with FDR < 0.05 in Mono cells Extinct vs Persistent: ", nrow(fcHurdleSig_DE), "\n")

### pDCs in Persistent vs extinct ##

cdr <-colSums(assay(per_ext_pDCs)>0)
colData(per_ext_pDCs)$cdr<- scale(cdr)

Treatment <- factor(colData(per_ext_pDCs)$Treatment)
Treatment<-relevel(Treatment,"persistent")
colData(per_ext_pDCs)$Treatment<-Treatment
#Run per_ext_pDCs
zlmCond <- zlm(formula = ~ Treatment  +cdr+ (1 | Sow), 
               sca = per_ext_pDCs, 
               method = 'glmer', 
               ebayes = FALSE, 
               strictConvergence = FALSE,
               fitArgsD = list(nAGQ = 0))
#summary(zlm(~ Treatment + (1|Sow), sca = ext_con_pDCs, method = "glmer", ebayes = FALSE, fitArgsD = list(nAGQ = 0)))$datatable
#only test the condition coefficient via contrast matrix
summaryCond <- summary(zlmCond, doLRT='Treatmentextinct') 
##print the top 4 genes by contrast using the logFC
print(summaryCond, n=4)
## by discrete Z-score
print(summaryCond, n=4, by='D')
## by continuous Z-score
print(summaryCond, n=4, by='C')
##Make data table of results
summaryDt <- summaryCond$datatable
fcHurdle <- merge(
  summaryDt[contrast == 'Treatmentextinct' & component == 'H', .(primerid, `Pr(>Chisq)`)],
  summaryDt[contrast == 'Treatmentextinct' & component == 'logFC', .(primerid, logFC = coef, ci.hi, ci.lo)],
  by = 'primerid'
)
#add fdr
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
#add gene names via primerid
fcHurdleSig <- merge(fcHurdle, as.data.table(mcols(per_ext_pDCs)), by='primerid')
#order by fdr
setorder(fcHurdleSig, fdr)
fcHurdleSig_DE<-fcHurdleSig[fdr < 0.05]
## number of significant genes
cat("# of genes with FDR < 0.05 in pDCs cells Extinct vs Persistent: ", nrow(fcHurdleSig_DE), "\n")

### CD8/NK in Persistent vs extinct ##

cdr <-colSums(assay(per_ext_CD8_NK)>0)
colData(per_ext_CD8_NK)$cdr<- scale(cdr)

Treatment <- factor(colData(per_ext_CD8_NK)$Treatment)
Treatment<-relevel(Treatment,"persistent")
colData(per_ext_CD8_NK)$Treatment<-Treatment
#Run per_ext_CD8_NK
zlmCond <- zlm(formula = ~ Treatment  +cdr+ (1 | Sow), 
               sca = per_ext_CD8_NK, 
               method = 'glmer', 
               ebayes = FALSE, 
               strictConvergence = FALSE,
               fitArgsD = list(nAGQ = 0))
#summary(zlm(~ Treatment + (1|Sow), sca = per_ext_CD8_NK, method = "glmer", ebayes = FALSE, fitArgsD = list(nAGQ = 0)))$datatable
#only test the condition coefficient via contrast matrix
summaryCond <- summary(zlmCond, doLRT='Treatmentextinct') 
##print the top 4 genes by contrast using the logFC
print(summaryCond, n=4)
## by discrete Z-score
print(summaryCond, n=4, by='D')
## by continuous Z-score
print(summaryCond, n=4, by='C')
##Make data table of results
summaryDt <- summaryCond$datatable
fcHurdle <- merge(
  summaryDt[contrast == 'Treatmentextinct' & component == 'H', .(primerid, `Pr(>Chisq)`)],
  summaryDt[contrast == 'Treatmentextinct' & component == 'logFC', .(primerid, logFC = coef, ci.hi, ci.lo)],
  by = 'primerid'
)
#add fdr
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
#add gene names via primerid
fcHurdleSig <- merge(fcHurdle, as.data.table(mcols(per_ext_CD8_NK)), by='primerid')
#order by fdr
setorder(fcHurdleSig, fdr)
fcHurdleSig_DE<-fcHurdleSig[fdr < 0.05]
## number of significant genes
cat("# of genes with FDR < 0.05 in CD8 & NK cells Extinct vs Persistent: ", nrow(fcHurdleSig_DE), "\n")


### NK in Persistent vs extinct ##

cdr <-colSums(assay(per_ext_NK)>0)
colData(per_ext_NK)$cdr<- scale(cdr)

Treatment <- factor(colData(per_ext_NK)$Treatment)
Treatment<-relevel(Treatment,"persistent")
colData(per_ext_NK)$Treatment<-Treatment
#Run per_ext_NK
zlmCond <- zlm(formula = ~ Treatment  +cdr+ (1 | Sow), 
               sca = per_ext_NK, 
               method = 'glmer', 
               ebayes = FALSE, 
               strictConvergence = FALSE,
               fitArgsD = list(nAGQ = 0))
#summary(zlm(~ Treatment + (1|Sow), sca = per_ext_NK, method = "glmer", ebayes = FALSE, fitArgsD = list(nAGQ = 0)))$datatable
#only test the condition coefficient via contrast matrix
summaryCond <- summary(zlmCond, doLRT='Treatmentextinct') 
##print the top 4 genes by contrast using the logFC
print(summaryCond, n=4)
## by discrete Z-score
print(summaryCond, n=4, by='D')
## by continuous Z-score
print(summaryCond, n=4, by='C')
##Make data table of results
summaryDt <- summaryCond$datatable
fcHurdle <- merge(
  summaryDt[contrast == 'Treatmentextinct' & component == 'H', .(primerid, `Pr(>Chisq)`)],
  summaryDt[contrast == 'Treatmentextinct' & component == 'logFC', .(primerid, logFC = coef, ci.hi, ci.lo)],
  by = 'primerid'
)
#add fdr
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
#add gene names via primerid
fcHurdleSig <- merge(fcHurdle, as.data.table(mcols(per_ext_NK)), by='primerid')
#order by fdr
setorder(fcHurdleSig, fdr)
fcHurdleSig_DE<-fcHurdleSig[fdr < 0.05]
## number of significant genes
cat("# of genes with FDR < 0.05 in NK cells Extinct vs Persistent: ", nrow(fcHurdleSig_DE), "\n")


### CD8T in Persistent vs extinct ##

cdr <-colSums(assay(per_ext_CD8)>0)
colData(per_ext_CD8)$cdr<- scale(cdr)

Treatment <- factor(colData(per_ext_CD8)$Treatment)
Treatment<-relevel(Treatment,"persistent")
colData(per_ext_CD8)$Treatment<-Treatment
#Run per_ext_CD8
zlmCond <- zlm(formula = ~ Treatment  +cdr+ (1 | Sow), 
               sca = per_ext_CD8, 
               method = 'glmer', 
               ebayes = FALSE, 
               strictConvergence = FALSE,
               fitArgsD = list(nAGQ = 0))
#summary(zlm(~ Treatment + (1|Sow), sca = per_ext_CD8, method = "glmer", ebayes = FALSE, fitArgsD = list(nAGQ = 0)))$datatable
#only test the condition coefficient via contrast matrix
summaryCond <- summary(zlmCond, doLRT='Treatmentextinct') 
##print the top 4 genes by contrast using the logFC
print(summaryCond, n=4)
## by discrete Z-score
print(summaryCond, n=4, by='D')
## by continuous Z-score
print(summaryCond, n=4, by='C')
##Make data table of results
summaryDt <- summaryCond$datatable
fcHurdle <- merge(
  summaryDt[contrast == 'Treatmentextinct' & component == 'H', .(primerid, `Pr(>Chisq)`)],
  summaryDt[contrast == 'Treatmentextinct' & component == 'logFC', .(primerid, logFC = coef, ci.hi, ci.lo)],
  by = 'primerid'
)
#add fdr
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
#add gene names via primerid
fcHurdleSig <- merge(fcHurdle, as.data.table(mcols(per_ext_CD8)), by='primerid')
#order by fdr
setorder(fcHurdleSig, fdr)
fcHurdleSig_DE<-fcHurdleSig[fdr < 0.05]
## number of significant genes
cat("# of genes with FDR < 0.05 in CD8 cells Extinct vs Persistent: ", nrow(fcHurdleSig_DE), "\n")
