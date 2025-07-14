library(Seurat)
library(tidyverse)
library(dplyr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(DESeq2)
library(ggplot2)
library(cowplot)
library(ggrepel)
set.seed(717)
source("~/Collection/code/scrna-seq.R")
source("~/Collection/code/plot.R")

path <- '~/PH/results/DEG/'
setwd(path)
outdir <- './files/'

#### 0.expr ----
seurat_obj0 <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
sce <- subset(seurat_obj0, idents = c("LEC.c12.FLT4","Endo.c13.NPR3"), invert = T)
DefaultAssay(sce) <- "RNA"
dim(sce)

sce$Group <- sce$disease
bs = split(colnames(sce),sce$sample)
ct = do.call(
  cbind,lapply(names(bs), function(x){ 
    # x=names(bs)[[1]]
    kp =colnames(sce) %in% bs[[x]]
    rowSums( as.matrix(sce@assays$RNA@counts[, kp]  ))
  })
)
colnames(ct) <- names(bs)
saveRDS(ct, file.path(outdir, 'bulk_sample_expr.rds'))

phe = unique(sce@meta.data[,c('sample','Group')])
phe
saveRDS(phe, file.path(outdir, 'bulk_metadata.rds'))

#### 1.deg ----
library(edgeR)

exp0 <- readRDS('./files/bulk_sample_expr.rds')
metadata <- readRDS('./files/bulk_metadata.rds')

List <- names(table(metadata$Group))[2:9] 
for(i in List){
  Sample_name <- metadata |> filter(Group %in% c('Normal', i)) |> pull(sample)
  length(intersect(colnames(exp0), metadata$sample))
  exp <- exp0[, Sample_name, drop = FALSE]
  dim(exp)
  head(exp)
  Group <- metadata |> filter(Group %in% c('Normal', i)) |> pull(Group)
  Group <- factor(Group, levels = c('Normal', i))
  
  d <- DGEList(counts = exp, group = Group)
  keep <- rowSums(cpm(d) > 1) >= 2
  d <- d[keep, , keep.lib.sizes = FALSE]
  d <- calcNormFactors(d)
  design <- model.matrix(~0 + Group)
  colnames(design) <- levels(Group)
  d <- estimateGLMCommonDisp(d, design)
  d <- estimateGLMTrendedDisp(d, design)
  d <- estimateGLMTagwiseDisp(d, design)
  fit <- glmFit(d, design)
  lrt <- glmLRT(fit, contrast = c(-1, 1))
  DEG_edgeR <- na.omit(as.data.frame(topTags(lrt, n = Inf)))
  DEG_edgeR$Gene_Symbol <- rownames(DEG_edgeR)
  
  markers <- DEG_edgeR %>%
    mutate(change = if_else(FDR > 0.05, 
                            "Stable",
                            if_else(abs(logFC) < 0.58, "Stable",
                                    if_else(logFC >= 0.58, "Upregulated", 
                                            "Downregulated"))))

  saveRDS(markers, file.path(outdir, paste0(i, "_deg.rds")))
  saveRDS(markers[markers$change == "Upregulated",], file.path(outdir, paste0(i, "_deg_up.rds")))
  saveRDS(markers[markers$change == "Downregulated",], file.path(outdir, paste0(i, "_deg_down.rds")))
}

all_markers_list <- list()
for (i in List) {
  tryCatch({
    markers <- readRDS(file.path(outdir, paste0(i, "_deg.rds")))
    markers$disease <- i
    all_markers_list[[length(all_markers_list) + 1]] <- markers
    
  }, error = function(e) {
    cat("Error processing ",i, ": ", e$message, "\n")
  })
}

all_markers_df <- do.call(rbind, all_markers_list)
saveRDS(all_markers_df, file.path(outdir, 'DEGs_disease_number.rds'))