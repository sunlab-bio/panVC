library(Seurat)
library(ggplot2)
library(tidyverse)
library(qs)
library(EnsDb.Hsapiens.v86)
source("cat_metacells.R")

#### 1.metacell ----
k <- 30
seurat_obj <- qs::qread("./results/files/filter/EC_trans_sub_filter.qs")
DefaultAssay(seurat_obj) <- 'RNA'
seurat_obj_metacell <- seurat_obj %>%
  cat_construct_metacells_sub_cell_type(k = k,
                                        name = i,
                                        min_cells = 100)

metacell_obj <- hdWGCNA::GetMetacellObject(seurat_obj_metacell)
metacell_obj$ID <- rownames(metacell_obj@meta.data)
metacell_obj$cell_type <- 'Endothelial cells'

####
inter <- intersect(names(table(seurat_obj$sample)),names(table(metacell_obj$sample)))
rownames(metacell_obj@meta.data) <- metacell_obj$ID

####
counts <- readRDS( file.path(outdir, paste0(k, '_counts.rds')))
metadata <- readRDS(file.path(outdir, paste0(k,'_metadata.rds')))

df <- as.data.frame(mapIds(EnsDb.Hsapiens.v86,
                           keys = rownames(counts), # the keys to select records for from the database.
                           keytype="SYMBOL",
                           column="GENEID" # the column to search on (for mapIds).
))
df$SYMBOL <- rownames(df)
colnames(df) <- c('GENEID', 'SYMBOL')
df <- na.omit(df) 
counts2 <- counts[df$SYMBOL, ]
rownames(counts2) <- df$GENEID

##
seurat_metacell_obj <- CreateSeuratObject(counts2)
seurat_metacell_obj@meta.data <- metadata
qs::qsave(seurat_metacell_obj, 'metacell_seuratobj_ID.qs')

#### 2.merge ----
files_list <- c(seurat_obj1, seurat_obj2, seurat_obj3, seurat_obj4, seurat_obj5, seurat_obj6,
                seurat_obj7, seurat_obj8, seurat_obj9, seurat_obj10, seurat_obj11)
for(i in 1:11){
  assign(paste0("gene", i), rownames(files_list[i]))}
gene_lists <- lapply(files_list, rownames)

common_genes <- Reduce(intersect, gene_lists)
for (i in 1:11) {
  files_list[[i]] <- files_list[[i]][common_genes, ]
}
seurat_obj0 <- merge(files_list[[1]],
                     y = files_list[-1])

####
df <- as.data.frame(mapIds(EnsDb.Hsapiens.v86,
                           keys = rownames(seurat_obj0), # the keys to select records for from the database.
                           keytype="GENEID",
                           column="SYMBOL" # the column to search on (for mapIds).
))
df$GENEID <- rownames(df)
colnames(df) <- c('SYMBOL', 'GENEID')
df <- na.omit(df) 

##
rownames(seurat_obj0@assays$RNA@counts) <- df$SYMBOL
rownames(seurat_obj0@assays$RNA@data) <- df$SYMBOL
rownames(seurat_obj0@assays$RNA@meta.features) <- df$SYMBOL
qs::qsave(seurat_obj0, file.path(outdir,paste0(case, "Endothelial_ori_merged.qs")))

####
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
seurat_obj <-
  CellCycleScoring(
    seurat_obj,
    s.features = s.genes,
    g2m.features = g2m.genes,
    set.ident = FALSE)
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", 'S.Score', 'G2M.Score'),
        ncol = 2, pt.size = 0)
qs::qsave(seurat_obj, file.path(outdir, paste0(case, "Endothelial_ok_merged.qs")))

#### 3.SCT ----
seurat_obj_list <- SplitObject(seurat_obj, split.by = 'dataset')
seurat_obj_list <- lapply(seurat_obj_list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = 2200)
  return(x)
})
seurat_obj_list <- lapply(X = seurat_obj_list, FUN = SCTransform)
seurat_obj_list <- lapply(
  X = seurat_obj_list,
  FUN = function(x) {
    x <- SCTransform(
      x,
      method = "glmGamPoi",
      vars.to.regress = c("nCount_RNA", "sample",'dis_score'),
      vst.flavor = "v2",
      verbose = T
    )})

features <- SelectIntegrationFeatures(seurat_obj_list,
                                      nfeatures = 2000)
dims.use=1:15
exclude_genes <- read.table('~/Collection/Data_P_H/files/exclude_genes.txt') |> pull()
HVGs <- features[!(features %in% exclude_genes)][1:1800]
seurat_obj_list <-  PrepSCTIntegration(object.list = seurat_obj_list, anchor.features = HVGs)
immune.anchors <- FindIntegrationAnchors(object.list = seurat_obj_list, normalization.method = "SCT",
                                         anchor.features = HVGs, dims = dims.use)
seurat_obj <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

DefaultAssay(seurat_obj) <- "integrated"
seurat_obj <- RunPCA(seurat_obj)

set.seed(016)
seurat_obj <- seurat_obj |>
  # RunTSNE(dims = dim.use, reduction = "pca") |>
  RunUMAP(dims = dim.use, reduction = "pca") |>
  FindNeighbors(dims = dim.use, reduction = "pca") |>
  FindClusters(resolution = 2)

DefaultAssay(seurat_obj) <- 'RNA'
VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", 
               "rb_score", "hs_score", "hb_score", 
               "dis_score", "rna_MALAT1"),
  pt.size = 0,ncol = 1,
  group.by = "seurat_clusters"
)
p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "source",
              cols = my36colors) + theme_bw()
p2 <-  DimPlot(seurat_obj, reduction = "umap", group.by = "dataset",
               cols = my36colors) + theme_bw() 
p3 <- DimPlot(seurat_obj, reduction = "umap", group.by = "sample") + theme_bw()
p4 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters",) + theme_bw()
p5 <- DimPlot(seurat_obj, cols = my36colors,  
              label = T, label.size = 3,label.box = T, repel = T,
              reduction = "umap", group.by = "seurat_clusters",) + theme_bw()
p6 <- FeaturePlot(seurat_obj,features= 'hs_score', max.cutoff = 0.75)+
  scale_colour_distiller(palette = "YlOrRd",direction = 1)
p7 <- FeaturePlot(seurat_obj,features= 'rna_MALAT1',
                  slot = 'data',min.cutoff = 4,max.cutoff = 8)+
  scale_colour_distiller(palette = "YlOrRd",direction = 1)
p8 <- FeaturePlot(seurat_obj,features= 'dis_score',max.cutoff = 5)+
  scale_colour_distiller(palette = "YlOrRd",direction = 1)
#
p <-p1 + p2 + p3 + p4 +plot_layout(nrow = 2, ncol = 2,guides = "keep")
ggsave( file.path(outdir2, paste0(case, "umap_data.pdf")),plot = p,height = 8,width = 12)
p <-p5 + p6 + p7 + p8 +plot_layout(nrow = 2, ncol = 2,guides = "keep")
ggsave(file.path(outdir2, paste0(case, "umap_score.pdf")),  plot = p,  height = 6,  width = 8)

p <- DimPlot(seurat_obj, 
             label = T, label.size = 3,label.box = T, repel = T,
             reduction = "umap", group.by = "seurat_clusters",) + theme_bw() +
  labs( x= "UMAP 1",y= "UMAP 2",title = "seurat_clusters") +
  theme(panel.grid = element_blank(),aspect.ratio = 1)
ggsave(file.path(outdir2,paste0(case, 'seurat_clusters.pdf')),p,height=7,width=7)

####
lisi_score <- lisi::compute_lisi(
  X = Embeddings(seurat_obj, reduction = "pca"),
  meta_data = seurat_obj@meta.data,
  label_colnames =  "dataset"
) |>
  dplyr::mutate(type = "pca") |>
  rownames_to_column("ID")
DimPlot(seurat_obj, reduction = "umap", group.by = "dataset",) + theme_bw() +
  labs( x= "UMAP 1",y= "UMAP 2",title = "dataset")

#
markers <- c('SEMA3G', 'PCSK5', 'NEBL',
             'RGCC', 'CA4', 'SLC9C1', 
             'NR2F2', 'TSHZ2', 'SLCO2A1', 
             'LYVE1', 'TBX1', 'PROX1')
for(i in 1:length(markers)){
  assign(paste0("p", i), FeaturePlot(seurat_obj,features= markers[i]))
}
p <-p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9+ p10+ p11 + p12 +
  plot_layout(ncol = 3)# ,guides = "collect"
ggsave(file.path(outdir2, paste0(case,"FeaturePlot.pdf")),
       plot = p,  height = 6,  width = 6)

#### 4.name ----
seurat_obj <- qread('./files/sample_15_1800_SCT_pca.qs')
DimPlot(seurat_obj, label = T,label.box = T,
        reduction = "umap", group.by = "seurat_clusters",)

dim.use <- 1:15
Resolution <- c(0.3,0.5,0.7,1,1.2,1.5)
seurat_obj <- seurat_obj |>
  # RunTSNE(dims = dim.use, reduction = "pca") |>
  # RunUMAP(dims = dim.use, reduction = "pca") |>
  FindNeighbors(dims = dim.use, reduction = "pca") |>
  FindClusters(resolution = Resolution)

for(i in as.character(Resolution)){
  p <- DimPlot(seurat_obj, cols = my36colors,label = T,label.box = T,
               reduction = "umap", group.by = paste0("integrated_snn_res.",i),) + theme_bw() +
    labs( x= "UMAP 1",y= "UMAP 2",title = paste0("integrated_snn_res.",i)) +
    theme(panel.grid = element_blank(),aspect.ratio = 1)+NoLegend() 
  ggsave(file.path(outdir2, paste0(case, paste0("UMAP_", i, ".pdf"))),
         p, height=5, width=5)
}

seurat_obj$seurat_clusters <- seurat_obj$integrated_snn_res.0.7
DimPlot(seurat_obj,label = T, label.box = T,
        group.by = 'integrated_snn_res.0.7',
        reduction = "umap") + theme_bw()

feature <- c('VWF', 'CDH5', 'FLT1',
             'SEMA3G', 'PCSK5', 'NEBL',
             'RGCC', 'CA4', 'SLC9C1', 
             'NR2F2', 'TSHZ2', 'SLCO2A1', 
             'CX3CL1', 'CCL2', 'IL6', 'CXCL3', 
             'POSTN', 'TMEM132C', 'CDH11',
             'FLT4','LYVE1', 'TBX1')
seurat_obj$classic1 <- "Capillary endo"
seurat_obj$classic1[seurat_obj$seurat_clusters%in% c('3', '11','12')] <- "Arterial endo"
seurat_obj$classic1[seurat_obj$seurat_clusters%in% c('9', '13')] <- 'Venous endo'
seurat_obj$classic1[seurat_obj$seurat_clusters%in% c('14')] <- 'Endocardial endo'
seurat_obj$classic1[seurat_obj$seurat_clusters%in% c('15')] <- 'Lymphatic endo'
DimPlot(seurat_obj, cols = classic1_colors, group.by = 'classic1',
        label = T, label.box = T, reduction = "umap")
#
feature <- c('DKK2','NEBL','VEGFA','PAPSS2','RGCC',
             'RDH10','CA2','TNFRSF4','RGS6','PLVAP', 
             'ACTA2','FLT4','CDH11')
DotPlot(seurat_obj, features = feature)  + coord_flip()+
  scale_colour_distiller(palette = "YlOrRd", direction = 1)
seurat_obj$type3 <- ''
new.cluster.ids <- c('DKK2','NEBL','VEGFA','PAPSS2','RGCC',
                     'RDH10','CA2','TNFRSF4','RGS6','PLVAP', 
                     'ACTA2','FLT4','CDH11')
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj,
                           new.cluster.ids)
DimPlot(seurat_obj,group.by = 'type3', label = T,reduction = 'umap',
        label.box = T,repel = T,label.size = 3 )
ggsave(file.path(outdir2, paste0(case,'Dimplot_type3.pdf')),
       ggplot2::last_plot(),height=5,width=5)
qs::qsave(seurat_obj, file.path(outdir,paste0(case, "named_type.qs")))
