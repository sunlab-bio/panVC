library(Seurat)
library(ggplot2)
library(tidyverse)
library(qs)
library(EnsDb.Hsapiens.v86)
source("~/Collection/code/scrna-seq.R")

#### 1.metacell ----
outdir <- '~/PH/results/named/files/'
outdir2 <- '~/PH/results/named/plots/'

k <- 30
seurat_obj <- qs::qread("./filter/EC_trans_sub_filter.qs")
DefaultAssay(seurat_obj) <- 'RNA'
seurat_obj_metacell <- seurat_obj %>%
  cat_construct_metacells_sub_cell_type(k = k,
                                        name = i,
                                        min_cells = 100)

metacell_obj <- hdWGCNA::GetMetacellObject(seurat_obj_metacell)
inter <- intersect(names(table(seurat_obj$sample)),names(table(metacell_obj$sample)))
rownames(metacell_obj@meta.data) <- metacell_obj$ID
saveRDS(metacell_obj,file.path(outdir, paste0(k, '_metacell_obj.rds')))

## IDtrans
counts <- metacell_obj@assays$RNA@counts
metadata <- metacell_obj@meta.data
df <- as.data.frame(mapIds(EnsDb.Hsapiens.v86,
                           keys = rownames(counts),
                           keytype="SYMBOL",
                           column="GENEID"
))
df$SYMBOL <- rownames(df)
colnames(df) <- c('GENEID', 'SYMBOL')
df <- na.omit(df) 
counts2 <- counts[df$SYMBOL, ]
rownames(counts2) <- df$GENEID

# CreateSeuratObject
seurat_metacell_obj <- CreateSeuratObject(counts2)
seurat_metacell_obj@meta.data <- metadata
qs::qsave(seurat_metacell_obj,file.path(outdir, paste0(k, 'metacell_seuratobj_ID.qs')))

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

df <- as.data.frame(mapIds(EnsDb.Hsapiens.v86,
                           keys = rownames(seurat_obj0),
                           keytype="GENEID",
                           column="SYMBOL"))
df$GENEID <- rownames(df)
colnames(df) <- c('SYMBOL', 'GENEID')
df <- na.omit(df) 

rownames(seurat_obj0@assays$RNA@counts) <- df$SYMBOL
rownames(seurat_obj0@assays$RNA@data) <- df$SYMBOL
rownames(seurat_obj0@assays$RNA@meta.features) <- df$SYMBOL
qs::qsave(seurat_obj0, file.path(outdir,"Endothelial_ori_merged.qs"))

#### 3.SCT ----
seurat_obj <- seurat_obj0
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
exclude_genes <- read.table('exclude_genes.txt') |> pull()
HVGs <- features[!(features %in% exclude_genes)][1:1800]
seurat_obj_list <-  PrepSCTIntegration(object.list = seurat_obj_list, anchor.features = HVGs)
immune.anchors <- FindIntegrationAnchors(object.list = seurat_obj_list, normalization.method = "SCT",
                                         anchor.features = HVGs, dims = dims.use)
seurat_obj <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
DefaultAssay(seurat_obj) <- "integrated"
seurat_obj <- RunPCA(seurat_obj)

set.seed(016)
seurat_obj <- seurat_obj |>
  RunUMAP(dims = dim.use, reduction = "pca") |>
  FindNeighbors(dims = dim.use, reduction = "pca") |>
  FindClusters(resolution = 2)
qs::qsave(seurat_obj, './files/sample_15_1800_SCT_pca.qs')

DefaultAssay(seurat_obj) <- 'RNA'
p <- DimPlot(seurat_obj, 
             label = T, label.size = 3,label.box = T, repel = T,
             reduction = "umap", group.by = "seurat_clusters",)
ggsave(file.path(outdir2, 'seurat_clusters.pdf'),p,height=7,width=7)

## lisi_score
lisi_score <- lisi::compute_lisi(
  X = Embeddings(seurat_obj, reduction = "pca"),
  meta_data = seurat_obj@meta.data,
  label_colnames =  "dataset") |>
  dplyr::mutate(type = "pca") |>
  rownames_to_column("ID")
round(mean(lisi_score$dataset), 2)
write_csv(lisi_score, file = "lisi_score.csv")


#### 4.named ----
## remove non-endothelial cells
seurat_obj <- qread('./files/sample_15_1800_SCT_pca.qs')
DimPlot(seurat_obj, label = T,label.box = T,
        reduction = "umap", group.by = "seurat_clusters",)
seurat_obj2 <- subset(seurat_obj, idents = c('18', '23') ,invert = TRUE)

## Re-cluster
dim.use <- 1:15
Resolution <- c(0.3,0.5,0.7,1)
seurat_obj <- seurat_obj2 |>
  FindNeighbors(dims = dim.use, reduction = "pca") |>
  FindClusters(resolution = Resolution)

for(i in as.character(Resolution)){
  p <- DimPlot(seurat_obj, cols = my36colors,label = T,label.box = T,
               reduction = "umap", group.by = paste0("integrated_snn_res.",i),) + theme_bw() +
    labs( x= "UMAP 1",y= "UMAP 2",title = paste0("integrated_snn_res.",i)) +
    theme(panel.grid = element_blank(),aspect.ratio = 1)+NoLegend() 
  ggsave(file.path(outdir2, paste0("sample_15_1800_", "UMAP_", i, ".pdf")),
         p, height=5, width=5)
}
DefaultAssay(seurat_obj) <- 'RNA'
seurat_obj$seurat_clusters <- seurat_obj$integrated_snn_res.0.7
DimPlot(seurat_obj,label = T, label.box = T,
        group.by = 'integrated_snn_res.0.7',
        reduction = "umap") + theme_bw()


## major celltype
seurat_obj$classic1 <- "CapECs"
seurat_obj$classic1[seurat_obj$seurat_clusters%in% c('3', '11','12')] <- "ArtECs"
seurat_obj$classic1[seurat_obj$seurat_clusters%in% c('9', '13')] <- 'VenECs'
seurat_obj$classic1[seurat_obj$seurat_clusters%in% c('14')] <- 'EndoECs'
seurat_obj$classic1[seurat_obj$seurat_clusters%in% c('15')] <- 'LECs'
DimPlot(seurat_obj, group.by = 'classic1',
        label = T, label.box = T, reduction = "umap")

## 13 subsets
seurat_obj$type3 <- ''
new.cluster.ids <- c('DKK2','NEBL','VEGFA','PAPSS2','RGCC',
                     'RDH10','CA2','TMEM163','RGS6','VCAM1', 
                     'ACTA2','FLT4','NPR3')
names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj,
                           new.cluster.ids)
DimPlot(seurat_obj,group.by = 'type3', label = T,reduction = 'umap',
        label.box = T,repel = T,label.size = 3 )
ggsave(file.path(outdir2, 'sample_15_1800_Dimplot_type3.pdf'),
       ggplot2::last_plot(),height=5,width=5)
qs::qsave(seurat_obj, file.path(outdir, "sample_15_1800_named_type.qs"))
