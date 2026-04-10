library(Seurat)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(qs)
library(pheatmap)
library(tidydr)
library(dplyr)
library(ggbreak)
library(msigdbr)
library(AUCell)
library(clusterProfiler)
source('~/Collection/code/scrna-seq.R')
source('~/Collection/code/plot.R')

path <- '~/PH/results/Fig4/'
setwd(path)
outdir <- './files/'
outdir2 <- './plots/'

#### 1.Fig4a ----
seurat_obj <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
seurat_obj <- subset(seurat_obj, classic1 %in% c('ArtECs','CapECs','VenECs'))
Idents(seurat_obj) <- 'classic1'
DimPlot(seurat_obj)+
  theme_void() +
  theme_dr(xlength = 0.2, ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"), ends = 'last', type = "closed"))+
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs( x= "UMAP 1",y= "UMAP 2", title = '')+
  theme(plot.margin = margin(5.5,15,5.5,5.5),
        panel.grid = element_blank(),aspect.ratio = 1) +
  scale_color_manual(values = c("#edd064","#f2ccac","#a1d5b9"))
ggsave(file.path(outdir2, 'VECsubtype_UMAP.pdf'),
       ggplot2::last_plot(), height=4, width=6)

#### 2.Fig4b ----
library(plotly)
#### 2.1 CytoTRACE ----
#### DiffusionMap ----
root_type <- "VenECs"
col_names <- c("#edd064","#f2ccac","#a1d5b9")
DCmap_umap <- function(seu_obj){
  library(Seurat)
  library(SeuratDisk)
  library(destiny)
  library(SingleCellExperiment)
  library(ggplot2)
  library(ggthemes)
  umap = Embeddings(seu_obj,"umap")
  dmm <- DiffusionMap(umap)
  starting_cell_type <- root_type
  starting_cell_ids <- WhichCells(seu_obj, ident = starting_cell_type)
  starting_cell_index <- which(rownames(umap) %in% starting_cell_ids)
  dpt <- DPT(dmm, tips = starting_cell_index[1])
  loc <- data.frame(eigenvectors(dmm)[,1:3], time = dpt$dpt, type = rownames(umap))
  return(list(dpt = dpt, loc = loc, dmm = dmm))
}

## umap
result_umap <- DCmap_umap(seurat_obj)
result_umap$loc$classic1 <- seurat_obj$classic1
saveRDS(result_umap, './files/DCmap_umap_classic1.rds')

## plot
result <- readRDS('./files/DCmap_umap_classic1.rds')

TYPE <- plot_ly(result$loc, x = ~DC1, y = ~DC2, z = ~DC3, color = ~classic1,
                colors = col_names,
                type = 'scatter3d', mode = 'markers') %>%
  layout(title = "3D Diffusion Map",
         scene = list(xaxis = list(title = 'DC 1'),
                      yaxis = list(title = 'DC 2'),
                      zaxis = list(title = 'DC 3')))
TIME <- plot_ly(result$loc, x = ~DC1, y = ~DC2, z = ~DC3, color = ~time,
                colors = colorRamp(c(colorRampPalette(colors = c("#2166AC","#4393C3","#92C5DE","#D1E5F0"))(50),
                                     colorRampPalette(colors = c("#FDDBC7","#F4A582","#D6604D","#B2182B"))(50))),
                type = 'scatter3d', mode = 'markers') %>%
  layout(title = "3D Diffusion Map",
         scene = list(xaxis = list(title = 'DC 1'),
                      yaxis = list(title = 'DC 2'),
                      zaxis = list(title = 'DC 3')))


#### 3.Fig4c ----
library(ggridges)
library(foreign)

seurat_obj0 <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
seurat_obj <- subset(seurat_obj0, idents = c('Art.c01.DKK2', 'Art.c02.NEBL', 'Art.c03.VEGFA',  
                                             'Cap.c04.PAPSS2', 'Cap.c05.RGCC','Cap.c06.RDH10',
                                             'Cap.c07.CA2', 'Cap.c08.TMEM163', 'Cap.c09.RGS6'))
type_colors <- c("#f2ccac","#edae11","#FE870D","#e1abbc","#F88A89","#c04932","#c1d5b9","#99CD91","#009E73")
#### 3.1 tip score ----
features <- read.table('./files/tip.txt')|> pull()
seurat_obj <- AddModuleScore(seurat_obj, features = list(features), name = "PW")
meta.data <- seurat_obj@meta.data[, c("cell_type", "PW1")]
meta.data <- meta.data %>% reshape2::melt(id.vars = c("cell_type"))
colnames(meta.data)[2:3] <- c("signature", "score")
Tip <- meta.data
saveRDS(meta.data, file.path(outdir, 'tip score CAPsubtype.rds'))

#### 3.2 plot ----
ggplot(meta.data, aes(x = score, y = cell_type,fill = cell_type)) +
  geom_density_ridges(rel_min_height = 0.005)+
  scale_y_discrete(expand = c(0.01, 0)) +
  scale_x_continuous(expand = c(0.01, 0))+
  theme_ridges()+ guides(fill=FALSE) +
  theme(axis.title = element_blank(),
        panel.grid = element_blank()) + 
  scale_fill_manual(values = type_colors[1:9])+
  labs(title = 'Tip score')
ggsave(file.path(outdir2, 'Tip_Ridges.pdf'), ggplot2::last_plot(), height=4, width=5)

#### 4.Fig4d ----
#### 4.1 stalk score ----
features <- read.table('./files/stalk.txt')|> pull()
seurat_obj <- AddModuleScore(seurat_obj, features = list(features), name = "PW")
meta.data <- seurat_obj@meta.data[, c("cell_type", "PW1")]
meta.data <- meta.data %>% reshape2::melt(id.vars = c("cell_type"))
colnames(meta.data)[2:3] <- c("signature", "score")
Stalk <- meta.data
saveRDS(meta.data, file.path(outdir, 'stalk score CAPsubtype.rds'))

#### 4.2 R(S/T) ----
data <- cbind(Tip[,c(1,3)], Stalk[,c(3)])
colnames(data) <- c("cell_type", "Tip", "Stalk")
data$RTS <- (data$Tip+1)/(data$Stalk+1)
saveRDS(data, file.path(outdir, 'R_TS_value.rds'))
data <- data |> filter(cell_type %in% c('Art.c01.DKK2', 'Art.c02.NEBL', 'Art.c03.VEGFA',  
                                        'Cap.c04.PAPSS2', 'Cap.c05.RGCC','Cap.c06.RDH10',
                                        'Cap.c07.CA2', 'Cap.c08.TMEM163', 'Cap.c09.RGS6'))
saveRDS(data, file.path(outdir, 'R_TS_value_sub.rds'))

#### 4.3 Vlnplot ----
mean.score <- mean(data$RTS)
median.score <- data %>% group_by(cell_type) %>% summarise(median_score = median(RTS))
data$RTS[data$RTS>2] = 2

library(ggpubr)
p1 <- data %>%
  ggplot(aes(x = cell_type, y = RTS, fill = cell_type)) +
  geom_violin(color = "NA") +
  geom_point( data = median.score,
              aes(x = cell_type, y = median_score),
              shape = 3, size = 2) +
  geom_hline(yintercept = 1,
             color = "darkgrey", linetype = "dashed") +
  scale_y_continuous("R T/S value") + 
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.text.x = element_text(colour = 'black', angle = 30, hjust = 1)
  )+
  scale_fill_manual(values = type_colors[1:9])+NoLegend() 
ggsave(file.path(outdir2, 'RTS Boxplot noCompar.pdf'), p1, height=3.5, width=4)

#### 4.4 Featureplot ----
rts_values <- readRDS(file.path(outdir, 'R_TS_value_sub.rds')) |> pull('RTS')
scale_rts <- function(x) {
  if (x <= 1) {
    return (x / 1 * 1)
  } else {
    return ((x - 1) / 2 * 1 + 1)
  }
}
scaled_rts_values <- sapply(rts_values, scale_rts)
seurat_obj$scaled_RTS <- scaled_rts_values
low_colors <- colorRampPalette(colors = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0"))(50)
high_colors <- colorRampPalette(colors = c("#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(50)
custom_colors <- c(low_colors, high_colors)
p2 <- FeaturePlot(seurat_obj, features = "scaled_RTS",pt.size = 0.2) +
  labs(title = 'R T/S value') +
  theme(panel.grid = element_blank(), aspect.ratio = 1) +
  scale_color_gradientn(colors = custom_colors, limits = c(0.5, 1.5), oob = scales::squish)
ggsave(file.path(outdir2, 'RTS FeaturePlot.pdf'), p2, height=4.5, width=4.5)

#### 5. Fig S7b + Fig4e ----
#### 5.1 diffusionMap + Fig4e ----
col_names <- c("#f2ccac","#edae11","#FE870D","#e1abbc","#F88A89","#c04932",
               "#c1d5b9","#99CD91","#009E73","#6a73cf","#916CB2","#237BB2","#8A9D5B")
seurat_obj <- qread('~/project/PH/filter_name/results/SCT15_1_P2/files/named/classic/Cap_Art.qs')
root_type <- "Cap.c08.TMEM163"
DCmap_umap <- function(seu_obj){
  library(Seurat)
  library(SeuratDisk)
  library(destiny)
  library(SingleCellExperiment)
  library(ggplot2)
  library(ggthemes)
  umap = Embeddings(seu_obj,"umap")
  dmm <- DiffusionMap(umap)
  starting_cell_type <- root_type
  starting_cell_ids <- WhichCells(seu_obj, ident = starting_cell_type)
  starting_cell_index <- which(rownames(umap) %in% starting_cell_ids)
  dpt <- DPT(dmm, tips = starting_cell_index[1])
  loc <- data.frame(eigenvectors(dmm)[,1:3], time = dpt$dpt, type = rownames(umap))
  return(list(dpt = dpt, loc = loc, dmm = dmm))
}

result_umap <- DCmap_umap(seurat_obj)
result_umap$loc$cell_type <- seurat_obj$cell_type
saveRDS(result_umap, 'DCmap_umap_sub.rds')

result <- readRDS('DCmap_umap_sub.rds')
TYPE <- plot_ly(result$loc, x = ~DC1, y = ~DC2, z = ~DC3, color = ~cell_type,
                colors = col_names,
                type = 'scatter3d', mode = 'markers') %>%
  layout(title = "3D Diffusion Map",
         scene = list(xaxis = list(title = 'DC 1'),
                      yaxis = list(title = 'DC 2'),
                      zaxis = list(title = 'DC 3')))
TYPE

TIME <- plot_ly(result$loc, x = ~DC1, y = ~DC2, z = ~DC3, color = ~time,
                colors = colorRamp(c(colorRampPalette(colors = c("#2166AC","#4393C3","#92C5DE","#D1E5F0"))(50),
                                     colorRampPalette(colors = c("#FDDBC7","#F4A582","#D6604D","#B2182B"))(50))),
                type = 'scatter3d', mode = 'markers') %>%
  layout(title = "3D Diffusion Map",
         scene = list(xaxis = list(title = 'DC 1'),
                      yaxis = list(title = 'DC 2'),
                      zaxis = list(title = 'DC 3')))
TIME

#### 5.2 slingshot ----
# seurat_obj <- subset(seurat_obj, idents = c('Art.c01.DKK2', 'Art.c02.NEBL', 'Art.c03.VEGFA',  
#                                              'Cap.c04.PAPSS2', 'Cap.c05.RGCC',
#                                              'Cap.c07.CA2', 'Cap.c08.TMEM163'))
strace <- function(seu_obj){
  library(Seurat)
  library(RColorBrewer)
  library(slingshot)
  library(SingleCellExperiment)
  library(ggplot2)
  library(ggthemes)
  seu_obj$idents = as.character(Idents(seu_obj))
  a.sce <- as.SingleCellExperiment(seu_obj)
  umap_embeddings <- Embeddings(seu_obj, "umap")
  reducedDims(a.sce)$UMAP <- umap_embeddings
  sim <- slingshot(a.sce, clusterLabels = 'cell_type', reducedDim = 'UMAP')
  pdf("slingshot_angio.pdf", height = 7, width = 6.5)
  col = col_names
  names(col) = unique(sim$cell_type)
  plot(reducedDims(sim)$UMAP, col = col[sim$cell_type], pch = 16, asp = 1.5)
  lines(SlingshotDataSet(sim), lwd = 2, type = 'lineages', col = 'black')
  dev.off()
  return(sim)
}
strace(seurat_obj)

#### 5.3 monocle3 ----
seu_obj <- seurat_obj
mono3 <- function(seu_obj ){
  library(Seurat)
  library(monocle3)
  library(tidyverse)
  library(patchwork)
  data <- GetAssayData(seu_obj, assay = 'RNA', slot = 'counts')
  cell_metadata <- seu_obj@meta.data
  gene_annotation <- data.frame(gene_short_name = rownames(data))
  rownames(gene_annotation) <- rownames(data)
  cds <- new_cell_data_set(data,
                           cell_metadata = cell_metadata,
                           gene_metadata = gene_annotation)
  cds <- preprocess_cds(cds, num_dim = 50)
  cds <- reduce_dimension(cds,preprocess_method = "PCA") #preprocess_method默认是PCA
  
  cds.embed <- cds@int_colData$reducedDims$UMAP
  int.embed <- Embeddings(seu_obj, reduction = "umap")
  int.embed <- int.embed[rownames(cds.embed),]
  cds@int_colData$reducedDims$UMAP <- int.embed
  
  cds <- cluster_cells(cds, reduction_method = "UMAP")
  cds <- learn_graph(cds)
  get_earliest_principal_node <- function(cds, cell_phenotype="broad_cell_type", root_types){
    root_pr_nodes <- lapply(root_types, function(root_type) {
      cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
      closest_vertex <-
        cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
      closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
      root_pr_nodes <-
        igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                  (which.max(table(closest_vertex[cell_ids,]))))]
      return(root_pr_nodes)
    })
    return(unlist(root_pr_nodes))
  }
  
  cds <- order_cells(cds, root_pr_nodes= get_earliest_principal_node(cds,
                                                                     cell_phenotype="cell_type",# celltype列名
                                                                     root_types=c('Cap.c08.TMEM163')))
  
  saveRDS(cds, "monocle3_sub.rds")
}

mono3(seurat_obj)

## 绘图
cds <- readRDS("monocle3_sub.rds")
#自定义函数
get_earliest_principal_node <- function(cds, cell_phenotype="broad_cell_type", root_types){
  root_pr_nodes <- lapply(root_types, function(root_type) {
    cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)
    closest_vertex <-
      cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <-
      igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                (which.max(table(closest_vertex[cell_ids,]))))]
    return(root_pr_nodes)
  })
  return(unlist(root_pr_nodes))
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds, cell_phenotype="seurat_clusters", root_types=c("1","4")))
cds2 <- order_cells(cds, root_pr_nodes=
                      get_earliest_principal_node(cds,
                                                  cell_phenotype="cell_type",# celltype列名
                                                  root_types=c('Cap.c08.TMEM163')))
saveRDS(cds2, "monocle3_sub_cds2.rds")

p2 <- plot_cells(cds2, reduction_method="UMAP",
                 show_trajectory_graph = T, color_cells_by="pseudotime",
                 group_label_size = 7, trajectory_graph_color = 'white',
                 label_groups_by_cluster=FALSE, label_branch_points=FALSE,
                 label_leaves=FALSE, label_cell_groups = FALSE) 
p2 <- rasterize(p2, layers = "Point", dpi = 600)
ggsave('./plots/monocle3_pseudotime2.pdf', p2, height=3, width=4)

#### 6.Fig4f ----
library(ggradar)
library(AUCell)

seurat_obj <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
types <- c('Art.c01.DKK2','Cap.c05.RGCC','Cap.c08.TMEM163')
AUCell_socre <- readRDS("~/PH/results/Fig3/files/C5_AUCell_socre.rds")
colnames(AUCell_socre) <- colnames(AUCell_socre)|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()

features <- c(
  # Angiogenesis
  'Sprouting angiogenesis',
  'Regulation of sprouting angiogenesis',
  'Vasculogenesis',
  # Inflammation
  'Nik nf kappab signaling',
  'Inflammasome complex assembly',
  'Interferon gamma mediated signaling pathway',
  # Vasculature development
  'Regulation of wnt signaling pathway',
  'Notch signaling pathway',
  'Regulation of vasculature development',
  # ECM organization
  'Regulation of extracellular matrix organization',
  'Extracellular matrix assembly',
  'Cell substrate adhesion',
  # Hypoxia
  'Response to oxygen levels',
  'Regulation of cellular response to hypoxia',
  'Cellular response to oxygen levels',
  # Platelet aggregation
  'Platelet activation',
  'Regulation of platelet derived growth factor receptor signaling pathway',
  'Response to platelet derived growth factor',
  # Immune response
  'Activation of immune response',
  'Immune response regulating signaling pathway',
  'Positive regulation of leukocyte cell cell adhesion',
  # Phospholipid metabolism
  'Glycerophospholipid biosynthetic process',
  'Glycerophospholipid metabolic process',
  'Phospholipid transport',
  # TGFβ signaling
  'Regulation of cellular response to transforming growth factor beta stimulus',
  'Transforming growth factor beta receptor signaling pathway',
  'Vascular endothelial growth factor receptor signaling pathway',
  # Cell proliferation and migration
  'Endothelial cell proliferation',
  'Endothelial cell migration',
  'Blood vessel endothelial cell migration')

AUC_score_f <- AUCell_socre[,colnames(AUCell_socre) %in% features]
AUC_score_f$cell_type <- seurat_obj$cell_type
meta <- AUC_score_f |> filter(AUC_score_f$cell_type %in% c('Art.c01.DKK2','Cap.c05.RGCC','Cap.c08.TMEM163'))

score <- aggregate(meta[, 1:30], list(meta$cell_type), FUN = mean)
rownames(score) <- score$Group.1
colnames(score)[1] <- 'cell_type'
score$cell_type <- as.character(score$cell_type)
data <- score[,-1]
data <- data[,features]

path_s <- c(
  # Angiogenesis
  'Regulation of sprouting angiogenesis',
  # Cell proliferation and migration
  'Endothelial cell migration',
  # TGFβ signaling
  'Vascular endothelial growth factor receptor signaling pathway',
  # Immune response
  'Positive regulation of leukocyte cell cell adhesion',
  # Inflammation
  'Nik nf kappab signaling',
  # Phospholipid metabolism
  'Phospholipid transport',
  # Platelet aggregation
  'Platelet activation',
  # Vasculature development(NOTCH pathway)
  'Notch signaling pathway',
  # Hypoxia
  'Regulation of cellular response to hypoxia',
  # ECM organization
  'Extracellular matrix assembly'
)
data_s <- data[,path_s]

pdata1 <- cbind(data.frame(cell_type=c('Art.c01.DKK2','Cap.c05.RGCC','Cap.c08.TMEM163')),
                data_s)
for (m in 2:11) {
  pdata1[,m] <- (pdata1[,m]-min(pdata1[,m]))/(max(pdata1[,m])-min(pdata1[,m]))
}
ggradar(pdata1, background.circle.transparency = 0,
        group.colours = c('#f2ccac',"#F88A89","#99CD91"),
        values.radar = c('0','0.5','1'),
        plot.extent.x.sf = 1.4,
        base.size = 3,axis.label.size = 4,grid.label.size = 4,
        group.point.size = 2,
        group.line.width = 1,
        legend.position = 'top',
        legend.text.size = 10)+
  theme(plot.title = element_text(hjust = 0.2,size = 12))
ggsave(paste0(outdir2,'radar_pathways.pdf'),last_plot(),width = 6,height = 3.5)

#### 7.Fig4g ----
seurat_obj <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
features <- c("Art.c01.DKK2","Art.c02.NEBL","Art.c03.VEGFA")

prop.table(table(seurat_obj$cell_type, seurat_obj$sample), margin = 2)
prop <- as.data.frame(prop.table(table(seurat_obj$cell_type, seurat_obj$sample), margin = 2))
colnames(prop) <- c('cell_type','sample','prop')
meta <- data.frame(sample=seurat_obj$sample,
                   group=seurat_obj$group,
                   disease=seurat_obj$disease)
meta$group <- as.character(meta$group)
meta$disease <- as.character(meta$disease)
meta <- meta %>% distinct(sample,group,disease, .keep_all = T)

## plot
prop$group <- ''
for (i in 1:length(rownames(prop))) {
  prop$group[i] <- meta[meta$sample == prop$sample[i],]$group
}

prop$group <- factor(prop$group, levels = c('Health', 'CVD'))

pdata <- prop[prop$cell_type %in% features,]
fdata <- pdata[pdata$prop > 0,]

ggplot(fdata, aes(x=cell_type, y=prop, color=group))+
  geom_boxplot(position = position_dodge(width = 0.8), size = 1.1,
               outlier.size = 0, width = 0.4, show.legend = T)+
  stat_compare_means(data = fdata, aes(x=cell_type, y=prop, color=group),
                     method = 'wilcox.test', label = "p.signif",
                     vjust = 1)+
  scale_color_manual(values = group_colors)+
  theme_classic() + labs(x = '', y = 'Proporation')+
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        axis.title.y = element_text(color = 'black', size = 12, vjust = 1.5))
ggsave(paste0(outdir2,'prop_art_group.pdf'),last_plot(),width = 6,height = 2.5)

#### 8.Fig4j+FigS8g ----
library(SeuratDisk)
library(destiny)
library(SingleCellExperiment)
library(ggthemes)

result <- readRDS('./files/DCmap_umap_celltype.rds')
genelist <- c('DLL4','HES5','NCSTN')
cols <- c("#B3DE69", "#FB8072", "#80B1D3", "#8DD3C7", "#FCCDE5", "#377EB8", "#1B9E77")
for(i in genelist){
  data <- result$loc
  data$GENE <- seurat_obj@assays$RNA@data[i,]
  data <- data %>% mutate(rank = rank(time, ties.method = "first"))
  data$disease <- seurat_obj$disease
  data <- data %>% filter(data$disease %in% c('Normal','ACM','DCM','ICM','AMI','CS','COV'))
  sampled_data <- data %>% sample_n(size = nrow(data) / 2)
  
  p <- ggplot(sampled_data, aes(x = time, y = GENE)) +
    geom_smooth(aes(color = disease), method = "loess", se = TRUE) +
    labs(title = i, x = "Pseudotime", y = "Expression") +
    theme_classic() +
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          axis.text = element_text(color = "black", size = 8),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(color = "black", size = 8))+
    scale_color_manual(values = cols)
  ggsave(file.path(outdir2,paste0(i, '_expr_perse_disease.pdf')), p, width = 3, height = 3)
}

#### 9.Fig4k ----
seurat_obj <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
seurat_obj <- subset(seurat_obj, classic1 == 'ArtECs')
sub_obj <- seurat_obj[, seurat_obj[["disease"]] == 'ICM']
selected_genes <- c('EFNB2','DLL4')
expr <- GetAssayData(sub_obj, slot = "data", assay = "RNA")[selected_genes, ] %>% as.data.frame()
expr <- as.data.frame(t(expr))

p <- expr %>%
  filter(DLL4 != 0, EFNB2 != 0) %>%  # 去除X或Y中含有0值的行
  ggplot(aes(x = DLL4, y = EFNB2)) +
  ggrastr::rasterise(ggpointdensity::geom_pointdensity(size = 0.5), dpi = 600) +
  geom_smooth(method = "lm", formula = y ~ x,
              color = "#6b76ae", fill = "#e5a323",
              size = 0.5, alpha = 0.2) +
  theme_cat() +
  theme(aspect.ratio = 1, legend.margin = margin(l= -8),
        plot.title = element_text(size = 12),
        axis.title = element_text(face = "italic", size = 12)) +
  guides(color = guide_colorbar(
    frame.colour = "black", frame.linewidth = 0.5,
    ticks.colour = "black", title = "Density")) +
  scale_color_viridis_c() + labs(title = 'ICM')+
  ggpubr::stat_cor(method = "pearson")  # 添加相关性
ggsave(file.path(outdir2, paste0("Cor_DLL4_Artc_EFNB2_ICM.pdf")),
       p3, height=3.5, width=3.5)

