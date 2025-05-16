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

path <- '~/PH/results/Fig5/'
setwd(path)
outdir <- './files/'
outdir2 <- './plots/'

#### 1.Fig5a ----
library(plotly)
library(Seurat)
library(qs)
library(dplyr)

seurat_obj <- qread("~/PH/results/named/files/Pan_all.qs")

seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
head(VariableFeatures(seurat_obj),10)

CM <- read.table('~/PH/results/named/files/CCA-SCT_2000rm_features.txt') |> pull()
FB <- read.table('~/PH/results/named/files/CCA-SCT_2000rm_features.txt') |> pull()
LYM <- read.table('~/PH/results/named/files/CCA-SCT_2000rm_features.txt') |> pull()
MYE <- read.table('~/PH/results/named/files/CCA-SCT_2000rm_features.txt') |> pull()
EC <- read.table('~/PH/results/named/files/sample_15_1800_1800rm_features.txt') |> pull()
HVGs <- unique(Reduce(intersect, list(CM,FB,LYM,MYE,EC)))

VariableFeatures(seurat_obj) <- HVGs
head(VariableFeatures(seurat_obj),10)

#### 1.1 harmony ----
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- seurat_obj |>
  harmony::RunHarmony(group.by.vars =c("orig.ident"),
                      dims.use = 1:20,
                      plot_convergence = TRUE)

seurat_obj <- RunUMAP(seurat_obj ,reduction = "harmony", dims = 1:20, n.components = 3L)
saveRDS(seurat_obj,"./files/seurat_all.rds")

#### 1.2 umap ----
head(seurat_obj[["umap"]]@cell.embeddings)
umap_1 <- seurat_obj[["umap"]]@cell.embeddings[,1]
umap_2 <- seurat_obj[["umap"]]@cell.embeddings[,2]
umap_3 <- seurat_obj[["umap"]]@cell.embeddings[,3]
plot.data <- FetchData(object = seurat_obj, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "cell_type"))
plot.data$label <- paste(rownames(plot.data))

# plot
type_names <- c(
  "Art.c01.DKK2", "Art.c02.NEBL", "Art.c03.VEGFA", "Cap.c04.PAPSS2", "Cap.c05.RGCC", "Cap.c06.RDH10",
  "Cap.c07.CA2", "Cap.c08.TMEM163", "Cap.c09.RGS6", "Ven.c10.VCAM1", "Ven.c11.ACTA2", "LEC.c12.FLT4", "Endo.c13.NPR3", # EC
  'Mono_CD14', 'Mono_CD16',  'Mac_SPP1', 'Mac_IFI44L', 'Mac_LYVE1',
  'Mac_NLRP3', 'Mac_THBS1', 'Mac_TMEM163', 'Mac_FSCN1','DC_HLA_DQA1', # Mye
  'BC_MS4A1', 'CD4_IL7R', 'CD4_THEMIS', 'CD4_IL2RA', 'CD8_GZMK', 'CD8_GZMH', 'NK_GNLY', 'NK_NCAM1', # Lym
  'FB_APOD', 'FB_PCOLCE2', 'FB_OSMR', 'FB_ACTA2', 'FB_POSTN', 'FB_NID2',  # FB
  'CM_ROR2', 'CM_MYH7', 'CM_PRELID2', 'CM_CNN1', 'CM_CRYAB'  # CM
)
type_colors2 <- c(colorRampPalette(c("#CD3F2C","#FFAC8C"))(3), colorRampPalette(c("#FF85A7","#FDE2E1"))(3),
                  colorRampPalette(c("#A2416E","#D8BFD8"))(3), colorRampPalette(c("#FFC75F","#F9F871"))(4), # EC
                  colorRampPalette(c("#845EC2","#CEE7F1"))(5), colorRampPalette(c("#9B89B3","#FCF8FF"))(5), # Mye
                  colorRampPalette(c("#5E7AB7","#74CCED"))(4), colorRampPalette(c("#3596B5","#C4F9FF"))(4), # Lym
                  colorRampPalette(c("#3A8200","#D9EDDF"))(6), # FB
                  "#7E5E12", "#B5AA99", colorRampPalette(c("#AF5C00","#F9E282"))(3) # CM
)
plot_ly(data = plot.data, 
        x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3, 
        color = ~cell_type, colors = type_colors2,
        type = "scatter3d", mode = "markers", 
        marker = list(size = 4, width=2),
        text=~label, hoverinfo="text")

#### 2.Fig5b ----
#### 2.1 cellchat ----
#### 2.1.1 Normal ----
library(CellChat)
library(patchwork)
library(Seurat)
library(qs)
library(future)
library(ggpubr)
options(future.globals.maxSize = 50 * 1024^3)

case <- "Normal_"
seurat_obj <- qread("~/PH/results/named/files/Pan_Normal_3000.qs")

data.input <- seurat_obj[["RNA"]]@data 
labels <- Idents(seurat_obj)
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)

future::plan("multisession", workers = 10) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

unique(cellchat@idents)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file = paste0(output.dir, case, "cellchat.rds"))

#### 2.1.2 CVD ----
case <- "CVD_"
seurat_obj <- qread("~/PH/results/named/files/Pan_CVD_3000.qs")

data.input <- seurat_obj[["RNA"]]@data 
labels <- Idents(seurat_obj)
meta <- data.frame(labels = labels, row.names = names(labels))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)

future::plan("multisession", workers = 10) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

unique(cellchat@idents)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file = paste0(output.dir, case, "cellchat.rds"))

#### 2.2 barplot_EC ----
cellchat.1 <- readRDS("./files/Normal_cellchat.rds")
cellchat.2 <- readRDS("./files/CVD_cellchat.rds")
object.list <- list(Normal = cellchat.1, CVD = cellchat.2)

classic1_colors <- c("#edd064","#f2ccac","#a1d5b9","#57C3F3",'#6778AE')
## 2.2.1 Bar_count ----
CVD_count <- cellchat.2@net[["count"]]
colnames(CVD_count)
CVD_count2 <- as.data.frame(CVD_count[1:11, 14:42])
CVD_count2$count <- rowSums(CVD_count2)
CVD_count2$Name <- rownames(CVD_count2)
CVD_count2$group <- c(rep('Arterial',3),rep('Capillary',6),rep('Venous',2))
p1 <- ggbarplot(CVD_count2, x="Name", y="count", fill = "group", color = "white", width = 0.5,
                palette = "aaas", legend = "right",  
                sort.by.groups=TRUE )+theme(axis.title.x = element_blank(),
                                            axis.text.x = element_blank(), 
                                            axis.ticks.x = element_blank(), 
                                            axis.line.x = element_blank(),axis.text.y = element_text(size = 6))+
  scale_fill_manual(values = classic1_colors[1:3])
## 2.2.2 Bar_weight ----
CVD_weight <- cellchat.2@net[["weight"]]
colnames(CVD_weight)
CVD_weight2 <- as.data.frame(CVD_weight[1:11, 14:42])
CVD_weight2$Weight <- rowSums(CVD_weight2)
CVD_weight2$Name <- rownames(CVD_weight2)
CVD_weight2$group <- c(rep('Arterial',3),rep('Capillary',6),rep('Venous',2))
p2 <- ggbarplot(CVD_weight2, x="Name", y="Weight", fill = "group", color = "white", width = 0.5,
                palette = "aaas", legend = "right",   
                sort.by.groups=TRUE )+theme(axis.title.x = element_blank(),
                                            axis.text.x = element_blank(), 
                                            axis.ticks.x = element_blank(), 
                                            axis.line.x = element_blank(),axis.text.y = element_text(size = 6))+
  scale_fill_manual(values = classic1_colors[1:3])
ggsave(file.path(outdir2,'Count_Weight_Barplot.pdf'), p1/p2, height=2, width=7)

#### 2.3 barplot_T4 ----
## 2.3.1 Bar_conut ----
CVD_count <- cellchat.2@net[["count"]]
colnames(CVD_count)
CVD_count2 <- as.data.frame(CVD_count[14:42, 1:11])
CVD_count2$count <- rowSums(CVD_count2)
CVD_count2$Name <- rownames(CVD_count2)
CVD_count2$group <- c(rep('Myeloid',10),rep('Lymphoid',8),rep('Fibroblast',6),rep('Cardiomyocyte', 5))
CVD_count2$group <- factor(CVD_count2$group, levels = c('Myeloid', 'Lymphoid', 'Fibroblast', 'Cardiomyocyte'))
p1 <- ggbarplot(CVD_count2, x="Name", y="count", fill = "group", color = "white", width = 0.7,
                palette = "aaas", legend = "right", 
                sort.by.groups=TRUE )+theme(axis.title.x = element_blank(),
                                            axis.text.x = element_blank(), 
                                            axis.ticks.x = element_blank(), 
                                            axis.line.x = element_blank(),axis.text.y = element_text(size = 6))+
  scale_fill_manual(values = c("#EFC000FF", "#B86886FF", "#0073C2FF", "#DD531CFF"))

## 2.3.2 Bar_weight ----
CVD_weight <- cellchat.2@net[["weight"]]
colnames(CVD_weight)
CVD_weight2 <- as.data.frame(CVD_weight[14:42, 1:11])
CVD_weight2$Weight <- rowSums(CVD_weight2)
CVD_weight2$Name <- rownames(CVD_weight2)
CVD_weight2$group <- c(rep('Myeloid',10),rep('Lymphoid',8),rep('Fibroblast',6),rep('Cardiomyocyte', 5))
CVD_weight2$group <- factor(CVD_weight2$group, levels = c('Myeloid', 'Lymphoid', 'Fibroblast', 'Cardiomyocyte'))
p2 <- ggbarplot(CVD_weight2, x="Name", y="Weight", fill = "group", color = "white", width = 0.7,
                palette = "aaas", legend = "right",
                sort.by.groups=TRUE )+theme(axis.title.x = element_blank(),
                                            axis.text.x = element_blank(), 
                                            axis.ticks.x = element_blank(), 
                                            axis.line.x = element_blank(),axis.text.y = element_text(size = 6))+
  scale_fill_manual(values = c("#EFC000FF", "#B86886FF", "#0073C2FF", "#DD531CFF"))
ggsave(file.path(outdir2,'Count_Weight_Barplot_T4.pdf'), p1/p2, height=3, width=15)

#### 2.4 heatmap ----
library(pheatmap)
count <- cellchat.2@net[["count"]]
count <- as.data.frame(t(count[1:11, 14:42]))
wright <- round(as.data.frame(t(cellchat.2@net[["weight"]][1:11, 14:42])),2)

pdf(file.path(outdir2, 'Heatmap_count.pdf'),width = 6,height = 6)
pheatmap(count, row_names_side = 'left',
         color = colorRampPalette(c("#21295C","white","#CD3333"))(1000),
         cellwidth = 15, cellheight = 8, border_color = "white", fontsize = 8, 
         show_rownames = T,show_colnames = T,cluster_rows = F, cluster_cols = F,
         scale = "row", treeheight_row = 10, treeheight_col =10,
         display_numbers = wright, fontsize_number = 4,
         gaps_row=c(10,18,24), clustering_method = "average") # 
dev.off()

#### 2.5 dotplot ----
seurat_obj <- qread('~/PH/results/named/files/T4_all.qs')
p <- DotPlot(seurat_obj, features= c('SPP1', 'CXCL2', 'CXCL10', 'CXCL12', 'IL1B', 'VEGFA', 'ANGPT1', 'ANGPT2'))+
  theme(axis.text.x = element_text(face = 'italic',angle=90,hjust=1,vjust=0.5,size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_color_distiller(palette= "RdBu")
ggsave(paste0(outdir2, 'LR_expr_dotplot.pdf'), ggplot2::last_plot(), height=7, width =3.8)

#### 3.Fig5c ----
endotypes <- c("Art.c01.DKK2","Art.c02.NEBL","Art.c03.VEGFA","Cap.c04.PAPSS2","Cap.c05.RGCC","Cap.c06.RDH10",
               "Cap.c07.CA2","Cap.c08.TMEM163","Cap.c09.RGS6","Ven.c10.VCAM1","Ven.c11.ACTA2","LEC.c12.FLT4","Endo.c13.NPR3")

#### 3.1 line ----
cellchat <- readRDS("./files/CVD_cellchat.rds") 
sources = c(14:42)
targets = c(1:11)
p <- netVisual_bubble(cellchat, sources.use = sources, targets.use = targets,
                      signaling = c('CCL', 'CXCL', "IL6", "VEGF", "ANGPT","SPP1" ),
                      angle.x = 90)
LR_data <- p[["data"]][,1:4]
LR_data <- LR_data[!is.na(LR_data$ligand), ]
table(LR_data$ligand)
table(LR_data$receptor)

LR_data$receptor <- gsub('FLT1_KDR', 'KDR', LR_data$receptor)
LR_data$receptor <- gsub('FLT4_KDR', 'KDR', LR_data$receptor)
LR_data$receptor <- gsub('ITGA5_ITGB1', 'ITGB1', LR_data$receptor)
LR_data$receptor <- gsub('ITGA8_ITGB1', 'ITGB1', LR_data$receptor)
LR_data$receptor <- gsub('ITGA9_ITGB1', 'ITGB1', LR_data$receptor)
LR_data$receptor <- gsub('ITGAV_ITGB1', 'ITGB1', LR_data$receptor)
LR_data$receptor <- gsub('ITGAV_ITGB5', 'ITGB5', LR_data$receptor)
write.table(LR_data, file = file.path(outdir, "LR_data.csv"),
            quote = F, sep = ",", row.names = T)

#
data <- data.frame(sou=LR_data$ligand,
                   x1=rep(2,length(LR_data$ligand)),
                   net_y1='',
                   tar=LR_data$receptor,
                   x2=rep(3,length(LR_data$receptor)),
                   net_y2='') |> unique()
LR_data_name <- cbind(unique(data$sou), c(1:length(unique(data$sou)))) |> as.data.frame()
colnames(LR_data_name) <- c('sou', 'num')
match_indices <- match(data$sou, LR_data_name$sou)
data$net_y1 <- LR_data_name$num[match_indices] |> as.numeric()

LR_data_name <- cbind(unique(data$tar), c(1:length(unique(data$tar)))) |> as.data.frame()
colnames(LR_data_name) <- c('tar', 'num')
match_indices <- match(data$tar, LR_data_name$tar)
data$net_y2 <- LR_data_name$num[match_indices] |> as.numeric()
write.table(data, file = file.path(outdir, "data_dotplot.csv"),
            quote = F, sep = ",", row.names = T)

add <- data.frame(sou=c("CXCL2","IL33","IL33"),
                  x1=c(2,2,2),
                  net_y1=c(11,12,12),
                  tar=c("LPAR1","IL1RAP","IL1R1"),
                  x2=c(3,3,3),
                  net_y2=c(10,11,12))
adata <- rbind(data,add)

## plot
p2 <- ggplot(adata)+
  geom_segment(aes(x1,net_y1,xend=x2,yend=net_y2),
               size=0.5,color=c('black','black','black','black','black',
                                'black','black','black','black','black',
                                'black','black','black','black','black',
                                'black','black','black'))+
  geom_point(aes(x=x1,y=net_y1),size=5, fill="#62d2a2", color="#62d2a2",
             stroke=1, shape = 21)+
  geom_point(aes(x=x2,y=net_y2),size=5, fill="#66bfbf", color="#66bfbf",
             stroke=1, shape = 21)+
  scale_y_continuous(limits = c(1, 12),expand = expansion(add=c(0.5,0.7)))+
  scale_x_continuous(expand = expansion(0,0.1))+
  theme_void()
ggsave(file.path(outdir2, 'mid_line2.pdf'), p2, width = 2, height = 4.5)

#### 3.2 Dotplot ----
## 3.2.1 left ----
features_1 <- c("ANGPT1","ANGPT2","CCL21","CCL8","CXCL12","PGF","SPP1","VEGFA","VEGFB","VEGFC","CXCL2","IL33")
seurat_obj <- qread("~/PH/results/named/files/Pan_all.qs")
sub_obj <- seurat_obj[, !(seurat_obj$cell_type %in% endotypes)]
p1 <- DotPlot(object = sub_obj, features = features_1, assay = "RNA") +
  guides(color = guide_colorbar(title = 'Average  Expression')) +
  scale_color_gradientn(colours = c('white',"grey","black"))+
  coord_flip() +
  theme(axis.title = element_blank(),
        legend.position = "left",
        plot.margin = unit(c(0.5,0,0.5,0.5), 'cm'),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.2),
        panel.border = element_rect(color = "black",fill=NA, size = 0.1),
        panel.grid.major = element_line(size = 0.1,color = 'darkgrey',linetype = 2),
        axis.text.x = element_text(angle = 60, vjust = 0, hjust = 0)) +
  scale_x_discrete(position = "top") + scale_y_discrete(position = "right")
p1
ggsave(file.path(outdir2, 'left_target2.pdf'), p1, width = 10, height = 4.5)

## 3.2.2 right ----
features_2 <- c("TEK","ITGB1","ACKR4","ACKR3","FLT1","ITGB5","CD44","KDR","FLT4","LPAR1","IL1RAP","IL1R1")
seurat_obj <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
seurat_obj <- subset(seurat_obj, cell_type != "Endo.c13.NPR3")
p3 <- DotPlot(object = seurat_obj, features = features_2, assay = "RNA") +
  guides(color = guide_colorbar(title = 'Average  Expression')) +
  scale_color_gradientn(colours = c('white',"grey","black"))+
  coord_flip() +
  theme(axis.title = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0), 'cm'),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.2),
        panel.border = element_rect(color = "black",fill=NA, size = 0.1),
        panel.grid.major = element_line(size = 0.1,color = 'darkgrey',linetype = 2),
        axis.text.x = element_text(angle = 60, vjust = 0, hjust=0)) +
  scale_y_discrete(position = "right")  + 
  NoLegend()
ggsave(file.path(outdir2, 'right_source2.pdf'), p3, width = 5, height = 4.5)

#### 4.FigS10a/b ----
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat@meta$labels <- factor(cellchat@meta$labels, levels = type_names)

pdf(file = paste0(outdir2, "DiffInteraction_heatmap_s.pdf"),width= 5,height= 5)
netVisual_heatmap(cellchat, 
                  sources.use = c('Mono_CD14', 'Mono_CD16',  'Mac_SPP1', 'Mac_IFI44L', 'Mac_LYVE1',
                                  'Mac_NLRP3', 'Mac_THBS1', 'Mac_TMEM163', 'Mac_FSCN1','DC_HLA_DQA1',
                                  'BC_MS4A1', 'CD4_IL7R', 'CD4_THEMIS', 'CD4_IL2RA', 'CD8_GZMK', 'CD8_GZMH', 'NK_GNLY', 'NK_NCAM1'),
                  targets.use = c("Art.c01.DKK2", "Art.c02.NEBL", "Art.c03.VEGFA", "Cap.c04.PAPSS2", "Cap.c05.RGCC", "Cap.c06.RDH10",
                                  "Cap.c07.CA2", "Cap.c08.TMEM163", "Cap.c09.RGS6", "Ven.c10.VCAM1", "Ven.c11.ACTA2"),
                  color.use = type_colors, measure = "count",width = 1, height = 1, remove.isolate = T)
netVisual_heatmap(cellchat, 
                  sources.use = c("Art.c01.DKK2", "Art.c02.NEBL", "Art.c03.VEGFA", "Cap.c04.PAPSS2", "Cap.c05.RGCC", "Cap.c06.RDH10",
                                  "Cap.c07.CA2", "Cap.c08.TMEM163", "Cap.c09.RGS6", "Ven.c10.VCAM1", "Ven.c11.ACTA2"),
                  targets.use = c('FB_APOD', 'FB_PCOLCE2', 'FB_OSMR', 'FB_ACTA2', 'FB_POSTN', 'FB_NID2',
                                  'CM_ROR2', 'CM_MYH7', 'CM_PRELID2', 'CM_CNN1', 'CM_CRYAB'),
                  color.use = type_colors, measure = "weight", width = 1, height = 1, remove.isolate = T)
dev.off()

#### 5.FigS10d ----
cellchat.1 <- netAnalysis_computeCentrality(cellchat.1)
cellchat.2 <- netAnalysis_computeCentrality(cellchat.2)
object.list <- list(Normal = cellchat.1, CVD = cellchat.2)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
##
pdf(file = paste0(outdir2, 'Dotplot_EC_Inflam_Angio.pdf'),width=6, height= 10)
sources = c(16:23)
targets = c(10)
## 
gg1 <- netVisual_bubble(cellchat, sources.use = sources, targets.use = targets,  
                        comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in CVD", 
                        signaling = c("CCL","CXCL", 'IL16', 'CDH', 'CDH5', 'IL1', 'IL6',
                                      'CX3C', "VISFATIN", "ANGPT", "THBS"), 
                        angle.x = 45, remove.isolate = F, vjust.x = 2)
##  
sources = c(16:23)
targets = c(8)
gg2 <- netVisual_bubble(cellchat, sources.use = sources, targets.use = targets,  
                        comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in CVD", 
                        signaling = c("VEGF", 'CDH5','ANGPT', 'NOTCH'), 
                        angle.x = 45, remove.isolate = F, vjust.x = 2)
gg2
gg1 / gg2
dev.off()

#### 6.FigS10f ----
library(nichenetr)
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(circlize)
library(qs)

#### 6.1 nichenet ----
seurat_obj <- qread("~/PH/results/named/files/Pan_all.qs")
Idents(seurat_obj) <- "cell_type"
clusters <- seurat_obj$cell_type

## load data
ligand_target_matrix = readRDS("./files/nichenet/ligand_target_matrix.rds")
lr_network = readRDS("./files/nichenet/lr_network.rds")
weighted_networks = readRDS("./files/nichenet/weighted_networks.rds")

## c10
sender_celltypes = c("Mac_SPP1")
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = seurat_obj, 
  receiver ='Ven.c10.VCAM1',
  condition_colname = "group", condition_oi = "CVD", condition_reference = "Health", 
  sender = sender_celltypes, 
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, 
  weighted_networks = weighted_networks, 
  top_n_ligands = 30,
  top_n_targets = 50,
  verbose = TRUE)
saveRDS(nichenet_output,file="./c10/SPP1/nichenet_output.rds")

## c08
sender_celltypes = c("Mac_THBS1")
nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = seurat_obj, 
  receiver ='Cap.c08.TMEM163',
  condition_colname = "group", condition_oi = "CVD", condition_reference = "Health", 
  sender = sender_celltypes, 
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, 
  weighted_networks = weighted_networks, 
  top_n_ligands = 30,
  top_n_targets = 50,
  verbose = TRUE)
saveRDS(nichenet_output,file="./c08/Mac_THBS1/nichenet_output.rds")

#### 6.2 Inflammation ----
nichenet_output <- readRDS("./c10/SPP1/nichenet_output.rds")
infgenes1 <- rev(c("IL33","CXCL2"))
infgenes2 <- c("CXCL12","CCL8","CCL21")
infgenes <- c(infgenes1,infgenes2)
active_ligand_target_links_df <- infgenes2 %>%
  lapply(get_weighted_ligand_target_links,
         geneset = nichenet_output[["geneset_oi"]],
         ligand_target_matrix = ligand_target_matrix,
         n = 50) %>%
  bind_rows() %>% drop_na()
active_ligand_target_links_df_circos <- nichenet_output$ligand_target_df %>%
  dplyr::filter(ligand %in% infgenes1)
active_ligand_target_links_df_circos <- rbind(active_ligand_target_links_df_circos,active_ligand_target_links_df)

ligands_to_remove <- setdiff(nichenet_output$ligand_target_df$ligand %>% unique(),
                             active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove <- setdiff(nichenet_output$ligand_target_df$target %>% unique(),
                             active_ligand_target_links_df_circos$target %>% unique())

circos_links <- active_ligand_target_links_df_circos %>%
  dplyr::filter(!target %in% targets_to_remove & !ligand %in% ligands_to_remove)

grid_col_ligand <- c("CXCL2" = "skyblue", "IL33" = "forestgreen", "CXCL12" = "purple","CCL8"="#ffc38d","CCL21"="#f4c40f")
targetgenes <- unique(circos_links$target)
grid_col_target <- setNames(rep("#ee6470", length(targetgenes)), targetgenes)

ligand_color <- circos_links %>%
  distinct(ligand, color_ligand_type = grid_col_ligand[ligand])
grid_ligand_color <- setNames(ligand_color$color_ligand_type, ligand_color$ligand)
target_color <- circos_links %>%
  distinct(target, color_target_type = grid_col_target[target])
grid_target_color <- setNames(target_color$color_target_type, target_color$target)

grid_col <- c(grid_ligand_color, grid_target_color)
transparency <- circos_links %>%
  mutate(weight = (weight - min(weight)) / (max(weight) - min(weight))) %>%
  mutate(transparency = 1 - weight) %>%
  .$transparency

target_order <- circos_links$target %>% unique()
ligand_order <- infgenes %>% intersect(circos_links$ligand)
# circos
pdf("./plots/inf_lt_circos2.pdf", width = 10, height = 10)
chordDiagram(circos_links, directional = 1, order = c(ligand_order, target_order),
             link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,
             transparency = transparency, diffHeight = 0.005,
             direction.type = c("diffHeight", "arrows"), link.arr.type = "big.arrow",
             annotationTrack = "grid", preAllocateTracks = list(track.height = 0.075))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA)

circos.clear()
dev.off()

#### 6.3 Angiogenesis ----
nichenet_output <- readRDS("./c08/Mac_THBS1/nichenet_output.rds")
#####关注配体活性-----
angiogenes <- rev(c("ANGPT1","VEGFA","ANGPT2"))
active_ligand_target_links_df_circos <- nichenet_output$ligand_target_df %>%
  filter(ligand %in% angiogenes)
ligands_to_remove <- setdiff(nichenet_output$ligand_target_df$ligand %>% unique(),
                             active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove <- setdiff(nichenet_output$ligand_target_df$target %>% unique(),
                             active_ligand_target_links_df_circos$target %>% unique())
circos_links <- active_ligand_target_links_df_circos %>%
  filter(!target %in% targets_to_remove & !ligand %in% ligands_to_remove)

grid_col_ligand <- c("ANGPT1" = "skyblue", "VEGFA" = "forestgreen", "ANGPT2" = "purple")
targetgenes <- unique(circos_links$target)
grid_col_target <- setNames(rep("#ee6470", length(targetgenes)), targetgenes) 

ligand_color <- circos_links %>%
  distinct(ligand, color_ligand_type = grid_col_ligand[ligand])
grid_ligand_color <- setNames(ligand_color$color_ligand_type, ligand_color$ligand)
target_color <- circos_links %>%
  distinct(target, color_target_type = grid_col_target[target])
grid_target_color <- setNames(target_color$color_target_type, target_color$target)

grid_col <- c(grid_ligand_color, grid_target_color)
transparency <- circos_links %>%
  mutate(weight = (weight - min(weight)) / (max(weight) - min(weight))) %>%
  mutate(transparency = 1 - weight) %>%
  .$transparency

target_order <- circos_links$target %>% unique()
ligand_order <- angiogenes %>% intersect(circos_links$ligand)
# circos
pdf("./plots/angio_lt_circos.pdf", width = 10, height = 10)
chordDiagram(circos_links, directional = 1, order = c(ligand_order, target_order),
             link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,
             transparency = transparency, diffHeight = 0.005,
             direction.type = c("diffHeight", "arrows"), link.arr.type = "big.arrow",
             annotationTrack = "grid", preAllocateTracks = list(track.height = 0.075))
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA)
circos.clear()
dev.off()

