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
  "Cap.c07.CA2", "Cap.c08.TMEM163", "Cap.c09.RGS6", "Ven.c10.VCAM1", "Ven.c11.TMSB4X", "LEC.c12.FLT4", "Endo.c13.NPR3", # EC
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

#### 3.Fig5c ----
endotypes <- c("Art.c01.DKK2", "Art.c02.NEBL", "Art.c03.VEGFA", "Cap.c04.PAPSS2", "Cap.c05.RGCC", "Cap.c06.RDH10",
               "Cap.c07.CA2", "Cap.c08.TMEM163", "Cap.c09.RGS6", "Ven.c10.VCAM1", "Ven.c11.TMSB4X", "LEC.c12.FLT4")

#### 3.1 line ----
sources = setdiff(type_names,endotypes)[2:39]
targets = endotypes
p <- netVisual_bubble(cellchat.2, sources.use = sources, targets.use = targets,
                      signaling = c("CCL", "CXCL", "VEGF", "ANGPT", "SPP1", "NOTCH"),
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

add <- data.frame(sou=c("DLL4"),
                  x1=c(2),
                  net_y1=c(8),
                  tar=c("NOTCH1"),
                  x2=c(3),
                  net_y2=c(5))
adata <- rbind(data, add)
rownames(adata) <- NULL

sou_sort <- rev(c('CXCL12', "CCL21", "CCL8", 'SPP1', "DLK1", "PGF", 
                  'JAG1', 'DLL1', 'DLL4', 'VEGFA', "VEGFB", "VEGFC", 'ANGPT1', 'ANGPT2'))
sou_all <- unique(adata$sou)
sou_sort_valid <- sou_sort[sou_sort %in% sou_all]
sou_remaining <- setdiff(sou_all, sou_sort_valid)
sou_final <- c(sou_sort_valid, sou_remaining)
sou_map <- data.frame(
  sou = sou_final,
  net_y1_new = seq_along(sou_final))
adata_new <- adata %>%
  left_join(sou_map, by = "sou")
tar_map <- adata %>%
  distinct(tar, net_y2) %>%
  arrange(net_y2) %>%
  mutate(net_y2_new = row_number())
adata_new <- adata_new %>%
  left_join(tar_map %>% select(tar, net_y2_new), by = "tar")
adata_new <- adata_new %>%
  mutate(
    net_y1 = net_y1_new,
    net_y2 = net_y2_new) %>%
  select(-net_y1_new, -net_y2_new)
adata_final <- adata_new %>% 
  arrange(net_y1) %>% 
  mutate(tar = factor(tar, levels = unique(tar))) %>% 
  mutate(net_y2 = as.integer(tar)) %>% 
  ungroup() 

adata_final
write.table(adata_final, file = file.path(outdir, "data_dotplot_end.csv"),
            quote = F, sep = ",", row.names = T)

## plot
p2 <- ggplot(adata_final)+
  geom_segment(aes(x1,net_y1,xend=x2,yend=net_y2),
               size=0.5,color=c(rep('black', nrow(adata))))+
  geom_point(aes(x=x1,y=net_y1),size=7, fill="#62d2a2", color="#62d2a2",
             stroke=1, shape = 21)+
  geom_point(aes(x=x2,y=net_y2),size=7, fill="#66bfbf", color="#66bfbf",
             stroke=1, shape = 21)+
  scale_y_continuous(limits = c(1, max(adata$net_y1)),expand = expansion(add=c(0.5,0.7)))+
  scale_x_continuous(expand = expansion(0,0.1))+
  theme_void()
p2
ggsave(file.path(outdir2, 'mid_line.pdf'), p2, width = 2, height = 4.5)

#### 3.2 Dotplot ----
adata_final <- read.csv(file.path(outdir, "data_dotplot_end.csv"))

## 3.2.1 left ----
seurat_obj1 <- qread('~/project/PH/alltype/merge/results/files/T5_ori_merged_normalized.qs')
features_1 <- sou_sort
p1 <- DotPlot(seurat_obj1, features = features_1, assay = "RNA") +
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
        axis.text.y = element_text(face = 'italic'),
        axis.text.x = element_text(angle = 60, vjust = 0, hjust = 0)) +
  scale_x_discrete(position = "top") + scale_y_discrete(position = "right")
p1
ggsave(file.path(outdir2, 'left_target.pdf'), p1, width = 10, height = 4.5)

## 3.2.2 right ----
features_2 <- unique(adata_final$tar)
seurat_obj2 <- qs::qread("~/project/PH/filter_name/results/SCT15_1_P2/files/named_new/Seurat_202504.qs")
seurat_obj2 <- subset(seurat_obj2, cell_type != "Endo.c13.NPR3")
table(seurat_obj2$cell_type)
p3 <- DotPlot(seurat_obj2, features = features_2, assay = "RNA") +
  guides(color = guide_colorbar(title = 'Average  Expression')) +
  scale_color_gradientn(colours = c('white',"grey","black"))+
  coord_flip() +
  theme(axis.title = element_blank(),
        plot.margin = unit(c(0.5,0.5,0.5,0), 'cm'),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.2),
        panel.border = element_rect(color = "black",fill=NA, size = 0.1),
        panel.grid.major = element_line(size = 0.1,color = 'darkgrey',linetype = 2),
        axis.text.y = element_text(face = 'italic'),
        axis.text.x = element_text(angle = 60, vjust = 0, hjust=0)) +
  scale_y_discrete(position = "right")  + 
  NoLegend()
p3
ggsave(file.path(outdir2, 'right_source.pdf'), p3, width = 3.5, height = 4)



