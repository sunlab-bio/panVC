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
strace <- function(seu_obj){
  library(RColorBrewer)
  library(slingshot)
  library(SingleCellExperiment)
  library(ggthemes)
  seu_obj$idents = as.character(Idents(seu_obj))
  a.sce <- as.SingleCellExperiment(seu_obj)
  umap_embeddings <- Embeddings(seu_obj, "umap")
  reducedDims(a.sce)$UMAP <- umap_embeddings
  sim <- slingshot(a.sce, clusterLabels = 'cell_type', reducedDim = 'UMAP')
  pdf("slingshot.pdf", height = 6, width = 6)
  col = col_names
  names(col) = unique(sim$cell_type)
  plot(reducedDims(sim)$UMAP, col = col[sim$cell_type], pch = 16, asp = 1)
  lines(SlingshotDataSet(sim), lwd = 2, type = 'lineages', col = 'black')
  dev.off()
  return(sim)
}
strace(seurat_obj)

#### 2.2 DiffusionMap ----
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
                                        'Cap.c04.PAPSS2', 'Cap.c05.RGCC', 'Cap.c07.CA2', 'Cap.c08.TMEM163'))
saveRDS(data, file.path(outdir, 'R_TS_value_sub.rds'))

#### 4.3 Vlnplot ----
mean.score <- mean(data$RTS)
median.score <- data %>% group_by(cell_type) %>% summarise(median_score = median(RTS))
data$RTS[data$RTS>2] = 2

library(ggpubr)
compar <- list(c('Cap.c08.TNFRSF4', 'Art.c01.DKK2'),
               c('Cap.c08.TNFRSF4', 'Cap.c04.PAPSS2'))
p1 <- data %>%
  ggplot(aes(x = cell_type, y = RTS, fill = cell_type)) +
  geom_violin(color = "NA") +
  geom_point(data = median.score, aes(x = cell_type, y = median_score), shape = 3, size = 2) +
  geom_hline(yintercept = 1, color = "darkgrey", linetype = "dashed") +
  scale_y_continuous("R T/S value") + 
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1)
  )+
  scale_fill_manual(values = type_colors[1:9])+
  stat_compare_means(comparisons = compar, label = "p.signif", vjust = 2)
ggsave(file.path(outdir2, 'RTS Boxplot.pdf'), p1, height=3.5, width=4)

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
seurat_obj <- subset(seurat_obj, idents = c('Art.c01.DKK2', 'Art.c02.NEBL', 'Art.c03.VEGFA', 
                                            'Cap.c04.PAPSS2', 'Cap.c05.RGCC', 'Cap.c07.CA2', 'Cap.c08.TMEM163'))
seurat_obj$scaled_RTS <- scaled_rts_values
low_colors <- colorRampPalette(colors = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0"))(50)
high_colors <- colorRampPalette(colors = c("#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(50)
custom_colors <- c(low_colors, high_colors)
p2 <- FeaturePlot(seurat_obj, features = "scaled_RTS",pt.size = 0.2) +
  labs(title = 'R T/S value') +
  theme(panel.grid = element_blank(), aspect.ratio = 1) +
  scale_color_gradientn(colors = custom_colors, limits = c(0.5, 1.5), oob = scales::squish)
ggsave(file.path(outdir2, 'RTS FeaturePlot.pdf'), p2, height=4.5, width=4.5)

#### 5.Fig4e ----
root_type <- "Cap.c08.TNFRSF4"
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
saveRDS(result_umap, './files/DCmap_umap_celltype.rds')

result <- readRDS('./files/DCmap_umap_celltype.rds')
TYPE <- plot_ly(result$loc, x = ~DC1, y = ~DC2, z = ~DC3, color = ~cell_type,
                colors = col_names,
                type = 'scatter3d', mode = 'markers') %>%
  layout(title = "3D Diffusion Map",
         scene = list(xaxis = list(title = 'DC 1'),
                      yaxis = list(title = 'DC 2'),
                      zaxis = list(title = 'DC 3')))

#### 6.FigS8d ----
library(scrat)
set.seed(016)
seurat_obj <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
seurat_obj <- subset(seurat_obj, idents = c('Art.c01.DKK2', 'Art.c02.NEBL', 'Art.c03.VEGFA', 
                                            'Cap.c04.PAPSS2', 'Cap.c05.RGCC', 'Cap.c07.CA2', 'Cap.c08.TMEM163'))
data_matrix <- seurat_obj@assays$RNA@counts
env <- scrat.new(list(dataset.name="Cap_Art_1000c"))
env$indata <- data_matrix
results <- scrat.run(env)
qsave(results, file.path(outdir, 'Cap_Art_1000c.qs'))

Names <- list.files(path = "./Cap_Art_1000c - Results/CSV Sheets/Spot Lists/overexpression spots/",
                    pattern = "\\.csv$", full.names = FALSE)
df_all <- data.frame()
for (i in 1:length(Names)) {
  df <- read.csv(Names[i], header = TRUE)
  df$group <- Names[i]
  df <- df[, c('ID', 'Rank', 'group')]
  df_all <- rbind(df_all, df)
}
df_all$group <- gsub('.csv', '', df_all$group)
write.csv(df_all, file = './files/Overexpression Spots All.csv')

#### 7.FigS8e ----
library(GSEA)
library(org.Hs.eg.db)

df_all <- read.csv('./files/Overexpression Spots All.csv', row.names = 1)
bp <- enrichGO(df_all$ID,
               OrgDb = org.Hs.eg.db,
               keyType = 'SYMBOL',
               ont = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.2)
term <- bp@result
saveRDS(term, file.path(outdir, "go_function.rds"))
write.table(term,  file = file.path(outdir, "go_function.csv"), row.names = F, quote = F, sep = ",")

term <- readRDS(file.path(outdir, "go_function.rds"))
term_ok <- term |> filter(pvalue <= 0.05)
pathway <- c(
  # Angiogenesis
  'sprouting angiogenesis',
  'regulation of angiogenesis',
  'vasculogenesis',
  # Inflammation
  'NIK/NF-kappaB signaling',
  'inflammasome complex assembly',
  'interferon-mediated signaling pathway',
  # Vasculature development
  'Wnt signaling pathway',
  'Notch signaling pathway',
  'regulation of vasculature development',
  # ECM organization
  'extracellular matrix organization',
  'extracellular matrix assembly',
  'cell-substrate adhesion',
  # Hypoxia
  'response to decreased oxygen levels',
  'response to hypoxia',
  'response to oxygen levels',
  # Platelet aggregation
  'platelet activation',
  'platelet-derived growth factor receptor signaling pathway',
  'response to platelet-derived growth factor',
  # Immune response
  'activation of immune response',
  'immune response-regulating signaling pathway',
  'positive regulation of leukocyte cell-cell adhesion',
  # Phospholipid metabolism
  'glycerophospholipid biosynthetic process',
  'glycerophospholipid metabolic process',
  'phospholipid transport',
  # TGFβ signaling
  'cellular response to transforming growth factor beta stimulus',
  'transforming growth factor beta receptor signaling pathway',
  'vascular endothelial growth factor receptor signaling pathway',
  # Cell proliferation and migration
  'endothelial cell proliferation',
  'endothelial cell migration',
  'blood vessel endothelial cell migration')

df <- metascape |> filter(Description %in% terms)
df <- df[,c("Category","CategoryID","GO","Description","LogP","Enrichment","GeneID", "Hits")]
df$Description <- factor(df$Description, levels = terms)
df <- df[order(df$Description), ]
G_names <- c('Angiogenesis','Inflammation','Vasculature development','ECM organization','Hypoxia','Platelet aggregation',
             'Immune response','Phospholipid metabolism','TGFβ signaling','Cell proliferation and migration')
df$group <- factor(rep(G_names, each = 3), levels = G_names)
df <- df %>% arrange(LogP) %>%
  mutate(labelx = 0, labely = 1:n())
saveRDS(df, file.path(outdir, 'metascape_DAPFs.rds'))

ggplot(data = df, aes(x = LogP,  fill = group,y = reorder(Description, LogP))) +
  geom_bar(stat = "identity", alpha = 1, width = 0.8) +
  geom_text(aes(x = labelx, y = labely, label = Description), size = 3.5, hjust = 0) +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12)) +
  xlab("-log10(pvalue)") + ggtitle('DAPFs') +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = c("#1B9E77", "#FB8072", "#57C3F3", "#FDB462", "#F0E442", 
                               "#FCCDE5", "#BC80BD", "#8DD3C7", "#377EB8", "#CCEBC5"))

ggsave(file.path(outdir2, 'Metascape_Barplot.pdf'), last_plot(), height=6, width=6)

#### 8.Fig4f ----
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
  # Vasculature development
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

#### 9.Fig4g ----
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

#### 10.FigS8g ----
seurat_obj0 <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
seurat_obj <- subset(seurat_obj0, idents = c('Art.c01.DKK2', 'Art.c02.NEBL', 'Art.c03.VEGFA',  
                                             'Cap.c04.PAPSS2', 'Cap.c05.RGCC',
                                             'Cap.c07.CA2', 'Cap.c08.TMEM163'))

PATHWAY_NOTCH <- c('DLL1', 'JAG1','JAG2','LFNG','NOTCH1','NOTCH2','NOTCH3', 'NOTCH4', 'NRARP',
                   'RBPJ', 'HES1', 'HES5', 'HEY1', 
                   'DVL1','DVL2','NUMB','ITCH','DTX1','DTX2','ADAM17','NCSTN','APH1A','APH1B','MAML1','MAML2','MAML3', 
                   'EP300', 'KAT2A','KAT2B', 'SNW1', 
                   'HR','CTBP1','CTBP2', 'TLE1','TLE2','TLE3','TLE4', 'NCOR2', 'CIR1',
                   'HDAC1','HDAC2','ATXN1') 
## UMAP
seurat_obj <- AddModuleScore(seurat_obj, features = list(PATHWAY_NOTCH), name = 'NOTCH',seed = 111)
data <- data.frame(seurat_obj@meta.data, seurat_obj@reductions$umap@cell.embeddings)
data$NOTCH1[data$NOTCH1 > 0.3] <- 0.3
ggplot(data,  aes(UMAP_1, UMAP_2, color=NOTCH1)) + 
  geom_point(size=0.5) + labs(title = "")+
  scale_color_distiller(palette= "RdBu")+
  theme_dr(xlength = 0.2, ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"), ends = 'last', type = "closed")) + 
  guides(color=guide_legend(title = "NOTCH score",override.aes = list(size = 5),reverse = T))+
  theme(aspect.ratio = 1, panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5,vjust = -1))
ggsave(file.path(outdir2, 'score_NOTCH_celltype_UMAP.pdf'),
       ggplot2::last_plot(),width = 4, height = 4)

#### 11.FigS8h ----
ggviolin(seurat_obj@meta.data, x='cell_type',y="NOTCH1",
         xlab = '', ylab = 'NOTCH pathway score',
         fill='cell_type',palette = type_colors,add = 'boxplot',
         add.params = list(fill='white', width=0.3),order = C_type)+
  theme(axis.text.x = element_blank())+NoLegend()
ggsave(file.path(outdir2, "score_NOTCH_celltype.pdf"),
       ggplot2::last_plot(), height=3, width=5)

#### 12.FigS8i ----
sub_obj <- subset(seurat_obj0, idents = c("Art.c01.DKK2","Art.c02.NEBL","Art.c03.VEGFA"))
s_genes <- c('NOTCH1', 'JAG1', 'JAG2', 'DLL4', 'HEY1', 
             'NCSTN','APH1A', 'KAT2A', 'SNW1', 'HES1', 'HES5') 
p2 <- DotPlot(sub_obj, group.by = 'disease', features= s_genes)+
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(angle=45,hjust=1,vjust=1,face = 'italic'),
        axis.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8))+
  scale_color_distiller(palette= "YlOrRd", direction = 1)
ggsave(file.path(outdir2, "G_Expre_celltype_disease.pdf"), p2, height=3, width=4.5)

#### 13.Fig4h ----
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

#### 14.Fig4j ----
seurat_obj <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
seurat_obj <- subset(seurat_obj, classic1 == 'ArtECs')
sub_obj <- seurat_obj[, seurat_obj[["disease"]] == 'ICM']
selected_genes <- c('EFNB2', 'NOTCH1','DLL4')
expr <- GetAssayData(sub_obj, slot = "data", assay = "RNA")[selected_genes, ] %>% as.data.frame()
expr <- as.data.frame(t(expr))

p <- expr %>%
  filter(NOTCH1 != 0, EFNB2 != 0) %>%  # 去除X或Y中含有0值的行
  ggplot(aes(x = NOTCH1, y = EFNB2)) +
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
ggsave(file.path(outdir2, paste0("Cor_NOTCH1_Artc_EFNB2_ICM.pdf")),
       p, height=3.5, width=3.5)

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

#### 15.FigS8j ----
#### 15.1 correlation ----
selected_genes <- c('GJA4', 'GJA5', 'EFNB2', 'NOTCH1','DLL4','HEY1', 'HES1', 'HES5', 'NCSTN', 'NRARP', 'NR2F2', 'UNC5B')
for(S_D in names(table(seurat_obj$disease))){
  sub_obj <- seurat_obj[, seurat_obj[["disease"]] == S_D]
  expr <- GetAssayData(sub_obj, layer = "data", assay = "RNA")[selected_genes, ] %>% as.data.frame()
  expr <- as.data.frame(t(expr))
  
  GENES <- c('NOTCH1','DLL4', 'NCSTN', 'NRARP', 'NR2F2', 'HES5', 'HES1', 'HEY1')
  plots <- list()
  cor_results_df <- data.frame(Gene = character(),
                               R_value = numeric(),
                               P_value = numeric(),
                               stringsAsFactors = FALSE)
  for (gene in GENES) {
    if (!(gene %in% colnames(expr))) {
      warning(paste("Gene", gene, "not found in expression data for disease:", S_D))
      next
    }
    filtered_data <- expr %>%
      filter(!!sym(gene) != 0, EFNB2 != 0)
    cor_test <- cor.test(filtered_data[[gene]], filtered_data$EFNB2)
    cor_results_df <- rbind(cor_results_df, data.frame(Gene = gene,
                                                       R_value = cor_test$estimate,
                                                       P_value = cor_test$p.value))
  }
  saveRDS(cor_results_df, file.path(outdir, paste0(S_D, "_cor.rds")))
}

#### 15.2 plot ----
Names <- list.files(path = outdir, pattern = "\\.rds$", full.names = FALSE)
df_all <- data.frame()
for (i in 1:length(Names)) {
  df <- readRDS(file.path(outdir, Names[i]))
  df$group <- Names[i]
  df$group <- gsub('_cor.rds', '', df$group)
  df_all <- rbind(df_all, df)
}
saveRDS(df_all, file.path(outdir, "Disease_cor.rds"))

cor_df <- df_all |> filter(Gene %in% c('NOTCH1', 'DLL4'))
cor_df$R_value[cor_df$R_value>0.6] = 0.6
cor_df$P_value[cor_df$P_value<1.5e-20] = 1.5e-20
cor_df <- cor_df %>%
  mutate(size = -log10(P_value), R_value = as.numeric(R_value)) 
cor_df$group <- factor(cor_df$group, levels = c('Normal', 'ACM', 'DCM', 'HCM', 'ICM', 'AMI', 'CHD', 'CS', 'COV' ))

ggplot(cor_df, aes(x = group, y = Gene)) +
  geom_point(aes(size = size, color = R_value)) +
  scale_size_continuous(range = c(2, 8)) + 
  scale_color_gradient2(low = "#fffbd5", mid = "#ee9ca7", high = '#B2182B', midpoint = 0.15) + 
  labs(x = NULL, y = NULL, size = "-log10(P-value)",
       color = "Correlation") +
  theme_classic2() +
  theme(axis.text.y = element_text(face = 'italic'),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
        axis.text = element_text(colour = 'black'))
ggsave(file.path(outdir2, "Cor_Dotplot_notch_art.pdf"),
       last_plot(), height=4, width=5.5)
