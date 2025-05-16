library(Seurat)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(qs)
library(tidydr)
library(dplyr)
library(UpSetR)
library(ggpubr)
library(ggbreak)
library(reshape2)
source('~/Collection/code/scrna-seq.R')
source('~/Collection/code/plot.R')

path <- '~/PH/results/Fig1/'
setwd(path)
outdir <- './files/'
outdir2 <- './plots/'
case <- 'Fig1_'

#### 1.Fig1b ----
seurat_obj <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
disease_colors <- c("#B3DE69","#FB8072", "#80B1D3", "#FDB462", "#8DD3C7",
                    "#FCCDE5","#CCEBC5", "#377EB8", "#1B9E77")
data_n <- as.data.frame(t(table(seurat_obj$disease)))
data_n <- data_n[,-1]
colnames(data_n) <- c('group', 'disease')
result <- seurat_obj@meta.data %>%
  group_by(disease) %>%
  summarize(unique_sample_count = n_distinct(sample))
data_n$sample <- result$unique_sample_count

for (i in types) {
  case <- i
  sub_obj <- subset(seurat_obj, disease == case)
  cell_name <- c()
  for (j in 1:nrow(sub_obj@meta.data)) {
    a <- strsplit(sub_obj@meta.data[j,5], ',')[[1]]
    cell_name <- c(cell_name, a)
  }
  cell_name <- unique(cell_name)
  print(length(cell_name))
}
data_n$disease <- c('48341','13208','116168','46785','20581','16802','25332','3826','11862')
data_n$disease <- as.numeric(data_n$disease)

p1 <- ggplot(data_n, aes(x = group, y = disease, fill = group)) +
  geom_bar(stat = "identity", color = "white", position = "dodge") +
  scale_fill_manual(values = disease_colors) +
  labs(title = "", x = "", y = "Number of cells") +
  theme_classic()+NoLegend()+
  geom_text(aes(label=disease),size=3,
            position = position_dodge(width = 0.8),
            vjust=0.7)+
  scale_y_break(breaks = c(50000, 110000), scales = 0.2, ticklabels = c(110000,115000))+
  theme(plot.margin = margin(100,5.5,5.5,5.5),
        axis.text =  element_text(color = 'black'),
        panel.grid = element_blank())
p1
p2 <- ggplot(data_n, aes(x = group, y = sample, fill = group)) +
  geom_bar(stat = "identity", color = "white", position = "dodge") +
  scale_fill_manual(values = disease_colors) +
  labs(title = "", x = "", y = "Number of samples") +
  theme_classic()+NoLegend()+
  geom_text(aes(label=sample),size=3,
            position = position_dodge(width = 0.8),
            vjust=0.7)+
  scale_y_break(breaks = c(60, 120), scales = 0.2, ticklabels = c(120,140))+
  theme(plot.margin = margin(10,5.5,5.5,5.5),
        axis.text =  element_text(color = 'black'),
        panel.grid = element_blank()) 
p1/p2 + plot_layout(ncol = 1, heights = c(0.5, 0.5))

ggsave(file.path(outdir2,paste0(case,'BarPlot_numbers.pdf')),
       ggplot2::last_plot(),
       height=6,
       width=5)

#### 2.Fig1c ----
type_colors <- c("#f2ccac","#edae11","#FE870D","#e1abbc","#F88A89","#c04932",
                 "#c1d5b9","#99CD91","#009E73","#6a73cf","#916CB2","#237BB2","#8A9D5B")

## UMAP
df = data.frame(clu=names(table(seurat_obj$cell_type)),
                per=sprintf("%1.2f%%", 100*table(seurat_obj$cell_type)/length(seurat_obj$cell_type)))
seurat_obj$type_per = df[match(seurat_obj$cell_type,df$clu),2]
seurat_obj$type_per = paste0(seurat_obj$cell_type," (",seurat_obj$type_per,")")
seurat_obj$type_per <- factor(seurat_obj$type_per, 
                              levels = c('Art.c01.DKK2 (3.51%)','Art.c02.NEBL (3.44%)','Art.c03.VEGFA (8.01%)','Cap.c04.PAPSS2 (12.96%)', 
                                         'Cap.c05.RGCC (21.46%)','Cap.c06.RDH10 (5.40%)','Cap.c07.CA2 (25.30%)','Cap.c08.TMEM163 (5.52%)','Cap.c09.RGS6 (7.08%)', 
                                         'Ven.c10.VCAM1 (3.87%)','Ven.c11.ACTA2 (1.90%)','LEC.c12.FLT4 (0.55%)','Endo.c13.NPR3 (0.99%)'))

DimPlot(seurat_obj, cols = type_colors, group.by = 'type_per')+
  theme_void() +
  theme_dr(xlength = 0.2, ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"),
                               ends = 'last', type = "closed")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs( x= "UMAP 1",y= "UMAP 2")+
  theme(plot.margin = margin(5.5,5.5,5.5,5.5),
        panel.grid = element_blank(),aspect.ratio = 1) 
ggsave(file.path(outdir2,paste0(case,'DimPlot_type_per.pdf')),
       ggplot2::last_plot(), height=6, width=6)

## dotplot
Name <- rev(names(table(seurat_obj$cell_type)))
seurat_obj$cell_type <- factor(seurat_obj$cell_type, levels = Name)
DotPlot(seurat_obj, group.by = 'cell_type',
        features= rev(c('DKK2', 'NEBL', 'VEGFA', 
                        'PAPSS2', 'RGCC','RDH10', 'CA2', 'TMEM163', 'RGS6', 'VCAM1', 'ACTA2', 
                        'FLT4', 'NPR3')))+
  labs(x= 'logFC' , y="-log10(p_val)")+ theme_minimal()+
  theme( axis.text = element_text(size = 8),
         axis.text.y = element_blank(),
         axis.text.x = element_text(angle=45,hjust=1,colour = 'black',vjust=1,face = 'italic'),
         axis.title = element_blank(),
         legend.text = element_text(size = 8),
         legend.title = element_text(size = 8),
         axis.ticks = element_blank())+
  scale_color_distiller(palette= "RdBu")
ggsave(file.path(outdir2,paste0(case,'DotPlot_type_marker.pdf')),
       ggplot2::last_plot(), height=3.5, width=4)

#### 3.Fig1d ----
cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@data) 
## 1.database
Dataset <- msigdbr(species = "Homo sapiens")
Dataset$gs_name <- gsub('GOBP_', '', Dataset$gs_name)
Dataset$gs_name <- gsub('GOMF_', '', Dataset$gs_name)
Dataset$gs_name <- gsub('GOCC_', '', Dataset$gs_name)
length(unique(Dataset$gs_name))
#挑选通路
vascular_GS <- read_csv('~/Collection/Data_P_H/files/615 vascular related gene sets.csv')
vascular_GS$Geneset <- gsub('GO_', '', vascular_GS$Geneset)
other <- c('CIRCADIAN_RHYTHM',
           'RHYTHMIC_BEHAVIOR',
           'RHYTHMIC_PROCESS',
           'CIRCADIAN_SLEEP_WAKE_CYCLE',
           'ENTRAINMENT_OF_CIRCADIAN_CLOCK',
           'REGULATION_OF_CIRCADIAN_RHYTHM',
           'KEGG_CIRCADIAN_RHYTHM_MAMMAL',
           'REGULATION_OF_CIRCADIAN_SLEEP_WAKE_CYCLE',
           'PID_CIRCADIAN_PATHWAY',
           'REACTOME_CIRCADIAN_CLOCK',
           'WP_CIRCADIAN_RHYTHM_GENES')
feature_list <- vascular_GS$Geneset
dataset <- Dataset %>% filter(Dataset$gs_name %in% c(feature_list, other))

## 2.AUCell
geneSets <- lapply(unique(dataset$gs_name), 
                   function(x){dataset$gene_symbol[dataset$gs_name == x]})
names(geneSets) <- unique(dataset$gs_name)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
colnames(AUCell_socre) <- colnames(AUCell_socre)
saveRDS(AUCell_socre,file.path(outdir, paste0(case,"AUCscore.rds")))

## 3.cell_type
# AUCell_socre <- readRDS(file.path(outdir, paste0(case,"AUCscore.rds")))
seurat_obj$ID <- rownames(seurat_obj@meta.data)
meta <- AUCell_socre
meta <- as.data.frame(meta)
#colnames(meta)
meta$cell_type <- seurat_obj$cell_type
object <-
  aggregate(meta[, -length(colnames(meta))], list(meta$cell_type), FUN = mean)
rownames(object) <- object$Group.1
object <-object[,-1]
write.table(object, file.path(outdir,paste0(case,'AUCscore_celltypes.csv')), sep=",",quote=F)
saveRDS(object,file.path(outdir, paste0(case,"AUCscore_celltypes.rds")))

## 4.heatmap
# object <- readRDS(file.path(outdir, paste0(case,"AUCscore_celltypes.rds")))
colnames(object) <- colnames(object)|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()

## 挑选通路绘图
object <- readRDS(file.path(outdir, paste0(case,"AUCscore_celltypes.rds")))
colnames(object) <- colnames(object)|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()
a <- as.data.frame(colnames(object))
pathway <- c(
  'Artery development',
  'Vasculogenesis',
  'Regulation of vasoconstriction',
  'Kim hypoxia',
  'Regulation of vasculature development',
  'Cellular lipid metabolic process',
  'Neutral lipid catabolic process',
  'Phospholipid metabolic process',
  'Kegg glycerolipid metabolism',
  'Membrane lipid metabolic process',
  'Triglyceride catabolic process',
  'Cell matrix adhesion',
  'Cell cell junction',
  'Biocarta il1r pathway',
  'Biocarta nfkb pathway',
  'Fulcher inflammatory response lectin vs lps up',
  'Karlsson tgfb1 targets up',
  'Reactome mhc class ii antigen presentation',
  'Petrova endothelium lymphatic vs blood up',
  'Petrova prox1 targets up',
  'Regulation of cell adhesion mediated by integrin',
  'Kegg circadian rhythm mammal',
  'Hormone metabolic process'
)

pathway <- pathway[!duplicated(pathway)]
pathway[!(pathway %in% colnames(object))]

object2 <- object[ ,colnames(object) %in% pathway]
object2 <- as.data.frame(t(object2))

object2$pathway <- rownames(object2)
object2$pathway <- factor(object2$pathway, levels = pathway)
object2 <- object2[order(object2$pathway), ]
saveRDS(object2,file.path(outdir, paste0(case,"AUCscore_celltypes_s.rds")))

p1 <- pheatmap(object2, treeheight_row = 5,
               cluster_cols = F, cluster_rows = F, scale="row",
               fontsize_number = 20, border="white",
               shape = "circle",cellwidth = 10,cellheight =10,
               gaps_row = c(4,12,18,21),
               col = c(colorRampPalette(colors = c("#2166AC","#4393C3","#92C5DE","#D1E5F0"))(50),
                       colorRampPalette(colors = c("#FDDBC7","#F4A582","#D6604D","#B2182B"))(50)),
               fontsize_row = 10,fontsize_col = 12) 
ggsave(file.path(outdir2,paste0(case,'AUCscpre_sub_row2.pdf')),
       p1,limitsize = FALSE, height=5, width=6)

#### 4.Fig1e ----
object <- readRDS(file.path(outdir, paste0(case,"AUCscore.rds")))
pathway <- c( "VASCULOGENESIS", 
              "VASCULATURE_DEVELOPMENT",
              'KIM_HYPOXIA',
              'KEGG_GLYCEROLIPID_METABOLISM',
              "BIOCARTA_NFKB_PATHWAY",
              "PETROVA_ENDOTHELIUM_LYMPHATIC_VS_BLOOD_UP", 
              "HORMONE_METABOLIC_PROCESS")
data <- object[, pathway]
# 绘图
for(i in 1:length(pathway)){
  seurat_obj$group <- seurat_obj$cell_type
  seurat_obj$AUC  <- data[,i]
  
  meta.data <- seurat_obj@meta.data[, c("group", "AUC")]
  head(meta.data)
  meta.data <- meta.data %>% reshape2::melt(id.vars = c("group"))
  head(meta.data)
  colnames(meta.data)[2:3] <- c("signature", "score")
  mean.score <- median(meta.data$score)
  median.score <- meta.data %>% group_by(group) %>%
    summarise(median_score = median(score))
  table(median.score$group)
  median.score <- median.score[order(median.score$median_score,decreasing = T),]
  p <- meta.data %>%
    ggplot(aes(x = group, y = score, fill = '',
               color =group)) +
    geom_violin(fill = 'NA') +
    geom_boxplot(data = meta.data, fill = 'NA',
                 aes(x = group, y = score, color =group),
                 width = 0.1,  alpha = 1)+
    geom_hline(yintercept = mean.score,
               color = "darkgrey", linetype = "dashed") +
    scale_y_continuous(pathway[i] |>
                         str_replace_all("_", " ")|>
                         str_to_sentence()) +
    theme_classic() + scale_color_manual(values = type_colors) +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1),
          axis.text = element_text(color = 'black'))
  p
  ggsave(file.path(outdir2,paste0(pathway[i],'_violin.pdf')),
         p, height = 2.5,width = 8)
}

#### 5.Fig1f ----
library("sscVis")
library("data.table")
library("grid")
library("cowplot")
library("ggrepel")
library("readr")
library("plyr")
library("ggpubr")
library("ggplot2")
library("dplyr")
library("tidyr")
library('pheatmap')
library('circlize')
library('dendextend')
cellInfo.tb <- seurat_obj@meta.data
A <- do.tissueDist(cellInfo.tb = cellInfo.tb,
                   meta.cluster = cellInfo.tb$cell_type,
                   colname.patient = "donor",
                   loc = cellInfo.tb$disease,
                   out.prefix = 'EC',
                   pdf.width=8,
                   pdf.height=5,
                   verbose=1)
saveRDS(A, file.path(outdir, 'OR_df.rds'))

#
A$OR.dist.mtx
A$p.dist.tb
A$OR.dist.tb
A$count.dist.melt.ext.tb

data <- as.data.frame(A$OR.dist.tb)
rownames(data) <- data$rid
data <- data[,-1]
data[data > 2] <- 2
data <- data[c('Art.c01.DKK2', 'Art.c02.NEBL', 'Art.c03.VEGFA', 'Cap.c04.PAPSS2', 'Cap.c05.RGCC', 
               'Cap.c06.RDH10', 'Cap.c07.CA2', 'Cap.c08.TMEM163', 'Cap.c09.RGS6', 'Ven.c10.VCAM1', 
               'Ven.c11.ACTA2', 'LEC.c12.FLT4' , 'Endo.c13.NPR3'), ]

##
data_p <- as.data.frame(A$p.dist.tb)
rownames(data_p) <- data_p$rid
data_p <- data_p[,-1]
data_p <- data_p[c('Art.c01.DKK2', 'Art.c02.NEBL', 'Art.c03.VEGFA', 'Cap.c04.PAPSS2', 'Cap.c05.RGCC', 
                   'Cap.c06.RDH10', 'Cap.c07.CA2', 'Cap.c08.TMEM163', 'Cap.c09.RGS6', 'Ven.c10.VCAM1', 
                   'Ven.c11.ACTA2', 'LEC.c12.FLT4' , 'Endo.c13.NPR3'), ]
data_p[data_p <= 0.01] <- "*"
data_p[data_p > 0.01] <- ""

##
data2 <- data
data2[data2 >= 1] <- "*"
data2[data2 !=  "*"] <- ""
p_text <- ifelse(data2 == "*" & data_p == "*", "*", "")
colnames(p_text) <- colnames(data2)
rownames(p_text) <- rownames(data2)

## plot
annotation_row = data.frame(Group = factor(rep(c("ArtECs","CapECs",'VenECs',
                                                 'LECs','EndoECs'), c(3, 6, 2, 1, 1))))
rownames(annotation_row) = c('Art.c01.DKK2', 'Art.c02.NEBL', 'Art.c03.VEGFA', 'Cap.c04.PAPSS2', 'Cap.c05.RGCC', 
                             'Cap.c06.RDH10', 'Cap.c07.CA2', 'Cap.c08.TMEM163', 'Cap.c09.RGS6', 'Ven.c10.VCAM1', 
                             'Ven.c11.ACTA2', 'LEC.c12.FLT4' , 'Endo.c13.NPR3')
ann_colors = list(Group = c("ArtECs"="#edd064","CapECs"="#f2ccac",'VenECs'="#a1d5b9",
                            'LECs'="#57C3F3",'EndoECs'='#6778AE'))

pdf(file.path(outdir2, 'OR_heatmap.pdf'),width = 5,height = 5)
pheatmap(data, scale = "none", cluster_rows=F, cluster_cols = F,
         display_numbers = p_text, fontsize_number = 10, number_color = 'white',
         color = colorRampPalette(c("white","#FFE4B5","#D55E00","darkred"))(256),
         cellwidth = 10, cellheight = 10, border_color = "white", fontsize = 8,
         annotation_colors = ann_colors,
         annotation_row = annotation_row, 
         treeheight_row = 10, treeheight_col =10,
         gaps_col=c(1),gaps_row = c(3,9,11,12),
         main = "Odds Ration",show_rownames = T,show_colnames = T)
dev.off()

#### 6.Fig1g ----
object <- readRDS(file.path(outdir, paste0(case,"AUCscore_celltypes.rds")))
colnames(object) <- colnames(object)|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()
object3 <- object[ ,colnames(object) %in% c('Endothelial cell migration', 'Endothelial cell proliferation')]
object3$type <- rownames(object3)
object3$type <- factor(object3$type, levels = object3$type)
ggplot(object3, aes(x = `Endothelial cell migration`, y = `Endothelial cell proliferation`, 
                    label = rownames(object3), color = type)) + 
  geom_point(size = 5) +
  geom_text(hjust = 1.1, vjust = 0.5, size = 3) + 
  theme_classic() + 
  labs(title = " ", x = "Migration Score", y = "Proliferation Score") +
  theme(axis.text = element_text(colour = 'black'), aspect.ratio = 1) +
  NoLegend() +
  scale_color_manual(values = type_colors) +
  geom_hline(yintercept = mean(object3[,2]), linetype = "dashed", size = 0.5, color = "grey") +
  geom_vline(xintercept = mean(object3[,1]), linetype = "dashed", size = 0.5, color = "grey") +
  scale_x_continuous(limits = c(0.155, 0.20)) +
  scale_y_continuous(limits = c(0.115, 0.155))
ggsave(file.path(outdir2,paste0(case,'AUCscpre_pointplot.pdf')),
       ggplot2::last_plot(), height=4, width=4)
