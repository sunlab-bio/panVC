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
library(ggpubr)
source('~/Collection/code/scrna-seq.R')
source('~/Collection/code/plot.R')

path <- '~/PH/results/Fig6/'
setwd(path)
outdir <- './files/'
outdir2 <- './plots/'

#### 1.Fig6a ----
seurat_obj0 <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
seurat_obj <- subset(seurat_obj0, cell_type == 'Endo.c13.NPR3')

seurat_obj$log2_NR1D1_ARNTL <- log2(seurat_obj[["RNA"]]@data["NR1D1", ] / seurat_obj[["RNA"]]@data["ARNTL", ])
seurat_obj$log2_NR1D1_ARNTL <- seurat_obj$log2_NR1D1_ARNTL

#### group
group <- seurat_obj@meta.data$group
plot_data <- data.frame(log2_NR1D1_ARNTL = seurat_obj$log2_NR1D1_ARNTL, group = group)

ggplot(plot_data, aes(x = group, y = log2_NR1D1_ARNTL, fill = group)) +
  geom_boxplot() +
  theme_classic2() +
  labs(y = "Log2(NR1D1/ARNTL)") +
  scale_fill_manual(values = c("#B3DE69","#b55f31")) +
  theme(axis.text = element_text(color = 'black'), axis.title.x = element_blank())+
  stat_compare_means(comparisons = list(c('Health', 'CVD')), vjust = 2)+NoLegend()
ggsave(file.path(outdir2, "Box_Log2.pdf"), last_plot(), height=3, width=2.5)

#### disease
disease_colors <- c("#B3DE69","#FB8072", "#80B1D3", "#FDB462", "#8DD3C7", "#FCCDE5","#CCEBC5", "#377EB8", "#1B9E77")

sub_obj <- subset(seurat_obj, disease %in% c('Normal','ACM','DCM','HCM','ICM','AMI','CHD'))
group <- sub_obj$disease
plot_data <- data.frame(log2_NR1D1_ARNTL = sub_obj$log2_NR1D1_ARNTL, group = group)
plot_data <- plot_data %>%
  filter(is.finite(log2_NR1D1_ARNTL) & !is.na(log2_NR1D1_ARNTL) & !is.na(group))
plot_data$group <- factor(plot_data$group)

ggplot(plot_data, aes(x = group, y = log2_NR1D1_ARNTL, fill = group)) +
  geom_boxplot() +
  theme_classic2() +
  labs(y = "Log2(NR1D1/ARNTL)") +
  scale_fill_manual(values = disease_colors) +
  theme(axis.text = element_text(color = 'black'), axis.title.x = element_blank())+
  stat_compare_means(comparisons = list(c('Normal', 'DCM'), c('Normal', 'CHD')), label = 'p.signif') + NoLegend()
ggsave(file.path(outdir2, "Box_Log2_disease.pdf"), last_plot(), height=5, width=5)

#### 2.Fig6b ----
library(ggraph)
CRGs <- read.table('./files/CRGs_35436363.txt') |> pull()
CRGs <- CRGs[(CRGs %in% rownames(seurat_obj))]

cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@counts) 
cells_AUC <- AUCell_calcAUC(CRGs, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
saveRDS(AUCell_socre, file.path(outdir, 'CRD score celltype.rds'))

#### 2.1 UMAP ----
seurat_obj$AUC  <- AUCell_socre
seurat_obj$AUC[seurat_obj$AUC> 0.25] <- 0.25
seurat_obj$AUC[seurat_obj$AUC< 0.20] <- 0.20
data <- data.frame(seurat_obj@meta.data, seurat_obj@reductions$umap@cell.embeddings)
ggplot(data,  aes(UMAP_1, UMAP_2, color=AUC)) + 
  geom_point(size=0.1) +
  scale_color_distiller(palette= "YlOrRd", direction = 1)+ 
  theme_dr(xlength = 0.2, ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"), ends = 'last', type = "closed")) +
  theme(aspect.ratio = 1, panel.grid = element_blank(),
        plot.title = element_text(hjust=0.5,vjust = -1))
ggsave(file.path(outdir2, 'CRD score UMAP.pdf'), last_plot(),width = 4, height = 4)

#### 2.2 Vlnplot ----
ggviolin(data, x='cell_type',y="AUC",
         xlab = '', ylab = "CRDGs score",
         fill='cell_type', add = 'boxplot',
         add.params = list(fill='white', width=0.15))+ 
  scale_fill_manual(values = type_colors)+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) + NoLegend()
ggsave(file.path(outdir2, 'CRD score vlnplot.pdf'),
       ggplot2::last_plot(),width = 6, height = 3.5)

#### 3.Fig6c ----
#### 3.1 cell_type ----
features <- c('PER1','PER2','CRY1','CRY2','CLOCK','CTNND2','SEMA5A','TIMELESS','BTRC', 'DDR2')
DotPlot(seurat_obj, features= features, group.by = 'cell_type')+
  theme(axis.text = element_text(size = 8),
        axis.text.x = element_text(angle=30, vjust=1,hjust=1),
        axis.text.y = element_text(angle=0,face = 'italic'),
        axis.title = element_blank(),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        axis.line = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        axis.ticks = element_blank())+coord_flip()+
  scale_color_distiller(palette= "RdBu")
ggsave(file.path(outdir2, 'CRD genes types dotplot.pdf'), last_plot(), width = 6, height = 2.5)

#### 3.2 disease ----
features <- c('PER1','PER2','CRY2','ARNTL', 'DIXDC1', 'ARHGAP26', 'CELF2',
              'SEMA5A','BTRC','CRY1', 'SLC12A2','SLC12A7','SLC12A4', 'BHLHE40')
expr <- AverageExpression(seurat_obj, assays = "RNA", group.by = "disease")[["RNA"]] %>%
  as.data.frame()
pdf(paste0(outdir2, 'CRD genes disease heatmap.pdf'), width = 3, height = 2.5)
pheatmap(expr[features,], scale = "row",
         color = c(colorRampPalette(colors = c("#2166AC","#4393C3","#92C5DE","#D1E5F0"))(50),
                   colorRampPalette(colors = c("#FDDBC7","#F4A582","#D6604D","#B2182B"))(50)), 
         cellwidth = 8, cellheight = 8, fontsize = 8, border_color = 'white',
         show_rownames = T,show_colnames = T, cluster_col = F, cluster_row = F,
         treeheight_row =2, clustering_method = "average")
dev.off()

#### 4.Fig6e ----
#### 4.1 DEG ----
diseases <- c("ACM","DCM","HCM","ICM","AMI","CHD","CS","COV")
DEG <- data.frame()
for(i in diseases){
  outdir3 <- paste0('./files/deg/', i, '/')
  logFC.cutoff <- 0.25
  pvalue.cutoff <- 0.05
  
  case <- i
  obj <- seurat_obj[ , seurat_obj$disease %in% c(i, "Normal")]
  Idents(obj) <- 'disease'
  markers <- FindMarkers(obj,
                         ident.1 = i,
                         ident.2 = "Normal",
                         min.pct = 0,
                         logfc.threshold = 0)
  markers$gene <- rownames(markers)
  markers$change <-
    as.factor(ifelse(
      markers$p_val_adj < pvalue.cutoff &
        abs(markers$avg_log2FC) > logFC.cutoff,
      ifelse(markers$avg_log2FC > logFC.cutoff , 'UP', 'DOWN'),
      'NOT'
    ))
  saveRDS(markers, file.path(outdir3, paste0(i,"_deg.rds")))
  markers$cluster <- case
  DEG <- rbind(DEG,markers)
}
deg_all_union <- DEG |> filter(DEG$change != 'NOT') |> pull(gene) |> unique()
write.table(deg_all_union, './files/deg_all_union.txt', sep = "\t", row.names = FALSE, col.names = FALSE)

#### 4.2 DRG ∩ CRGs ----
DEGs <- deg_all_union
DEGs_CRGs <- intersect(DEGs,CRGs)
write.table(DEGs_CRGs, file = file.path(outdir, 'DEGs_CRGs.txt'), row.names = F, col.names = F, quote = F, sep = "\t")

library(clusterProfiler)
library(org.Hs.eg.db) 

GO <- enrichGO(
  DEGs_CRGs,
  OrgDb = org.Hs.eg.db,
  keyType = 'SYMBOL', ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2)
term <- GO@result

df <- term[1:10,]
df$pvalue <- -log10(df$pvalue)
df$labelx=rep(0,nrow(df))
df$labely=seq(nrow(df),1)
ggplot(data = df, aes(Count, reorder(Description,Count), fill = pvalue)) +
  geom_bar(stat="identity", alpha=1, width = 0.8) + 
  scale_fill_gradient(low = "#6dd5fa", high = "#2980b9" ) +
  geom_text(aes(x = labelx, y = Description, label = Description), size=3.5,  hjust =0)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,size=14),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(colour = 'black',size = 10,vjust = 3),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  labs(fill = "-log10(pvalue)", x="Count", y='TOP 10 terms')+ 
  scale_y_discrete(labels=function(x) str_wrap(x, width=40))
ggsave(file.path(outdir2, 'DEGs_CRGs_GOBP.pdf'), last_plot(), height=3, width=4.5)

#### 5.Fig6f ----
seurat_obj0 <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
seurat_obj <- subset(seurat_obj0, cell_type == 'Endo.c13.NPR3')

cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@counts) 
cells_AUC <- AUCell_calcAUC(DEGs_CRGs, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
seurat_obj$AUC  <- AUCell_socre
seurat_obj$AUC[seurat_obj$AUC> 0.30] <- 0.30
seurat_obj$AUC[seurat_obj$AUC< 0.2] <- 0.2
data <- data.frame(seurat_obj@meta.data, seurat_obj@reductions$umap@cell.embeddings)
N_score <- mean(data |> filter(disease == 'Normal') |> pull(AUC))
data2 <- data |> filter(disease %in% c('Normal', 'ACM', 'DCM', 'ICM', 'CHD', 'CS', 'COV'))
ggviolin(data2, x='disease', y="AUC",
         xlab = '', ylab = "CRDGs score",
         fill='disease', add = 'boxplot',
         add.params = list(fill='white', width=0.1))+ 
  scale_fill_manual(values = disease_colors)+
  geom_hline(yintercept = N_score,
             linetype = "dashed", size = 0.5, color = "black")+
  NoLegend()+
  stat_compare_means(comparisons = list(c('Normal', 'ACM'),
                                        c('Normal', 'DCM'),
                                        c('Normal', 'ICM'),
                                        c('Normal', 'CHD'),
                                        c('Normal', 'CS'),
                                        c('Normal', 'COV')), label = 'p.signif')
ggsave(file.path(outdir2, 'DEGs_CRGs vlnplot Psignif.pdf'), last_plot(),width = 5, height = 3)

#### 6.Fig6g ----
seurat_obj0 <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
seurat_obj <- subset(seurat_obj0, cell_type == 'Endo.c13.NPR3')
Idents(seurat_obj) <- 'disease'

cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@counts) 
cells_AUC <- AUCell_calcAUC(DEGs_CRGs, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
seurat_obj$AUC  <- AUCell_socre

colnames(seurat_obj@meta.data)[length(colnames(seurat_obj@meta.data))] <- 'DEGs_CRGs'
mean(seurat_obj$DEGs_CRGs)
seurat_obj$CVD <- if_else(seurat_obj$DEGs_CRGs > mean(seurat_obj$DEGs_CRGs),
                          "CRDhigh", "CRDlow")
saveRDS(seurat_obj, file.path(outdir, 'seurat_obj_CVD.rds'))

## 
CRD_colors <- c( "#D55E00", "#194A55")
DimPlot(seurat_obj, cols = CRD_colors, group.by = 'CVD', split.by = 'group')+
  theme_void() +
  theme_dr(xlength = 0.2, ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"), ends = 'last', type = "closed")) + 
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs(title = '', x= "UMAP 1",y= "UMAP 2")+
  theme(plot.margin = margin(5.5,5.5,5.5,5.5),legend.position = c(0.5,0.8),
        panel.grid = element_blank(),aspect.ratio = 1) 
ggsave(file.path(outdir2, 'DimPlot_group.pdf'), last_plot(), height=3, width=5)

#### 7.FigS11d ----
CVD_group <- readRDS('./files/seurat_obj_CVD.rds')
colnames(CVD_group@meta.data)
seurat_obj$DEGs_CRGs <- CVD_group$DEGs_CRGs
seurat_obj$CVD <- CVD_group$CVD
het_data <- aggregate(seurat_obj$DEGs_CRGs, by=list(type=seurat_obj$sample),mean)
colnames(het_data) <- c('sample', 'score')

colnames(seurat_obj@meta.data)
info <- seurat_obj@meta.data[, c("sample", "disease", 'age_group', 'sex', 'tissue')]
info <- info[!duplicated(info$sample), ]

het_data <- merge(het_data, info, by = "sample")
het_data <- het_data[order(het_data$score),]
het_data$sample <- factor(het_data$sample, levels = het_data$sample)
het_data$score2 <- (het_data$score - min(het_data$score)) / (max(het_data$score) - min(het_data$score))
saveRDS(het_data, file.path(outdir, 'heat_data.rds'))

## barplot
p1 <- ggplot(het_data, aes(x=sample,y=score2,fill=sample))+
  geom_bar(stat = "identity",width = 0.5)+
  theme_classic()+scale_fill_manual(values = rep('#CC79A7', 171))+
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text = element_text(color = 'black'),axis.text.x = element_blank(),
        axis.title.x = element_blank(),legend.position = "none")+
  labs(title = '',y='CRDG score') 
ggsave(paste0(outdir2, 'Barplot_CRDGscore.pdf'),last_plot(), height=3, width=7)

## heatmap
annotation_col <-
  as.data.frame(het_data[, c("disease", 'age_group', 'sex', 'tissue')])
rownames(annotation_col) <- het_data$sample

annotation_colors <- list(
  disease=c('Normal' = "#B3DE69", 'ACM' = "#FB8072", 'DCM'="#80B1D3", 'HCM'="#FDB462", 'ICM'="#8DD3C7",
            'AMI'="#FCCDE5", 'CHD'="#CCEBC5", 'CS'="#377EB8", 'COV'="#1B9E77"),
  tissue=c('AX'="#E69F00", 'IVS'="#56B4E9", 'LA'="#BC80BD", 'LV'= "#009E73", 
           'RA'="#F0E442", 'RV'="#0072B2", 'RVOT'="#D55E00", 'SP'="#c74732"),
  age_group = c('<10'='#8DD3C7', '10-20'='#B3DE69', "20-30"='#FB8072', "30-40"='#1B9E77',
                "40-50"='#80B1D3', "50-60"='#BC80BD', "60-70"='#FDB462', ">70"='#377EB8'),
  sex = c('Female'="#ee6470", 'Male'="#00a6e1")) 

df <- het_data[,c('sample','score')]
rownames(df) <- df$sample
df <- df[, 2, drop = F]
df <- as.data.frame(t(df))

pdf(paste0(outdir2, "heatmap_CRDscore.pdf"), width = 8,height = 8)
pheatmap(df, scale = "none",
         annotation_col = annotation_col, annotation_colors = annotation_colors,
         color = c(colorRampPalette(colors = c("#2166AC","#4393C3","#92C5DE","#D1E5F0"))(50),
                   colorRampPalette(colors = c("#FDDBC7","#F4A582","#D6604D","#B2182B"))(50)),#红到蓝梯度
         treeheight_row = 10, treeheight_col = 10,
         border_color = "white",
         cellwidth = 2, cellheight = 5,
         fontsize = 5, 
         cluster_rows = F, cluster_cols = F,
         show_rownames = T, show_colnames = F)
dev.off()

#### 8.FigS11e ----
## UMAP
DimPlot(seurat_obj, cols = CRD_colors, group.by = 'CVD')+
  theme_void() +
  theme_dr(xlength = 0.2, ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"), ends = 'last', type = "closed"))
guides(color = guide_legend(override.aes = list(size = 5))) +
  labs(title = '', x= "UMAP 1",y= "UMAP 2")+
  theme(plot.margin = margin(5.5,5.5,5.5,5.5),
        panel.grid = element_blank(),aspect.ratio = 1) 
ggsave(file.path(outdir2, 'DimPlot_UMAP.pdf'), last_plot(), height=3, width=4)

## barplot
ggplot(data = seurat_obj@meta.data, aes(x =group, fill =seurat_obj$CVD)) +
  geom_bar(position = "fill", width = 0.6, color = "white") +
  labs(title = "", x = "", y = "Percentage of cells (%)") + 
  theme(legend.position = "top",
        legend.text = element_text(size = 8),
        legend.title = element_blank(),
        legend.background = element_blank(),
        plot.margin = margin(4, 4, 5, 4, "mm"),
        axis.title.x = element_text(size = 9),
        axis.line.x = element_line(size = 0.5),
        axis.ticks.x = element_line(size = 0.5, colour = "black"),
        axis.ticks.length = unit(.5, "lines"),
        axis.text = element_text(size = 9, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 9),
        strip.background = element_rect(color = "white", fill = "white"),
        strip.text.y = element_text(size = 9)) +
  scale_fill_manual(values = CVD_colors) + 
  scale_y_continuous(expand = c(0, 0))+coord_flip()
ggsave(file.path(outdir2, 'prpo_BarPlot_group.pdf'), last_plot(), height=3, width=5)

#### 9.Fig6h ----
#### 9.1 DEG ----
seurat_obj <- readRDS('./files/seurat_obj_CVD.rds')
Idents(seurat_obj) <- 'CVD'

top.marker <- FindMarkers(
  seurat_obj, 
  logfc.threshold = 0, min.pct = 0, 
  only.pos = F, test.use = "wilcox",
  ident.1= 'CRDhigh',ident.2= 'CRDlow' )
top.marker$gene=rownames(top.marker)
top.marker$condition=ifelse(top.marker$avg_log2FC > 0, 'CRDhigh', 'CRDlow')

top.marker=as.data.frame(top.marker)
top.marker=top.marker%>%arrange(desc(avg_log2FC))
marker_condition=data.frame()
marker_condition=marker_condition%>%rbind(top.marker)

top.marker$cluster='CVD'
marker_condition$sig=""
marker_condition$sig[abs(marker_condition$avg_log2FC) > 0.25 & marker_condition$p_val < 0.05] = "sig"
marker_condition$sig2=paste(marker_condition$cluster,marker_condition$sig,sep = "_")
marker_condition$sig2[str_detect(marker_condition$sig2,"_$")]="not_sig"
marker_condition$sig2=str_replace(marker_condition$sig2,"_sig","")

pvalue.cutoff <- 0.05
logFC.cutoff <- 0.25

marker_condition$change <-
  as.factor(ifelse(
    marker_condition$p_val < pvalue.cutoff &
      abs(marker_condition$avg_log2FC) > logFC.cutoff,
    ifelse(marker_condition$avg_log2FC > logFC.cutoff , 'UP', 'DOWN'),
    'NOT' ))
saveRDS(marker_condition, file.path(outdir, "DEG_CVD_group.rds"))

#### 9.2 volcano plot ----
library(ggrepel)
library(ggrastr)

s_genes0 <- c('NLGN1', 'PPP1CB', 'EGR1', 'PER1', 'NAMPT', 'KLF9',
              'BMP6', 'BMPR1B', 'GATA6', 'GDF7', 'GATA4', 'TGFBR3',
              'NRG1', 'SERPINE1','CDH11','SLC9C1','C22orf34','CD36', 'PELI2', 'ITGA6',
              'PLA2G5', 'PLA2R1', 'PLA2G4A', 'ANXA1', 'PTGS2', 'PTGS1', 
              'EDN1', 'FOS', 'BNIP3L', 
              'PDGFB', 'DLL4', 'KDR', 'RGCC', 'VEGFC', 'NRP1',
              'PLA2G2A', 'PITPNM2', 'PIK3C2B','PIK3R3', 'PLCB1','PLCB4')
s_genes <- unique(s_genes0)
adata <- marker_condition
adata %>% glimpse()
adata$label <- ""
adata$label[adata$gene %in% s_genes] <-
  adata[adata$gene %in% s_genes, ]$gene

p <- catvolcano(adata,
                x = avg_log2FC,
                y = p_val,
                log2FC = 0.25,
                p_value = 0.05,
                text = s_genes)
ggsave(file.path(outdir2, 'DEG_volcano.pdf'), p, width = 5, height = 5)

#### 10.Fig6i ----
library(clusterProfiler)
library(org.Hs.eg.db)
DEG <- readRDS(file.path(outdir, "DEG_CVD_group.rds"))

#### 10.1 up ----
deg <- DEG[DEG$change == "UP",] 
top.genes <- deg[order(deg$avg_log2FC, decreasing = T),]$gene
bp <- enrichGO(
  top.genes,
  OrgDb = org.Hs.eg.db,
  keyType = 'SYMBOL',
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2)
term <- bp@result
saveRDS(term, file = file.path(outdir, "go_up.rds"))

s_term <- c('response to decreased oxygen levels',
            'response to hypoxia',
            'icosanoid metabolic process',
            'prostaglandin metabolic process',
            'BMP signaling pathway',
            'transforming growth factor beta production')
data=term[term$Description %in% s_term, ]

data$labelx=rep(0,nrow(data))
data$labely=seq(nrow(data),1)
ggplot(data, aes(x = -log10(pvalue), y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity", alpha=1, fill= "#FD9AA0", width = 0.8) +
  geom_text(aes(x=labelx, y=labely, label = data$Description), size=3.5, hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(pvalue)")+
  scale_x_continuous(expand = c(0,0))
ggsave(file = file.path(outdir2, "go_up_s.pdf"), last_plot(), height=2.5, width=3.5)

#### 10.2 down ----
deg <- DEG[DEG$change == "DOWN",] 
top.genes <- deg[order(deg$avg_log2FC, decreasing = T),]$gene
bp <- enrichGO(
  top.genes,
  OrgDb = org.Hs.eg.db,
  keyType = 'SYMBOL',
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2)
term <- bp@result
saveRDS(term, file = file.path(outdir, "go_down.rds"))

s_term <- c('endothelium development',
            'phosphatidylinositol-mediated signaling',
            'inositol lipid-mediated signaling',
            'glycerophospholipid metabolic process',
            'blood vessel endothelial cell migration',
            'sprouting angiogenesis')
data=term[term$Description %in% s_term, ]

data$labelx=rep(0,nrow(data))
data$labely=seq(nrow(data),1)
ggplot(data, aes(x = -log10(pvalue), y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity", alpha=1, fill= "#6DCCFD", width = 0.8) +
  geom_text(aes(x=labelx, y=labely, label = data$Description), size=3.5, hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  xlab("-log10(pvalue)")+
  scale_x_continuous(expand = c(0,0))
ggsave(file = file.path(outdir2, "go_donw_s.pdf"), last_plot(), height=2.5, width=3.5)

#### 11.FigS11f ----
library(GSVA)
library(GSEABase)

## GSEA
gsea.input <- readRDS(file.path(outdir, "DEG_CVD_group.rds"))
category <- "H"
genesets <- msigdbr(species = "Homo sapiens",  category = category)
genesets <- subset(genesets, select = c("gs_name", "gene_symbol"))
rownames(gsea.input) <- gsea.input$gene

gsea.input <- gsea.input[order(gsea.input$avg_log2FC, decreasing = T),]
genelist <- structure(gsea.input$avg_log2FC,names=rownames(gsea.input))
head(genelist)
tail(genelist)
res <- GSEA( genelist, TERM2GENE = genesets , pvalueCutoff=1, eps=0)
gsea_result<-as.data.frame(res@result)
gsea_result <- gsea_result[order(gsea_result$pvalue, decreasing = F),]
gsea_result$Description<- gsub('_',' ',gsea_result$Description)
gsea_result$Description<- str_to_sentence(gsea_result$Description)
# save
saveRDS(res, file.path(outdir, "H_gsea_res.rds"))
saveRDS(gsea_result, file.path(outdir, "H_gsea.rds"))

####
res <- readRDS(file.path(outdir,'H_gsea_res.rds')) 
res@result$ID <- res@result$ID|>
  str_replace_all("HALLMARK_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()
p <- res |>
  cat_gseaplot(
    c('HALLMARK_ANDROGEN_RESPONSE',
      'HALLMARK_ESTROGEN_RESPONSE_EARLY'),
    subplots = c(1, 2),
    pvalue_table = T) 
ggsave(file.path(outdir2, 'ANDROGEN_ESTROGEN.pdf'), p, height = 3, width = 5)

#### 12.Fig6j ----
#### 12.1 AUCell ----
seurat_obj <- readRDS('./files/seurat_obj_CVD.rds')
## DEGs_CRGs
DEGs_CRGs <- read.table(file = file.path(outdir, 'DEGs_CRGs.txt'))|> pull()
cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@counts) 
cells_AUC <- AUCell_calcAUC(DEGs_CRGs, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
saveRDS(AUCell_socre, './files/DEGs_CRGs_AUCscore.rds')

## C5
seurat_obj <- seurat_obj0
table(seurat_obj$disease)
cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@counts) 
category <- "C5"
genesets <- msigdbr(species = "Homo sapiens", category = category)
Dataset <- genesets %>% filter(str_detect(genesets$gs_name,"GOBP_"))
Dataset$gs_name <- Dataset$gs_name|>
  str_replace_all("GOBP_", "") |>
  str_replace_all("_", " ") |>
  str_to_sentence()
s_term <- c('Entrainment of circadian clock',
            'Regulation of bmp signaling pathway',
            'Regulation of glucose import',
            'Icosanoid biosynthetic process',
            "Prostanoid metabolic process",
            'Vascular endothelial cell proliferation',
            "Phosphatidylinositol metabolic process",
            'Endothelial cell migration',
            "Phosphatidylinositol 3 kinase signaling",
            'Oxidative phosphorylation'
)
dataset <- Dataset %>% filter(Dataset$gs_name %in% s_term)
geneSets <- lapply(unique(dataset$gs_name), 
                   function(x){dataset$gene_symbol[dataset$gs_name == x]})
names(geneSets) <- unique(dataset$gs_name)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
colnames(AUCell_socre) <- colnames(AUCell_socre)
# 
DEGs_CRGs <- readRDS('./files/DEGs_CRGs_AUCscore.rds')
AUCell_socre$DEGs_CRGs <- DEGs_CRGs$geneSet
saveRDS(AUCell_socre,file.path(outdir, paste0(category,"_AUCscore.rds")))

#### 12.2 cor ----
adata <- readRDS('./files/C5_AUCscore.rds')
feature_x_list <- c('DEGs_CRGs')
feature_y_list <- colnames(adata)[1:83]
method <- "pearson"
pathway_cor_test <- data.frame(
  feature_x=NULL,
  feature_y=NULL,
  p_value=NULL,
  estimate=NULL,
  num=NULL,
  method=NULL)

####计算
for (feature_x in feature_x_list) {
  for (feature_y in feature_y_list) {
    cor.test.res <-
      cor.test(adata[, feature_x],
               adata[, feature_y],
               method = method)
    p.value <- cor.test.res$p.value
    estimate <- cor.test.res$estimate
    pathway_cor_test <- 
      rbind(
        pathway_cor_test,
        data.frame(
          feature_x = feature_x,
          feature_y = feature_y,
          p_value = p.value,
          estimate = estimate,
          num = nrow(adata),
          method = method
        ))
  }
}
saveRDS(pathway_cor_test, "./files/C5_pathway_cor.rds")

## plot
pathway_cor_test0 <- pathway_cor_test
s_pathway <- c('Entrainment of circadian clock',
               'Regulation of bmp signaling pathway',
               'Regulation of glucose import',
               'Icosanoid biosynthetic process',
               "Prostanoid metabolic process",
               'Vascular endothelial cell proliferation',
               "Phosphatidylinositol metabolic process",
               'Endothelial cell migration',
               "Phosphatidylinositol 3 kinase signaling",
               'Oxidative phosphorylation'
)
pathway_cor_test <- pathway_cor_test0 |> filter(feature_y %in% s_pathway)

## heatmap
pathway_cor_test$text=ifelse(pathway_cor_test$p_value<0.001,"***",ifelse(pathway_cor_test$p_value<0.01,"**",ifelse(pathway_cor_test$p_value<0.05,"*","")))

library(forcats)
pathway_cor_test$feature_x <- fct_reorder(pathway_cor_test$feature_x, pathway_cor_test$estimate, .desc = TRUE)
pathway_cor_test$feature_y <- fct_reorder(pathway_cor_test$feature_y, pathway_cor_test$estimate, .desc = TRUE)
pathway_cor_test$estimate[pathway_cor_test$estimate> 0.5] = 0.5
pathway_cor_test$estimate[pathway_cor_test$estimate< -0.5] = -0.5

pdf(file = file.path(outdir2, str_c(category, "_cor_heatmap2.pdf")), width = 5.5, height = 3)
ggplot(pathway_cor_test, aes(feature_y, feature_x)) + 
  geom_tile(aes(fill = estimate), colour = "white", size = 1) + 
  scale_fill_gradient2(low = "#0099CC", mid = "white", high = "#CC0033") + 
  geom_text(aes(label = text), col = "black", size = 3) +
  theme_minimal() +    
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(size = 10, colour = "black")) +
  labs(fill = paste0("***  p<0.001", "\n", "**  p<0.01", "\n", "*  p<0.05", "\n", "\n", "Correlation")) +  
  scale_x_discrete(position = "bottom")      
dev.off()

#### 13.Fig6k ----
adata0 <- readRDS('./files/C5_AUCscore.rds')
for (i in c('ICM','DCM')) {
  C_disease <- seurat_obj@meta.data |> filter(disease == i) |> pull(ID)
  adata <- adata0[C_disease,]
  p <- adata %>%
    ggplot(aes(x = DEGs_CRGs, y = `Phosphatidylinositol metabolic process`)) +
    ggrastr::rasterise(ggpointdensity::geom_pointdensity(size = 0.5), dpi = 600) +
    geom_smooth(method = "lm", formula = y ~ x,
                color = "#6b76ae", fill = "#e5a323",
                size = 0.5, alpha = 0.2) +
    theme_cat() +labs(x='CRDGs score')+
    theme(aspect.ratio = 1, legend.margin = margin(l= -8),
          plot.title = element_text(size = 12),
          axis.title = element_text(size = 12)) +
    guides(color = guide_colorbar(
      frame.colour = "black", frame.linewidth = 0.5,
      ticks.colour = "black", title = "Density")) +
    scale_color_viridis_c() + labs(title = i)+
    ggpubr::stat_cor(method = "pearson")
  ggsave(file.path(outdir2, paste0("Phospha_",i, "_cor.pdf")), p, height=3.5, width=3.5)
}
for (i in c('ICM','DCM')) {
  C_disease <- seurat_obj@meta.data |> filter(disease == i) |> pull(ID)
  adata <- adata0[C_disease,]
  p <- adata %>%
    ggplot(aes(x = DEGs_CRGs, y = `Prostanoid metabolic process`)) +
    ggrastr::rasterise(ggpointdensity::geom_pointdensity(size = 0.5), dpi = 600) +
    geom_smooth(method = "lm", formula = y ~ x,
                color = "#6b76ae", fill = "#e5a323",
                size = 0.5, alpha = 0.2) +
    theme_cat() +labs(x='CRDGs score')+
    theme(aspect.ratio = 1, legend.margin = margin(l= -8),
          plot.title = element_text(size = 12),
          axis.title = element_text(size = 12)) +
    guides(color = guide_colorbar(
      frame.colour = "black", frame.linewidth = 0.5,
      ticks.colour = "black", title = "Density")) +
    scale_color_viridis_c() + labs(title = i)+
    ggpubr::stat_cor(method = "pearson")
  ggsave(file.path(outdir2, paste0("Prostanoid_",i, "_cor.pdf")), p, height=3.5, width=3.5)
}
for (i in c('ICM','DCM')) {
  C_disease <- seurat_obj@meta.data |> filter(disease == i) |> pull(ID)
  adata <- adata0[C_disease,]
  p <- adata %>%
    ggplot(aes(x = DEGs_CRGs, y = `Vascular endothelial cell proliferation`)) +
    ggrastr::rasterise(ggpointdensity::geom_pointdensity(size = 0.5), dpi = 600) +
    geom_smooth(method = "lm", formula = y ~ x,
                color = "#6b76ae", fill = "#e5a323",
                size = 0.5, alpha = 0.2) +
    theme_cat() +labs(x='CRDGs score')+
    theme(aspect.ratio = 1, legend.margin = margin(l= -8),
          plot.title = element_text(size = 12),
          axis.title = element_text(size = 12)) +
    guides(color = guide_colorbar(
      frame.colour = "black", frame.linewidth = 0.5,
      ticks.colour = "black", title = "Density")) +
    scale_color_viridis_c() + labs(title = i)+
    ggpubr::stat_cor(method = "pearson")
  ggsave(file.path(outdir2, paste0("ECprolif_",i, "_cor.pdf")), p, height=3.5, width=3.5)
}

#### 14.FigS11i ----
AUCell_socre <- readRDS('./files/DEGs_CRGs_AUCscore.rds')
seurat_obj$DEGs_CRGs  <- AUCell_socre$geneSet

selected_genes <- c('BMP6','BMPR1B','GATA6','GDF7','GATA4','TGFBR3', 
                    'NRG1','EDN1','FOS','BNIP3L','STC1','CXCL12','LGALS3',
                    'PLA2G5','PLA2R1','PLA2G4A','ANXA1','PTGS2','PTGS1',
                    'PDGFB','DLL4','KDR','RGCC','VEGFC','NRP1',
                    'PITPNM2','PIK3C2B','PLCB1','PLA2G2A','PIK3R3','PLCB4'
)
expr <- t(as.data.frame(seurat_obj[selected_genes,]@assays$RNA@data))
adata <- cbind(expr,seurat_obj$DEGs_CRGs)
colnames(adata)[32] <- 'CRDGs'

## correlation
method = "pearson"
pathway_cor_test <- data.frame(
  feature_x=NULL,
  feature_y=NULL,
  p_value=NULL,
  estimate=NULL,
  num=NULL,
  method=NULL)

for (feature_x in colnames(expr)) {
  for (feature_y in 'CRDGs') {
    cor.test.res <-
      cor.test(adata[, feature_x],
               adata[, feature_y],
               method = method)
    p.value <- cor.test.res$p.value
    estimate <- cor.test.res$estimate
    pathway_cor_test <- 
      rbind(
        pathway_cor_test,
        data.frame(
          feature_x = feature_x,
          feature_y = feature_y,
          p_value = p.value,
          estimate = estimate,
          num = nrow(adata),
          method = method
        ))
  }
}

expr <- pathway_cor_test
rownames(expr) <- expr$feature_x
expr <- expr[order(expr$estimate,decreasing = T),]
expr$feature_x <- factor(expr$feature_x,levels = rownames(expr))
expr <- expr |> filter(p_value < 0.05)
group_by(expr, feature_x) %>%
  summarise(estimate = max(estimate)) %>%
  ggplot(aes(feature_x, estimate)) +
  geom_segment(aes(y = 0, yend = estimate, 
                   x = feature_x, xend = feature_x)) +
  geom_point(shape = 21, colour = "black", fill="#f46f20", size = 3)+
  labs(y='Correlation with CRDGs score', x='')+
  theme(axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(face = 'italic',angle = 90, vjust = 0.5, hjust = 1,colour = "black"),
        panel.grid = element_blank())
ggsave(file.path(outdir2,'geom_CRDG_kGenes_pointplot.pdf'), last_plot(), width = 7,height = 3.5)
