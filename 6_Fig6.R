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

#### 2.Fig6b Vlnplot ----
library(ggraph)
CRGs <- read.table('./files/CRGs_35436363.txt') |> pull()
CRGs <- CRGs[(CRGs %in% rownames(seurat_obj))]

cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@counts) 
cells_AUC <- AUCell_calcAUC(CRGs, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
saveRDS(AUCell_socre, file.path(outdir, 'CRD score celltype.rds'))
seurat_obj$AUC  <- AUCell_socre
seurat_obj$AUC[seurat_obj$AUC> 0.25] <- 0.25
seurat_obj$AUC[seurat_obj$AUC< 0.20] <- 0.20
data <- data.frame(seurat_obj@meta.data, seurat_obj@reductions$umap@cell.embeddings)
ggviolin(data, x='cell_type',y="AUC",
         xlab = '', ylab = "CRDGs score",
         fill='cell_type', add = 'boxplot',
         add.params = list(fill='white', width=0.15))+ 
  scale_fill_manual(values = type_colors)+
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) + NoLegend()
ggsave(file.path(outdir2, 'CRD score vlnplot.pdf'),
       ggplot2::last_plot(),width = 6, height = 3.5)

#### 3.Fig6c disease ----
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

#### 6.1 barplot ----
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

#### 7.Fig6h ----
group <- seurat_obj@meta.data$CVD
plot_data <- data.frame(log2_NR1D1_ARNTL = seurat_obj$log2_NR1D1_ARNTL, group = group)

ggplot(plot_data, aes(x = group, y = log2_NR1D1_ARNTL, color = group)) +
  geom_boxplot(width = 0.6, fill = NA,
               outlier.color = NA, linewidth = 1) +
  theme_classic2() +
  labs(y = "Log2(NR1D1/ARNTL)") +
  # geom_jitter(width = 0.2, size = 0.6, alpha = 0.5) +
  scale_color_manual(values = c("#D55E00", "#194A55")) +
  theme(axis.text = element_text(color = 'black'), axis.title.x = element_blank())+
  stat_compare_means(comparisons = list(c('CRDGs low', 'CRDGs high')), vjust = 2)+NoLegend()

ggsave(file.path(outdir2, "Box_Log2_CRD.pdf"),
       ggplot2::last_plot(), height=3, width=2.5)

#### 8.Fig6i ----
#### 8.1 DEG ----
seurat_obj0 <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
sce <- subset(seurat_obj0, idents = c("LEC.c12.FLT4","Endo.c13.NPR3"), invert = T)
DefaultAssay(sce) <- "RNA"
Idents(sce) <- 'CVD'
sce$Group <- sce$CVD
bs = split(colnames(sce),sce$sample)
ct = do.call(
  cbind,lapply(names(bs), function(x){ 
    kp =colnames(sce) %in% bs[[x]]
    rowSums( as.matrix(sce@assays$RNA@counts[, kp]  ))
  })
)
colnames(ct) <- names(bs)
saveRDS(ct, file.path(outdir, 'bulk_sample_expr.rds'))
phe = unique(sce@meta.data[,c('sample','Group')])
saveRDS(phe, file.path(outdir, 'bulk_metadata.rds'))

##
exp0 <- readRDS(file.path(outdir, 'bulk_sample_expr.rds'))
metadata <- readRDS(file.path(outdir, 'bulk_metadata.rds'))
exp <- exp0[, metadata$sample]
Group <- metadata |> pull(Group)
Group <- factor(Group, levels = c('CVDlow', 'CVDhigh'))

##
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

pvalue.cutoff <- 0.05
logFC.cutoff <- 0.58
DEG_edgeR$change <-
  as.factor(ifelse(
    DEG_edgeR$FDR < pvalue.cutoff &
      abs(DEG_edgeR$logFC) > logFC.cutoff,
    ifelse(DEG_edgeR$logFC > logFC.cutoff , 'UP', 'DOWN'),
    'NOT' ))
DEG_edgeR$cluster <- 'CVD'
saveRDS(DEG_edgeR, file.path(outdir, "DEG_CVD_group.rds"))

#### 8.2 volcano plot ----
library(ggrepel)
library(ggrastr)

s_genes <- c('NLGN1','PPP1CB','PER1','NAMPT','MMP19', 
             'BMP6','BMPR1B','GATA6','GDF7','GATA4','TGFBR3',
             'PLA2G5','GGT7','PLA2G4A','ANXA1','PTGS2','PTGS1',
             'HIF1A','EDN1','LMNA','DLL4','KDR','RGCC','NRP1',
             'HLA-DRB1','CD36','PLA2G2A','PLCB1','PITPNM2','UNC5B','PDGFB','KDR','CA8')
adata <- readRDS(file.path(outdir, "DEG_CVD_group.rds"))
s_genes <- intersect(s_genes, adata|> filter(change %in% c('UP', 'DOWN'))|> pull(Gene_Symbol))
adata %>% glimpse()
adata$label <- ""
adata$label[adata$Gene_Symbol %in% s_genes] <- adata[adata$Gene_Symbol %in% s_genes, ]$Gene_Symbol
p <- catvolcano(adata,
                x = logFC,
                y = FDR,
                log2FC = 0.58,
                p_value = 0.05,
                text = s_genes)
ggsave(file.path(outdir2, 'DEG_volcano.pdf'), p, width = 5, height = 5)

#### 9.Fig6j ----
library(clusterProfiler)
library(org.Hs.eg.db)
DEG <- readRDS(file.path(outdir, "DEG_CVD_group.rds"))
DEG$gene <- DEG$Gene_Symbol

## UP
deg <- DEG[DEG$change == "UP",] 
top.genes <- deg[order(deg$logFC, decreasing = T),]$gene
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

s_term <- c('transforming growth factor beta production',
            'response to decreased oxygen levels',
            'icosanoid biosynthetic process',
            'extrinsic apoptotic signaling pathway',
            'prostaglandin metabolic process',
            'response to BMP')
data=term_up[term_up$Description %in% s_term, ]
data$labelx=rep(0,nrow(data))
data$labely=seq(nrow(data),1)
ggplot(data, aes(x = -log10(pvalue),
                 y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity", alpha=1,
           fill= "#FD9AA0", width = 0.8) +
  geom_text(aes(x=labelx, y=labely,
                label = data$Description),
            size=3.5, hjust =0)+
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
ggsave(file = file.path(outdir2, "go_up_s.pdf"),
       ggplot2::last_plot(),height=2.5, width=3.5)

## DOWN
deg <- DEG[DEG$change == "DOWN",] 
head(deg)
top.genes <-
  deg[order(deg$logFC, decreasing = T),]$gene
bp <-
  enrichGO(
    top.genes,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
term <- bp@result
saveRDS(term, file = file.path(outdir, "go_down.rds"))

s_term <- c('activation of innate immune response',
            'endothelium development',
            'endothelial cell differentiation',
            'phosphatidylinositol-mediated signaling',
            'endothelial cell migration',
            'regulation of endothelial cell proliferation')
data=term_down[term_down$Description %in% s_term, ]
data$labelx=rep(0,nrow(data))
data$labely=seq(nrow(data),1)
ggplot(data, aes(x = -log10(pvalue),
                 y = reorder(Description,-log10(pvalue))))  +
  geom_bar(stat="identity", alpha=1,
           fill= "#6DCCFD", width = 0.8) +
  geom_text(aes(x=labelx, y=labely,
                label = data$Description),
            size=3.5, hjust =0)+
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
ggsave(file = file.path(outdir2, "go_donw_s.pdf"),
       ggplot2::last_plot(),height=2.5, width=3.5)

#### 10.Fig6k ----
category <- 'C5'
outdir5 <- paste0('./files/group/')
dir.create(outdir5,recursive = T)
outdir6 <- paste0("./plots/group/")
dir.create(outdir6,recursive = T)
seurat_obj <- qread('~/project/PH/filter_name/results/SCT15_1_P2/files/named/classic/Endocardial endo.qs')
Idents(seurat_obj) <- 'disease'

adata <- readRDS('~/project/PH/signature/CVD/results/files/C5_AUCscore.rds')
colnames(adata)
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
saveRDS(pathway_cor_test, "~/project/PH/signature/CVD/results/files/group/C5_pathway_cor.rds")

pathway_cor_test0 <- readRDS("~/project/PH/signature/CVD/results/files/group/C5_pathway_cor.rds")
s_pathway <- c('Entrainment of circadian clock',# positive
               'Regulation of bmp signaling pathway',
               'Regulation of glucose import',
               'Icosanoid biosynthetic process',
               "Prostanoid metabolic process",
               'Vascular endothelial cell proliferation',# negtive
               "Phosphatidylinositol metabolic process",
               'Endothelial cell migration',
               "Phosphatidylinositol 3 kinase signaling",
               'Oxidative phosphorylation'
)
pathway_cor_test <- pathway_cor_test0 |> filter(feature_y %in% s_pathway)
pathway_cor_test$text=ifelse(pathway_cor_test$p_value<0.001,"***",ifelse(pathway_cor_test$p_value<0.01,"**",ifelse(pathway_cor_test$p_value<0.05,"*","")))
library(forcats)
pathway_cor_test$feature_x <- fct_reorder(pathway_cor_test$feature_x, pathway_cor_test$estimate, .desc = TRUE)
pathway_cor_test$feature_y <- fct_reorder(pathway_cor_test$feature_y, pathway_cor_test$estimate, .desc = TRUE)
pathway_cor_test$estimate[pathway_cor_test$estimate> 0.5] = 0.5
pathway_cor_test$estimate[pathway_cor_test$estimate< -0.5] = -0.5

pdf(file = file.path(outdir6, str_c(category, "_cor_heatmap2.pdf")), width = 5.5, height = 3)
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

#### 11.Fig6l + FigS11g ----
outdir5 <- paste0('./files/group/disease/')
dir.create(outdir5,recursive = T)
outdir6 <- paste0("./plots/group/disease/")
dir.create(outdir6,recursive = T)
source("/home/xingwl/share/20220802_scrna-m6A/custom_plot_function.R")
seurat_obj <- qread('~/project/PH/filter_name/results/SCT15_1_P2/files/named/classic/Endocardial endo.qs')
adata0 <- readRDS('~/project/PH/signature/CVD/results/files/C5_AUCscore.rds')
disease <- names(table(seurat_obj$disease))[2:10]

for(i in c('ICM', 'DCM', 'COV')){
  C_disease <- seurat_obj@meta.data |> filter(disease == i) |> pull(ID)
  adata <- adata0[C_disease,]
  p3 <- adata %>%
    ggplot(aes(x = DEGs_CRGs, y = `Phosphatidylinositol metabolic process`)) +
    ggrastr::rasterise(ggpointdensity::geom_pointdensity(size = 0.5), dpi = 600) +
    geom_smooth(method = "lm", formula = y ~ x,
                color = "#6b76ae", fill = "#e5a323",
                size = 0.5, alpha = 0.2) +
    theme_cat() +labs(x='CRDGs score')+
    theme(aspect.ratio = 1, legend.margin = margin(l= -8),
          axis.text = element_text(size = 12),
          plot.title = element_text(size = 12),
          axis.title = element_text(size = 12)) +
    guides(color = guide_colorbar(
      frame.colour = "black", frame.linewidth = 0.5,
      ticks.colour = "black", title = "Density")) +
    scale_color_viridis_c() + labs(title = i)+
    ggpubr::stat_cor(method = "pearson")
  ggsave(file.path(outdir6, paste0("Phosphatidy_",i, "_cor.pdf")), 
         p3, height=3.5, width=3.5)
}

# Prostanoid metabolic process
# 'Vascular endothelial cell proliferation',
# "Phosphatidylinositol metabolic process",


