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

path <- '~/PH/results/Fig3/'
setwd(path)
outdir <- './files/'
outdir2 <- './plots/'
case <- 'Fig3_'

#### 1.Fig3a ----
#### 1.1 AUCell ----
seurat_obj <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
library(msigdbr)
cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@counts) 
category <- "C5"
Dataset <- msigdbr(species = "Homo sapiens", category = category)
Dataset <- Dataset %>% filter(str_detect(Dataset$gs_name,"GOBP_"))

dataset <- Dataset
geneSets <- lapply(unique(dataset$gs_name), 
                   function(x){dataset$gene_symbol[dataset$gs_name == x]})
names(geneSets) <- unique(dataset$gs_name)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
saveRDS(cells_AUC, paste0(outdir,category, '_cells_AUC.rds'))

AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))

colnames(AUCell_socre) <- colnames(AUCell_socre)|>
  str_replace_all("_", " ") |>
  str_to_sentence()
saveRDS(AUCell_socre, paste0(outdir,category, '_AUCell_socre.rds'))

#### 1.2 radar plot ----
AUCell_socre[1:5,1:5]
s_term1 <- c("Positive regulation of leukocyte adhesion to vascular endothelial cell")
dataset1 <- AUCell_socre[,s_term1,  drop = F]

## 血管\胶原
library(data.table)
pathway_gene <- read.csv(paste0(outdir,"pathway_gene_PMC11429526.csv"))

for (i in new_name) {
  genes <- pathway_gene[pathway_gene$pathway == i,]$gene
  cells_ranking <- AUCell_buildRankings(seurat_obj@assays$RNA@counts)
  cells_AUC <- AUCell_calcAUC(genes, cells_ranking, aucMaxRank=nrow(cells_ranking)*0.1)
  AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
  metadata <- cbind(metadata,AUCell_socre)
}

dataset <- cbind(dataset1[,1], metadata[,3:4])
colnames(dataset) <- new_name
saveRDS(dataset[,-4], '~/scrna/PH/results/DEG/files/AUCell/AUC_pathways.rds')

## plot
library(ggradar)

seurat_obj$ID <- rownames(seurat_obj@meta.data)
meta <- dataset
meta <- as.data.frame(meta)
meta$cell_type <- seurat_obj$classic1
object <-
  aggregate(meta[, -length(colnames(meta))], list(meta$cell_type), FUN = mean)
colnames(object)[1] <- 'cell type'

for (i in 2:4) {
  object[,i] <- (object[,i]-min(object[,i]))/(max(object[,i])-min(object[,i]))
}
rownames(object) <- object[,1]
object <- object[,-1]
object <- as.data.frame(t(object))
object$pathways <- rownames(object)
object <- object[,c(6,1:5)]

ggradar(object, background.circle.transparency = 0,
        plot.extent.x.sf = 1.5,
        base.size = 3,axis.label.size = 4,grid.label.size = 4,
        group.point.size = 2,
        group.line.width = 1,
        legend.position = 'right',
        legend.text.size = 10) +
  theme(plot.title = element_text(hjust = 0.2,size = 12))
ggsave(paste0(outdir2,'radar_classic1_pathways_new.pdf'),last_plot(),width = 10,height = 5)

#### 2.Fig3b ----
type_colors <- c("#f2ccac","#edae11","#FE870D","#e1abbc","#F88A89","#c04932",
                 "#c1d5b9","#99CD91","#009E73","#6a73cf","#916CB2","#237BB2","#8A9D5B")
meta <- dataset
meta <- as.data.frame(meta)
meta$cell_type <- seurat_obj$cell_type
object <- aggregate(meta[, -length(colnames(meta))], list(meta$cell_type), FUN = mean)
colnames(object)[1] <- 'cell_type'
object[,1] <- as.character(object[,1])
object <- object[c(1:11),]

for (i in 2:4) {
  case <- new_name[i-1]
  data <- object[,c(1,i)]
  colnames(data)[2] <- 'val'
  
  p <- ggplot(data, aes(x=reorder(cell_type,val), y=round(val,4), fill=cell_type)) +
    geom_bar(stat = 'identity', width = 0.7, fill = type_colors) +
    labs(y = case,x = '') +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    theme(axis.text.x = element_text(colour = 'black', size = 10),
          axis.text.y = element_text(colour = 'black', size = 10))+
    coord_flip()
  ggsave(paste0(outdir2,case,'_bar_pathways_type.pdf'),last_plot(),width = 5,height = 5)
}

#### 3.Fig3c ----
seurat_obj0 <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
seurat_obj <- seurat_obj0 |> subset(cell_type %in% c("Art.c01.DKK2", "Art.c02.NEBL", "Art.c03.VEGFA",
                                                     "Cap.c04.PAPSS2", "Cap.c05.RGCC", "Cap.c06.RDH10",
                                                     "Cap.c07.CA2", "Cap.c08.TMEM163", "Cap.c09.RGS6",
                                                     "Ven.c10.VCAM1", "Ven.c11.ACTA2"))
prop0 <- as.data.frame(prop.table(table(seurat_obj$cell_type, seurat_obj$group),
                                  margin = 2) * 100)
prop <- data.frame(Normal=prop0$Freq[1:11],
                   CVD=prop0$Freq[14:24])
rownames(prop) <- prop0$Var1[1:11]
prop$ratio <- (prop$CVD)/(prop$Normal)
prop$group <- 'CVD'
prop$group[prop$ratio < 1] <- 'Normal'
prop$ratio2 <- ifelse(prop$ratio >= 1, prop$ratio, 1 / prop$ratio)
prop$celltype <- rownames(prop)

prop <- prop %>%
  arrange(group, desc(ratio)) %>%
  mutate(celltype = factor(celltype, levels = unique(celltype)))
prop <- prop |> mutate(ratio2_adj = ifelse(group == "CVD", ratio2, 2 - ratio2)) |> na.omit()

## 绘制图表
ggplot(prop, aes(celltype, ratio2_adj, fill = celltype)) +
  labs(y = 'Relative ratio of CVD') +
  geom_segment(aes(y = 1, yend = ratio2_adj, x = celltype, xend = celltype)) +
  geom_point(shape = 21, colour = "white", size = 3) +
  scale_y_continuous(breaks = c(-0.5, 0, 0.5, 1, 1.5, 2, 2.5),
                     labels = c('2.5', "2", "1.5", "1", "1.5","2", '2.5'),
                     limits = c(-0.5, 2.5)) +
  theme_classic() +
  theme(aspect.ratio = 0.7,
        axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("#237BB2","#6a73cf","#009E73","#e1abbc","#99CD91","#8A9D5B",
                               "#F88A89","#c1d5b9","#edae11","#c04932","#FE870D","#f2ccac","#916CB2")) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+NoLegend()
ggsave(file.path(outdir2,'Ratio_N_CVD_Pointplot.pdf'),
       ggplot2::last_plot(), height=3.5,  width=3.5)

#### 4.Fig3d ----
seurat_obj <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
Idents(seurat_obj) <- seurat_obj$classic1
seurat_obj <- subset(seurat_obj, idents = 'VenECs')
table(seurat_obj$cell_type)

av_expr_v10 <- AverageExpression(
  seurat_obj[,seurat_obj$cell_type== 'Ven.c10.VCAM1'], assays = "RNA",
  group.by = "disease" )[["RNA"]] %>% 
  as.data.frame()
av_expr_v11 <- AverageExpression(
  seurat_obj[,seurat_obj$cell_type== 'Ven.c11.ACTA2'], assays = "RNA",
  group.by = "disease" )[["RNA"]] %>% 
  as.data.frame()

features= c('NFKB1', 'CXCL10', 'SELP','VCAM1','JAM2','VEGFC','IL1R1','IL6ST','NAIP', 
            'SEMA6A','GAB1','NRP2','TNXB', 'TGFBR2', 'SMAD1', 'FBN1' )
expr <- cbind(av_expr_v10[features,], av_expr_v11[features,])
annotation_col <-  data.frame(Disease = factor(colnames(expr)), Celltype = c(rep('Ven_c10_VCAM1', 9), rep('Ven_c11_ACTA2', 9)))
colnames(expr) <- c(1:18)

rownames(annotation_col) <- c(1:18)
ann_colors = list(
  Celltype = c(Ven_c10_VCAM1 = "#6a73cf", Ven_c11_ACTA2 = "#916CB2"),
  Disease = c(Normal="#B3DE69",ACM="#FB8072", DCM="#80B1D3", HCM="#FDB462", ICM="#8DD3C7",
              AMI="#FCCDE5", CHD="#CCEBC5", CS="#377EB8", COV="#1B9E77"))

pdf(file = paste0(outdir2, "v10_gene_heatmap2.pdf"), width = 4.5,height = 3.5)
pheatmap(expr, annotation_colors = ann_colors,
         color = colorRampPalette(c("#3980b8","white","#ef3b3c"))(256),
         scale = 'row',annotation_col = annotation_col,
         cellwidth = 8, cellheight = 8, 
         border_color = "white", fontsize = 8, 
         cluster_rows = T, cluster_cols = F,
         show_rownames = T, show_colnames = F,
         treeheight_row = 10,
         clustering_method = "average",gaps_col = 10)
dev.off()

#### 5.FigS7d ----
disease_colors <- c( "#FB8072", "#80B1D3", "#FDB462", "#8DD3C7", "#FCCDE5",  "#CCEBC5", "#377EB8", "#1B9E77")
names(disease_colors) <- c('ACM', 'DCM', 'HCM', 'ICM', 'AMI', 'CHD', 'CS', 'COV')
group_colors <- c("high" = "#f46f20", "low" = "#156077")

## EC
Endo <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
prop <- as.data.frame(prop.table(table(Endo$cell_type, Endo$sample),margin = 2) * 100)
colnames(prop) <- c('cell_type_EC', 'sample', 'prop_EC')

prop2 <- prop[prop$cell_type_EC == "Ven.c10.VCAM1", ]
c10_median <- median(prop2$prop_EC)
prop2$prop_group <- "high"
prop2$prop_group[prop2$prop_EC < c10_median] <- "low"
write.table(prop2, file.path(outdir, 'Endo_prop.csv'), sep = ',')

## TY4
TY4 <- qread('~/PH/results/files/Cell4_CVD.qs')
prop_A <- as.data.frame(prop.table(table(TY4$cell_type, TY4$sample), margin = 2) * 100)
colnames(prop_A) <- c('cell_type_TY4', 'sample', 'prop_TY4')
prop_B <- prop_A |> filter(sample %in% prop2$sample)
prop_C <- merge(prop_B, prop2, by = 'sample', all.x = T)
info <- Endo@meta.data %>%
  distinct(sample, .keep_all = TRUE)
prop_D <- prop_C %>%
  left_join(info %>% select(sample, disease), by = "sample")
prop_D$prop_group <- factor(prop_D$prop_group, levels = c('low', 'high'))
prop_D$disease <- factor(prop_D$disease, levels = c('ACM', 'DCM', 'HCM', 'ICM', 'AMI', 'CHD', 'CS', 'COV'))
write.table(prop_D, file.path(outdir, 'TY4_prop.csv'), sep = ',')

## plot
library(ggpubr)
ctype <- names(table(prop_D$cell_type_TY4))
for(i in ctype){
  data <- prop_D[prop_D$cell_type_TY4 == i, ]
  ggplot(data, aes(x = prop_group, y = prop_TY4)) +
    geom_jitter(aes(color = disease), width = 0.2, size = 1) +
    geom_violin(aes(color = prop_group), width = 0.6, fill = NA) +
    geom_boxplot(aes(color = prop_group), width = 0.05, outlier.shape = NA) +
    scale_color_manual(values = c(disease_colors, group_colors)) + 
    labs(title = i, x = "", y = "Frequence") +
    theme_classic() +
    theme(axis.text = element_text(colour = 'black'), legend.position = "none") +
    stat_compare_means(method = "wilcox.test")
  ggsave(file.path(outdir2, paste0(i, '_group.pdf')),
         last_plot(), width = 3,height = 3)
}

#### 6.FigS7f_monocle3 ----
seurat_obj <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
seurat_obj <- subset(seurat_obj, classic1 == "VenECs")
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
  cds <- reduce_dimension(cds, preprocess_method = "PCA")
  
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
                                                                     cell_phenotype="cell_type",
                                                                     root_types=c('Ven.c10.VCAM1')))
  saveRDS(cds, file.path(outdir, "monocle3.rds"))
}
mono3(seurat_obj)

## plot
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
cds <- readRDS(file.path(outdir, "monocle3.rds"))
cds2 <- order_cells(cds, root_pr_nodes= 
                      get_earliest_principal_node(cds,
                                                  cell_phenotype="cell_type",
                                                  root_types=c('Ven.c10.VCAM1')))
p <- plot_cells(cds2, reduction_method="UMAP",
                show_trajectory_graph = T, color_cells_by="pseudotime",
                group_label_size = 7, trajectory_graph_color = 'white',
                label_groups_by_cluster=FALSE, label_branch_points=FALSE,
                label_leaves=FALSE, label_cell_groups = FALSE) + 
  ggtitle('Monocle3')+theme(aspect.ratio = 1)
ggsave(file.path(outdir2, "monocle3_pseudotime.pdf"), p, height=4, width=4)

#### 7.Fig3f ----
#### 7.1 CytoTRACE ----
library(CytoTRACE)

table(is.na(seurat_obj$cell_type))
table(seurat_obj$cell_type)
phe <- seurat_obj$cell_type
phe = as.character(phe)
names(phe) <- rownames(seurat_obj@meta.data)

mat_3k <- as.matrix(seurat_obj@assays$RNA@counts)
mat_3k[1:4,1:4]
results <- CytoTRACE(mat = mat_3k)
saveRDS(results, file.path(outdir, 'CytoTRACE_results.rds'))

phenot <- seurat_obj$cell_type
phenot <- as.character(phenot)
names(phenot) <- rownames(seurat_obj@meta.data)
emb <- seurat_obj@reductions[["umap"]]@cell.embeddings
plotCytoTRACE(results, phenotype = phenot, emb = emb, outputDir = outdir2)
plotCytoGenes(results, numOfGenes = 10, outputDir = outdir2)
plotCytoTRACE(results, phenotype = phenot,emb = emb,gene = "VEZF1",outputDir = outdir2)
##
data <- data.frame(TRACE=results$CytoTRACE,
                   Rank=results$CytoTRACErank,
                   GCS=results$GCS,
                   Counts=results$Counts)
data$Cells <- rownames(data)
data$group <- seurat_obj$group
data$disease <- seurat_obj$disease
genes <- c("HOXB2", "KLF16", "SRF", "VEZF1", "EGR1",
           'NR2F2', 'TEK', 'ANGPT1', 'ANGPT2',
           "RELA",  "RFX5", "CEBPB", "TGIF1",  "ATF4",
           "SOX4","CEBPG", 'LBX2', 'IRF5', 'NR1H4')

#### 7.2 TF ----
data <- subset(seurat_obj, idents = 'Ven.c10.VCAM1')
data <- hp_run_pyscenic(x = data,
                        species = "human",
                        outdir = './files/TF/Ven.c10.VCAM1/')

####
c10 <- subset(seurat_obj, cell_type == "Ven.c10.VCAM1")
TF_G <- read.table('./files/TF/Ven.c10.PLVAP/tfs_target.tsv',sep = ',',header = T) 
TF_G2 <- TF_G[!duplicated(TF_G[c('tf')]), ]
regulonAUC <- importAUCfromText('./files/TF/Ven.c10.VCAM1/auc_mtx.csv')
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
cellInfo <- c10@meta.data

##
c10[["scenic"]] <- CreateAssayObject(counts = t(getAUC(regulonAUC)))
Idents(c10) <- 'disease'
regulonActivity_byCellType <-
  sapply(split(rownames(cellInfo), Idents(c10)),
         function(cells)
           rowMeans(t(getAUC(regulonAUC))[, cells]))
rownames(regulonActivity_byCellType) <- TF_G2$symbol

##
tf.name <- rownames(regulonActivity_byCellType)
data <- regulonActivity_byCellType[tf.name,]
data_clean <- data[rowSums(data != 0) > 0, ]
TF_s <- c("VEZF1 (386g)", "NR1H4 (69g)", "HOXB2 (1257g)", "KLF16 (403g)", "SRF (69g)",
          "RFX5 (2100g)", "TGIF1 (7g)",  "ATF4 (858g)", "SOX4 (20g)", "CEBPG (34g)")
pdf(file.path(outdir2, "TF_gene_disease.pdf"), width = 3, height = 2)
pheatmap(data_clean[TF_s,-7],
         scale = "row", cluster_cols = F, cluster_rows = F,
         color = c(colorRampPalette(colors = c("#2166AC","#4393C3","#92C5DE","#D1E5F0"))(50),
                   colorRampPalette(colors = c("#FDDBC7","#F4A582","#D6604D","#B2182B"))(50)),
         border_color = "white", fontsize = 5, angle_col = "45",
         cellheight = 8, cellwidth = 12,
         main = 'Transcription Factor')
dev.off()

#### 8.FigS7g ----
regulonAUC <- importAUCfromText('./files/TF/Ven.c10.VCAM1/auc_mtx.csv')
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]

cellTypes <- c10@meta.data[,c('disease',"orig.ident")]
rss <- calcRSS(AUC=t(getAUC(regulonAUC)), 
               cellAnnotation=cellTypes[rownames(regulonAUC), 'disease'])
rownames(rss) <- gsub(',','',rownames(rss))
rownames(rss) <- TF_G2$symbol
rss <- na.omit(rss)

rssPlot <- plotRSS(rss,zThreshold = 0.5)
plotly::ggplotly(rssPlot$plot)

## Rank plot
data <- as.data.frame(rss)
data$TFs <- rownames(data)
data$gene_name = rownames(data)
data$rss <- data$Normal
data <- data[order(data$Normal,decreasing = T),]
data$rank <- 1:length(data$gene_name)
data <- data[,c('rss','rank')]
head(data, n=50) |> rownames()
pointed_out <- data |> filter(rownames(data) %in% 
                                c("HOXB2 (1257g)", "E2F1 (85g)", "KLF16 (403g)", "SRF (69g)", "VEZF1 (386g)", "EGR1 (2003g)",
                                  "RELA (1685g)",  "RFX5 (2100g)", "CEBPB (583g)", "TGIF1 (7g)",  "ATF4 (858g)", "ELK1 (316g)",
                                  "ZNF200 (627g)", "SOX4 (20g)", "RARG (519g)", "IRF5 (178g)","CEBPG (34g)"))
ggplot(data, aes(x=rank, y=rss))+ 
  geom_point(size=3,color="grey50",alpha =0.4)+theme_minimal()+
  geom_point(data = pointed_out, stroke = 0.5, size=4, shape=16, color="#DC050C",alpha =0.6)+
  ggrepel::geom_text_repel(data=pointed_out, aes(label=rownames(pointed_out)), color="black",
                           size=4, fontface="italic", segment.size=0.5, nudge_x=60,
                           direction="y", hjust=0,ylim = c(0.06, 0.50))+ theme_classic()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title = element_text(colour = 'black', size = 15),
        axis.text = element_text(colour = 'black', size = 12))+
  labs(x='Rank', y='RSS')
ggsave(file.path(outdir2,'TF_rank_all.pdf'), ggplot2::last_plot(), width = 6,height = 6)

#### 9.Fig3i ----
s_gene <- TF_G |> filter(symbol == 'VEZF1 (386g)') |> pull('target_gene')
bp <-
  enrichGO(
    s_gene,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2)
term <- bp@result
s_term <- c('heart process',
            'cytoplasmic translation',
            'endothelial cell development',
            'endothelial cell differentiation',
            'stem cell differentiation'
)
df <- term |> filter(Description %in% s_term)
df$pvalue <- as.numeric(df$pvalue)
df$Count <- as.numeric(df$Count)
df$labelx=rep(0,nrow(df))
df$labely=seq(nrow(df),1)
ggplot(data = df, aes(Count, reorder(Description,Count), fill = pvalue)) +
  geom_bar(stat="identity", alpha=1, width = 0.8) + 
  scale_fill_gradient(low = "#ee6470", high = "#ee9ca7" ) +
  geom_text(aes(x = labelx, y = Description, label = Description),
            size=3.5,  hjust =0)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,size=14),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(colour = 'black',size = 10,vjust = 3),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  labs(fill = "-log10(pvalue)")+xlab("Count")+
  scale_y_discrete(labels=function(x) str_wrap(x, width=40))
ggsave(file.path(outdir2, 'VEZF1_Barplot_s.pdf'), last_plot(), height=2.5, width=4)

#### 10.Fig3h ----
results <- readRDS(paste0(outdir, 'CytoTRACE_results.rds'))
data <- data.frame(TRACE=results$CytoTRACE,
                   Rank=results$CytoTRACErank,
                   GCS=results$GCS,
                   Counts=results$Counts)
data$Cells <- rownames(data)
data$group <- seurat_obj$group
data$disease <- seurat_obj$disease

seurat_obj$GENE <- seurat_obj@assays$RNA@data['VEZF1',]
data$GENE <- seurat_obj$GENE
head(data)
## group
ggplot(data, aes(x = -Rank, y = GENE, color = group)) +
  geom_smooth(method = "loess", se = T) + 
  labs(title = 'VEZF1', x = "Pseudotime", y = "Expression") +
  theme_classic() + scale_color_manual(values = c("darkgreen","#A20E04"))+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(color = "black", size = 8)) 
ggsave(file.path(outdir2,paste0('VEZF1', '_expr_tra.pdf')), last_plot(), width = 4, height = 4)
## disease
data <- data |> filter(data$disease %in% c('Normal','DCM','HCM','ICM','AMI','CHD','COV'))
ggplot(data, aes(x = -Rank, y = GENE, color = disease)) +
  geom_smooth(method = "loess", se = T) + 
  labs(title = 'VEZF1', x = "Pseudotime", y = "Expression") +
  theme_classic() + scale_color_manual(values = disease_color)+
  theme(aspect.ratio = 1,
        panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(color = "black", size = 8)) 
ggsave(file.path(outdir2,paste0('VEZF1', '_expr_tra.pdf')), last_plot(), width = 4, height = 4)
