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
source("~/Collection/code/scrna-seq.R")
source("~/Collection/code/plot.R")

path <- '~/PH/results/Fig2/'
setwd(path)
outdir <- './files/'
outdir2 <- './plots/'
case <- 'Fig2_'

#### 1.Fig2a ----
DEG <- readRDS('~/PH/results/DEG/files/DEGs_disease_number.rds') |> 
  filter(change == 'Upregulated')
DEG$cluster <- DEG$disease
DEG$gene <- DEG$Gene_Symbol
DEG_split <- split(DEG, DEG$cluster)
for(i in 1:length(table(DEG$cluster))){
  DEG_split[[i]] <- DEG_split[[i]]$gene
}

upsetData=fromList(DEG_split)
color_ct <- c("#FCCDE5", "#8DD3C7", "#80B1D3", "#FDB462",
              "#CCEBC5", "#1B9E77", "#FB8072", "#377EB8")
pdf(file=file.path(outdir2, 'UP_upset.pdf'),onefile = FALSE,width=5.5,height=4)
upset(upsetData,
      nsets = length(DEG_split),
      nintersects = 24,
      order.by = "freq",
      show.numbers = "yes",
      number.angles = 0,
      point.size = 3,
      matrix.color="#b0b9b8",
      line.size = 0.8,
      mainbar.y.label = "Intersections of up-regulated genes",
      sets.x.label = "Set Size",
      sets.bar.color = color_ct,
      main.bar.color ="black",
      queries = list(
        list(query=intersects, params=list("ACM","DCM","HCM","ICM", "AMI","CHD","CS","COV"), color= "#fce38a", active=T),
        list(query=intersects, params=list("AMI", 'ICM'), color="#fce38a", active=T),
        list(query=intersects, params=list("AMI"), color="#FCCDE5", active=T)
      )
)  
dev.off()

#### barplot
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)
library(enrichplot)

genelist <- Reduce(intersect, DEG_split)
title <- paste0('geneset (', length(genelist), ')')
bp <-
  enrichGO(
    genelist,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2)
term <- bp@result
saveRDS(term, paste0(outdir, title, '_terms_GO.rds'))

terms <- c('histone modification',
           'dephosphorylation',
           'small GTPase mediated signal transduction',
           'transforming growth factor beta receptor signaling pathway',
           'miRNA metabolic process',
           'Wnt signaling pathway',
           'phosphatidylinositol dephosphorylation',
           'stress-activated MAPK cascade')
data <- term[term$Description %in% terms,]
data$labelx=rep(0,nrow(data))
data$labely=seq(nrow(data),1)
ggplot(data, aes(x = -log10(pvalue),
                 y = reorder(Description,-log10(pvalue)))) +
  geom_bar(stat="identity", alpha=1,
           fill= "#fce38a", width = 0.8) +
  geom_text(aes(x=labelx, y=labely,
                label = data$Description),
            size=4, hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  labs(title = paste0('geneset (', length(genelist), ')'), 
       x="-log10(pvalue)")+
  scale_x_continuous(expand = c(0,0))
ggsave(file.path(outdir2, paste0(title, '_Barplot.pdf')),
       last_plot(), height=3.5, width=4.5)

#### 2.Fig2b ----
DEG <- readRDS('~/PH/results/DEG/files/DEGs_disease_number.rds') |> 
  filter(change == 'Downregulated')
DEG$cluster <- DEG$disease
DEG$gene <- DEG$Gene_Symbol
DEG_split <- split(DEG, DEG$cluster)
for(i in 1:length(table(DEG$cluster))){
  DEG_split[[i]] <- DEG_split[[i]]$gene
}

upsetData=fromList(DEG_split)
color_ct <- c("#FDB462", "#CCEBC5", "#1B9E77","#FB8072", "#80B1D3",
              "#FCCDE5", "#8DD3C7", "#377EB8" )
pdf(file=file.path(outdir2, 'DOWN_upset.pdf'),onefile = FALSE,width=5.5,height=4)
upset(upsetData,
      nsets = length(DEG_split),
      nintersects = 24,
      order.by = "freq",
      show.numbers = "yes",
      number.angles = 0,
      point.size = 3,
      matrix.color="#b0b9b8",
      line.size = 0.8,
      mainbar.y.label = "Intersections of down-regulated genes",
      sets.x.label = "Set Size",
      sets.bar.color = color_ct,
      main.bar.color ="black",
      queries = list(
        list(query=intersects, params=list("ACM","DCM","HCM","ICM", "AMI","CHD","CS","COV"), color= "#fce38a", active=T),
        list(query=intersects, params=list("HCM"), color="#FDB462", active=T)
      )
)  
dev.off()

#### barplot
genelist <- Reduce(intersect, DEG_split)
title <- paste0('geneset (', length(genelist),')')
bp <-
  enrichGO(
    genelist,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2)
term <- bp@result
saveRDS(term, paste0(outdir, title, '_terms_GO.rds'))

terms <- c('ribonucleoprotein complex biogenesis',
           'mitochondrial ATP synthesis coupled electron transport',
           'glycolytic process',
           'aerobic respiration',
           'electron transport chain',
           'regulation of protein stability',
           'NADH dehydrogenase complex assembly',
           'RNA splicing')
data <- term[term$Description %in% terms,]
data$labelx=rep(0,nrow(data))
data$labely=seq(nrow(data),1)
ggplot(data, aes(x = -log10(pvalue),
                 y = reorder(Description,-log10(pvalue)))) +
  geom_bar(stat="identity", alpha=1,
           fill= "#fce38a", width = 0.8) +
  geom_text(aes(x=labelx, y=labely,
                label = data$Description),
            size=4, hjust =0)+
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  labs(title = paste0('geneset (', length(genelist),')'), 
       x="-log10(pvalue)")+
  scale_x_continuous(expand = c(0,0))
ggsave(file.path(outdir2, paste0(title, '_Barplot.pdf')),
       last_plot(), height=3.5, width=4.5)

#### 3.Fig2c ----
seurat_obj0 <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
sce <- subset(seurat_obj0, idents = c("LEC.c12.FLT4","Endo.c13.NPR3"), invert = T)
Idents(seurat_obj) <- 'classic1'

features <- c("NDRG1", "PFKFB3", "PDK3",
              'PTPRM', 'PTPN4', 'INPP5D', 'MTMR3','PPP6R3',
              'GSK3B','AXIN1','BTRC','FBXW11','TCF7L2',
              'CD44','ITGA4','CD276','EMILIN2',
              'CD28','SYK', "NFKB1", "IL1R1",
              'MYH6','TNNT2','DSP','MYBPC3','SCN5A',
              'NOTCH3','DLL4','HES1','HEY1','JAG2',
              'NDUFS7', 'SDHB', 'CYCS', 'COX4I1', 'UQCRFS1',
              'ACSL4',"GPX4",'AIFM2','DHODH','SLC40A1', 'FTH1', 'FTL')
features <- features[!duplicated(features)]

types <- c('ACM', 'DCM', 'HCM', 'ICM', 'AMI', 'CHD', 'CS', 'COV' )
DEG <- readRDS("~/PH/results/DEG/files/DEGs_disease_number.rds")
DEG$gene <- DEG$Gene_Symbol
DEG$avg_log2FC <- DEG$logFC
deg_s <- subset(DEG, gene %in% features)

deg_s$types <- ""
deg_s$types[deg_s$gene%in% c("NDRG1", "PFKFB3", "PDK3")] <- "Hypoxia"
deg_s$types[deg_s$gene%in% c('PTPRM', 'PTPN4', 'INPP5D', 'MTMR3','PPP6R3')] <- "Dephosphorylation"
deg_s$types[deg_s$gene%in% c('GSK3B','AXIN1','BTRC','FBXW11','TCF7L2')] <- "WNT signaling pathway"
deg_s$types[deg_s$gene%in% c('MYH6','TNNT2','DSP','MYBPC3','SCN5A')] <- "Cardiac contraction"
deg_s$types[deg_s$gene%in% c('CD44','ITGA4','CD276','EMILIN2')] <- "Cell adhesion"
deg_s$types[deg_s$gene%in% c('CD28','SYK', "NFKB1", "IL1R1")] <- "Immune response"
deg_s$types[deg_s$gene%in% c('NOTCH3','DLL4','HES1','HEY1','JAG2')] <- "NOTCH signaling"
deg_s$types[deg_s$gene%in% c('NDUFS7', 'SDHB', 'CYCS', 'COX4I1', 'UQCRFS1')] <- "Energy metabolism"
deg_s$types[deg_s$gene%in% c('ACSL4',"GPX4",'AIFM2','DHODH','SLC40A1', 'FTH1', 'FTL')] <- "Ferroptosis"

deg_s$types <- factor(deg_s$types, levels = c("Dephosphorylation", "Hypoxia", "WNT signaling pathway", 
                                              'Cardiac contraction',"Cell adhesion", "Ferroptosis",
                                              "Immune response","Energy metabolism", 'NOTCH signaling'))

p <- DotPlot(seurat_obj,features=features[!duplicated(features)], group.by = 'disease')
data <- p[['data']]
colnames(data) <- c("avg.exp", "pct.exp", "gene", "disease", "avg.exp.scaled")
deg_s <-merge(deg_s,data, by=c('gene','disease'))
colnames(deg_s)
deg_s$avg_log2FC[deg_s$avg_log2FC < -2] = -2
deg_s <- deg_s |> filter(abs(avg_log2FC) >= 0.58)|> filter(FDR <= 0.05)

## plot
deg_s %>%
  catdotplot(
    x = factor(gene, level = features),
    y = factor(disease, level = c('ACM', 'DCM', 'HCM', 'ICM', 'AMI', 'HF', 'CHD',  'CS', 'COV' )),
    color = avg_log2FC, size = pct.exp, dot_scale = 5) +
  facet_grid(. ~ types,
             space = "free", scales = "free") +
  scale_color_gradient2(name="log2 FC",
                        breaks=c(-1,0,0.5,1,1.5),
                        low = "navy", mid = "white", high = 'darkred')+
  theme(plot.title = element_text(hjust = 0,vjust = -12, size=8),
        axis.text.x = element_text(face = "italic",size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = "top", panel.spacing.x = unit(0, "pt"))+
  labs(title=' ')
ggsave(file.path(outdir2,'Dotplot_pathway_genes_all.pdf'), p, height=3, width=10.5)

#### 4.Fig2d ----
library(msigdbr)
library(AUCell)

selected_genes <- c('PTPRM','PTPRK','PTPN4','PTPN14','DUSP16','INPP5D','MTMR3','PPP6R3','PPP1R16B','PPP2R5A')
category <- "C5"
Dataset <- msigdbr(species = "Homo sapiens", category = category)
Dataset <- Dataset %>% filter(str_detect(Dataset$gs_name,"GOBP_"))
Dataset$gs_name <- gsub('GOBP_', '', Dataset$gs_name)
Dataset$gs_name2 <- str_to_lower(Dataset$gs_name)
Dataset$gs_name2 <- gsub('_', ' ', Dataset$gs_name2)

feature_x_list <- c(
  'atp metabolic process',
  'atp synthesis coupled electron transport',
  'electron transport chain',
  "mitochondrial respiratory chain complex assembly",
  "aerobic respiration",
  "cellular respiration",
  "oxidative phosphorylation",
  "respiratory electron transport chain",
  "dephosphorylation",                               
  "phosphatidylinositol dephosphorylation")
dataset <- Dataset %>% filter(Dataset$gs_name2 %in% feature_x_list)
geneSets <- lapply(unique(dataset$gs_name2), 
                   function(x){dataset$gene_symbol[dataset$gs_name2 == x]})
names(geneSets) <- unique(dataset$gs_name2)

Idents(sce) <- 'group'
seurat_obj <- subset(sce, idents = 'CVD')
cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@data) 
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
saveRDS(AUCell_socre, file = file.path(outdir, "CVD_VECs_AUCscore.rds"))

##
expr <- t(as.data.frame(seurat_obj[selected_genes,]@assays$RNA@data))
adata <- cbind(expr,AUCell_socre)

method = "pearson"
pathway_cor_test <- data.frame(
  feature_x=NULL,
  feature_y=NULL,
  p_value=NULL,
  estimate=NULL,
  num=NULL,
  method=NULL)

for (feature_x in colnames(expr)) {
  for (feature_y in colnames(AUCell_socre)) {
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

key_gene <- c('PTPRM', 'PTPRK', 'PTPN4', 'PTPN14', 'DUSP16','INPP5D', 'MTMR3', 'PPP6R3', 'PPP1R16B', 'PPP2R5A' )
key_term <- c('atp metabolic process',
              'atp synthesis coupled electron transport',
              'electron transport chain',
              "mitochondrial respiratory chain complex assembly",
              "aerobic respiration",
              "cellular respiration",
              "oxidative phosphorylation",
              "respiratory electron transport chain",
              "dephosphorylation",                               
              "phosphatidylinositol dephosphorylation"                      
)
pathway_cor_test2 <- pathway_cor_test |>
  filter(feature_x %in% key_gene)|>
  filter(feature_y %in% key_term)|>
  filter(p_value <= 0.05)
pathway_cor_test2$feature_x <- factor(pathway_cor_test2$feature_x, levels = key_gene)
pathway_cor_test2$feature_y <- factor(pathway_cor_test2$feature_y, levels = key_term)
pathway_cor_test2$p_value[pathway_cor_test2$p_value < 1.0e-8] <-  1.0e-8

pathway_cor_test2 %>%
  catdotplot(
    x = feature_x,
    y = feature_y,
    size = -log10(p_value),
    color = estimate) +
  coord_fixed() +
  theme(panel.grid = element_line(size = 0.2, color = "lightgrey"),
        axis.text.x = element_text(face = "italic"))

#### 5.FigS5d ----
library(ggrastr)

AUCell_socre <- readRDS(paste0(outdir, "CVD_VECs_AUCscore.rds"))
df <- as.data.frame(AUCell_socre[,c("electron transport chain")])
rownames(df) <- rownames(AUCell_socre)
colnames(df) <- 'ETC'
df2 <- df |> filter(ETC > 0.1)
expr <- AUCell_socre[,c("dephosphorylation", 
                        "electron transport chain")]
colnames(expr) <- gsub(' ', '_', colnames(expr))

method = "pearson"
oknames <- setdiff(rownames(df), rownames(df2))
expr2 <- expr[oknames,]
cor.test.res <- cor.test(expr2[,1],
                         expr2[,2],
                         method = method)
p.value <- cor.test.res$p.value
estimate <- cor.test.res$estimate

expr2 %>%
  ggplot(aes(x = expr2[,1], 
             y = expr2[,2])) +
  ggrastr::rasterise(ggpointdensity::geom_pointdensity(size = 0.5),dpi = 600) +
  geom_smooth(method = "lm",  formula = y ~ x,
              color = "#6b76ae", fill = "#e5a323",
              size = 0.5, alpha = 0.2) +
  theme_bw()+
  theme(aspect.ratio = 1, axis.text = element_text(colour = 'black'),
        legend.margin = margin(l=-8),
        panel.grid = element_blank()) +
  guides(color = guide_colorbar(
    frame.colour = "black", frame.linewidth = 0.5,
    ticks.colour = "black", title = "Density"))+
  scale_color_viridis_c() +
  labs(x = 'dephosphorylation', y = colnames(expr)[i], 
       title = paste0("R = ",round(estimate, 2),',', 'P =', round(p.value, 2)))
ggsave(file.path(outdir2,paste0('low_',colnames(expr)[i],'.pdf')),
       last_plot(),width = 5,height = 5)

#### 6.FigS5e ----
#### 6.0 TF ----
library(Seurat)
library(tidyverse)
library(qs)
library(AUCell)
library(SCENIC)
library(Seurat)
library(pheatmap)
library(patchwork)
library(ggplot2)
library(stringr)
library(IRkernel)
library(tictoc)
library(data.table)

outdir3 <- './TF/files/'
dir.create(outdir3, recursive = T)

seurat_obj0 <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
types <- c('ACM','DCM','HCM','ICM','AMI','CHD','CS','COV')
for(i in types){
  ## UP
  seurat_obj <- seurat_obj0[, seurat_obj0[["disease"]] == i] 
  DEG <- readRDS(paste0('~/PH/results/DEG/files/',i, '_deg_up.rds')) |> pull(Gene_Symbol)
  
  seurat_obj@assays$RNA@counts <- seurat_obj@assays$RNA@counts[DEG, ]
  seurat_obj@assays$RNA@data <- seurat_obj@assays$RNA@data[DEG, ]
  seurat_obj@assays$RNA@meta.features <- seurat_obj@assays$RNA@meta.features[DEG, ]
  length(rownames(seurat_obj))
  saveRDS(seurat_obj, file = paste0(outdir3,i,"_obj_U.rds"))
  
  seurat_obj <- hp_run_pyscenic(x = seurat_obj,
                                species = "human",
                                outdir = paste0(outdir, i, '_U'))
  
  ## DOWN
  seurat_obj <- seurat_obj0[, seurat_obj0[["disease"]] == i]
  DEG <- readRDS(paste0('~/PH/results/DEG/files/',i, '_deg_down.rds')) |> pull(Gene_Symbol)

  seurat_obj@assays$RNA@counts <- seurat_obj@assays$RNA@counts[DEG, ]
  seurat_obj@assays$RNA@data <- seurat_obj@assays$RNA@data[DEG, ]
  seurat_obj@assays$RNA@meta.features <- seurat_obj@assays$RNA@meta.features[DEG, ]
  length(rownames(seurat_obj))
  saveRDS(seurat_obj, file = paste0(outdir3,i,"_obj_D.rds"))
  
  seurat_obj <- hp_run_pyscenic(x = seurat_obj,
                                species = "human",
                                outdir = paste0(outdir, i, '_D'))
}

#### 6.1 heatmap ----
TF_G <- data.frame('symbol'='',
                   'tf'='', 
                   'target_gene'='', 
                   'disease'='', 
                   'change'='')
for(i in types){
  TF_UP <- read.table(paste0('./TF/files/',i,"_U/tfs_targer.tsv"),header=TRUE,sep=",")
  TF_UP$disease <- i
  TF_UP$change <- 'UP'
  length(unique(TF_UP$tf))
  
  TF_DOWN <- read.table(paste0('./TF/files/',i,"_D/tfs_targer.tsv"), header=TRUE,sep=",")
  TF_DOWN$disease <- i
  TF_DOWN$change <- 'DOWN'
  data <- rbind(TF_UP,TF_DOWN)
  TF_G <- rbind(TF_G,data)
}
TF_G <- TF_G[-1,]
saveRDS(TF_G,file.path(outdir,'ALL_TF_G.rds'))

TF_G$Gene <- paste(TF_G$target_gene, TF_G$disease,sep = '_')

tf_s <- c("WT1","SREBF2","SMAD3","NFKB1","TCF4","ERG","ETV6","FOXO1","FOXO3","NFATC2",
          "CEBPB","GATA3","BHLHE41","FOXO4","ERF","HEY1","VEZF1")
TF_p <- TF_G |> filter(TF_G$tf %in% tf_s)
pdata <- TF_p[,c(2,4,5)]
pdata <- unique(pdata)

pivot_table <- pdata %>%
  group_by(disease, tf) %>%
  summarise(change = ifelse(any(change == "UP"), "UP", "DOWN"), .groups = "drop") %>%
  pivot_wider(names_from = tf, values_from = change) %>%
  as.data.frame()
names <- pivot_table$disease

pivot_table <- pivot_table[,-1]
pivot_table[is.na(pivot_table)] <- 0
pivot_table[pivot_table == 'UP'] <- 1
pivot_table[pivot_table == 'DOWN'] <- -1
pivot_table <- as.data.frame(lapply(pivot_table, as.numeric))
rownames(pivot_table) <- names

pivot_table <- pivot_table[, tf_s]
pivot_table <- pivot_table[c('ACM','DCM','HCM','ICM','AMI','CHD','CS','COV'),]
test2 <- as.matrix(pivot_table)

annotation_row = data.frame(factor(disease, levels = disease))
rownames(annotation_row) = disease
rownames(annotation_row) = factor(rownames(annotation_row) ,levels = disease)
colnames(annotation_row) = 'disease'
ann_colors = list(
  disease = c('ACM'="#FB8072", 'DCM'="#80B1D3", 'HCM'="#FDB462", 'ICM'="#8DD3C7",
              'AMI'="#FCCDE5", 'CHD'="#CCEBC5", 'CS'="#377EB8", 'COV'="#1B9E77"))

pdf(file.path(outdir2, 'TF_disease_heatmap.pdf'),width = 6,height = 3)
pheatmap(test2, cluster_rows = F,cluster_cols = F,legend = F,scale='none',
         main='UP/DOWN TFs', show_rownames = F,
         annotation_colors = ann_colors, annotation_row = annotation_row, 
         cellwidth = 12, cellheight = 12, border_color = "grey", fontsize = 8, 
         color = c("#00a6e1","white","#ee6470"))
dev.off()

#### 7.FigS5f ----
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)
library(enrichplot)

features <- c('FOXO1','FOXO3','FOXO4')
for (i in features) {
  genelist <- TF_G[TF_G$tf %in% i, ]$target_gene |> unique()
  bp <-
    enrichGO(
      genelist,
      OrgDb = org.Hs.eg.db,
      keyType = 'SYMBOL',
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2)
  term <- bp@result
  saveRDS(term, paste0(outdir, i, '_terms_GO.rds'))
  write.csv(term, paste0(outdir, i, '_terms_GO.csv'))
}

term1 <- readRDS("./files/FOXO1_terms_GO.rds")
path1 <- c('response to transforming growth factor beta',
           'regulation of autophagy')
term3 <- readRDS("./files/FOXO3_terms_GO.rds")
path3 <- c('regulation of mitotic metaphase/anaphase transition',
           'regulation of G0 to G1 transition')
term4 <- readRDS("./files/FOXO4_terms_GO.rds")
path4 <- c('positive regulation of binding',
           'epithelial cell-cell adhesion')

data <- rbind(term4[term4$Description%in%path4,],
              term3[term3$Description%in%path3,],
              term1[term1$Description%in%path1,])
data$group <- c(rep('FOXO4',2),rep('FOXO3',2),rep('FOXO1',2))
data$Count <- as.numeric(data$Count)
data$pvalue <- as.numeric(data$pvalue)
data$s_pval <- -log10(data$pvalue)
data$labelx=rep(0,nrow(data))
data$labely=seq(1,nrow(data))

ggbarplot(data, x="Description", y="s_pval", fill = "group", color = "white", 
          orientation = "horiz",legend = "right",sort.by.groups=TRUE,
          palette = c('#f1bbba','#92c051','#a3daff'))+
  geom_text(aes(x=labely, y=labelx, label = Description),
            size=3.5, hjust =0)+
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  scale_y_continuous(expand=c(0, 0)) +
  scale_x_discrete(expand=c(0,0))+
  ylab("-log10(pvalue)")+xlab("")+ggtitle("")
ggsave(paste0(outdir2, 'FOXO_GO_barplot.pdf'),
       last_plot(), height=2.5, width=4.5)

#### 8.FigS5g ---
deg_s <- readRDS('~/PH/results/DEG/files/DEGs_disease_number.rds')
genes <- unique(deg_s$Gene_Symbol)
DEG <- readRDS("~/PH/results/DEG/files/DEGs_disease_number.rds")|>filter(change=='Upregulated')

data <- as.data.frame(table(DEG$Gene_Symbol))
count <- data.frame(gene=genes,
                    freq='')
for (i in genes) {
  count[count$gene==i,]$freq <- ifelse(i %in% data$Var1, data[data$Var1 == i,]$Freq, 0)
}
count <- count[order(count$gene),]

df_deg <- data.frame()
for (i in types) {
  deg <- deg_s[deg_s$disease == i,]
  deg <- deg[order(deg$logFC),]
  deg$rank <- 1:nrow(deg)
  df_deg <- rbind(df_deg,deg)
}

df_rank <- data.frame()
for (i in genes) {
  df <- df_deg[df_deg$Gene_Symbol == i,]
  a <- sum(df$rank)/length(genes)
  df_rank <- rbind(df_rank,data.frame(gene=i,
                                      avg_rank=a))
}
df_rank <- df_rank[order(df_rank$gene),]

library(gcookbook)
adata <- cbind(count,df_rank)
adata <- adata[,c(3,2,4)]
cut <- adata[adata$freq>=6 & adata$avg_rank>=6,]$gene 
saveRDS(cut,paste0(outdir, "cutgene_cut6.rds"))

genes <- c('DOCK1','PPARG','ELMO1','CBLB','SPRED2','IGF1R','NR3C2','ZBTB46')
adata$cutoff <- 'none'
adata[adata$gene %in% cut,]$cutoff <- 'label'
adata$label <- ''
adata[adata$gene %in% genes,]$label <- adata[adata$gene %in% genes,]$gene
ggplot(adata,aes(x=freq,y=avg_rank)) +
  geom_point(aes(colour = factor(cutoff)))+
  scale_color_manual(values = c('red','grey'))+
  geom_vline(xintercept = 7,
             lty = 2, col = "black", lwd = 0.5) +
  geom_hline(yintercept = 6,
             lty = 2, col = "black", lwd = 0.5) +
  ggrepel::geom_text_repel(
    aes(label = label),
    # max.overlaps = 1000000,
    max.overlaps = Inf,
    color = "black",
    size = 3,
    fontface = "italic"
  )+
  xlab('Number of disease types') + ylab('Average rank score')+
  scale_y_continuous(limits = c(0,11))+
  theme_classic()+NoLegend()+
  theme(axis.text.x = element_text(colour = 'black'),
        axis.text.y = element_text(colour = 'black'),
        aspect.ratio = 1.2) 
ggsave(file.path(outdir2, 'rank_point_adj_2.pdf'),
       last_plot(), height=5, width=6)

#### 8.FigS5h ----
#### 8.1 inter genes ----
library(readxl)
folder_path <- './database/'
file_list <- list.files(folder_path)

data1 <- read_excel(paste0(folder_path,file_list[1]), sheet = 1)
genelist1 <- data1[,c('disease','gene','scoreGDA')]
genelist1$disease <- "Aortic dissections"
#
data2 <- read_excel(paste0(folder_path,file_list[2]), sheet = 1)
genelist2 <- data2[,c('disease','gene','scoreGDA')]
genelist2$disease[genelist2$disease %in% c('Aneurysms,Aortic','Aortic Aneurysm,Abdominal','Aortic Aneurysm,Thoracic')] <- 
  "Aortic Aneurysm"
genelist2 <- genelist2[genelist2$disease %in% "Aortic Aneurysm",]
#
data3 <- read_excel(paste0(folder_path,file_list[3]), sheet = 1)
genelist3 <- data3[,c('disease','gene','scoreGDA')]
genelist3 <- genelist3[genelist3$disease %in% c('High blood pressure', 'Pulmonary arterial hypertension'),]
#
data4 <- read_excel(paste0(folder_path,file_list[4]), sheet = 1)
genelist4 <- data4[,c('disease','gene','scoreGDA')]
genelist4$disease <- "Atherosclerosis"
#
data5 <- read_excel(paste0(folder_path,file_list[5]), sheet = 1)
genelist5 <- data5[,c('disease','gene','scoreGDA')]
genelist5$disease[genelist5$disease %in% c('Arterioscleroses,Coronary')] <- "Coronary Arterioscleroses"
genelist5$disease[genelist5$disease %in% c('Defect,Congenital Heart')] <- "Congenital Heart Defect"
genelist5$disease[genelist5$disease %in% c('Disease,Ischemic Heart')] <- "Ischemic Heart Disease"
genelist5 <- genelist5[genelist5$disease != 'cardiac toxicity',]
#
data6 <- read_excel(paste0(folder_path,file_list[6]), sheet = 1)
genelist6 <- data6[,c('disease','gene','scoreGDA')]
genelist6$disease[genelist6$disease %in% c('Cardiomyopathies,Dilated', 'Familial dilated cardiomyopathy')] <- 
  "Dilated cardiomyopathy"
genelist6$disease[genelist6$disease %in% c('Cardiomyopathy,Hypertrophic,Familial', 'Hypertrophic cardiomyopathy')] <- 
  "Hypertrophic cardiomyopathy"
#
data7 <- read_excel(paste0(folder_path,file_list[7]), sheet = 1)
genelist7 <- data7[,c('disease','gene','scoreGDA')]
genelist7$disease <- "Coronary Disease"
#
data8 <- read_excel(paste0(folder_path,file_list[8]), sheet = 1)
genelist8 <- data8[,c('disease','gene','scoreGDA')]
genelist8$disease <- "COVID-19"
#
data9 <- read_excel(paste0(folder_path,file_list[9]), sheet = 1)
genelist9 <- data9[,c('disease','gene','scoreGDA')]
genelist9$disease[genelist9$disease %in% 
                    c('Acute myocardial infarction','Infarctions,Myocardial','ST segment elevation myocardial infarction')] <- 
  "Myocardial infarction"
genelist9$disease[genelist9$disease %in% 
                    c('Congestive heart failure','Heart failure','Heart Failure,Diastolic','Heart Failure,Systolic')] <- 
  "Heart failure"
#
data10 <- read_excel(paste0(folder_path,file_list[10]), sheet = 1)
genelist10 <- data10[,c('disease','gene','scoreGDA')]
genelist10$disease <- "Stroke"
#
data11 <- read_excel(paste0(folder_path,file_list[11]), sheet = 1)
genelist11 <- data11[,c('disease','gene','scoreGDA')]
genelist11$disease <- "Thrombosis"
##
genelist <- rbind(genelist1,genelist2,genelist3,genelist4,genelist5,genelist6,
                  genelist7,genelist8,genelist9,genelist10,genelist11)
result <- genelist %>%
  group_by(disease, gene) %>%
  slice_max(scoreGDA, with_ties = FALSE) %>%
  ungroup()

disgenes <- unique(result$gene)
saveRDS(disgenes, paste0(outdir, 'disgenet_disgenes.rds'))

genes_inter_6 <- intersect(disgenes, cut6)
saveRDS(genes_inter_6,paste0(outdir,'cut6_disgene_inter.rds'))

DG_list <- result[result$gene %in% inter_genes_6,]
write.csv(DG_list,"./files/DG_list.csv")

#### 8.2 AUCell ----
ggviolin(data, x='disease',y="AUC",
         xlab = '', ylab = "AUC score",
         color='disease', width = 0.8,
         add = 'boxplot', add.params = list(fill='white', width=0.1))+
  scale_color_manual(values = disease_colors)+
  geom_hline(yintercept = mean(data$AUC),
             linetype = "dashed", size = 0.5, color = "black")+
  theme(axis.text.x = element_text(colour = 'black', size = 10),
        axis.text.y = element_text(colour = 'black', size = 10),
        legend.position = 'none')
ggsave(file.path(outdir2, 'Score_vlnplot_dis.pdf'),
       last_plot(),width = 6, height = 2.5)

##
data_D <- data[data$group == 'CVD',]
sub_data <- data_D %>%
  dplyr::filter(age_group %in% c('40-50','50-60','60-70','>70'))
sub_data$age_group <- factor(sub_data$age_group, levels = c('40-50','50-60','60-70','>70'))
ggviolin(sub_data, x='age_group',y="AUC",
         xlab = '', ylab = "AUC score",
         color='age_group', width = 0.8,
         add = 'boxplot', add.params = list(fill='white', width=0.1))+
  scale_color_manual(values = age_colors)+
  geom_hline(yintercept = mean(data$AUC),
             linetype = "dashed", size = 0.5, color = "black")+
  theme(axis.text.x = element_text(colour = 'black', size = 10),
        axis.text.y = element_text(colour = 'black', size = 10),
        legend.position = 'none')
ggsave(file.path(outdir2, 'Score_vlnplot_age_D_40up.pdf'),
       ggplot2::last_plot(),width = 4, height = 2.5)

##
library(ggsignif)
ggviolin(data_D, x='sex',y="AUC",
         xlab = '', ylab = "AUC score",
         color='sex', width = 0.8,
         add = 'boxplot', add.params = list(fill='white', width=0.1))+
  scale_color_manual(values = sex_colors)+
  geom_hline(yintercept = mean(data[data$sex %in% "Female",]$AUC),
             linetype = "dashed", size = 0.5, color = "black")+
  theme(axis.text.x = element_text(colour = 'black', size = 10),
        axis.text.y = element_text(colour = 'black', size = 10),
        legend.position = 'none') +
  geom_signif(comparisons = list(c("Female", "Male")),
              map_signif_level = TRUE,
              textsize = 4, vjust = 0.5,
              size = 0.5)
ggsave(file.path(outdir2, 'Score_vlnplot_sex_D.pdf'),
       ggplot2::last_plot(),width = 2.5, height = 2.5)

#### 9.FigS5i ----
seurat_obj <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
selected_genes <- genes_inter_6
sub_obj <- subset(seurat_obj0, Article %in% c('Data5_HF_AndrewLKoenig_2022'))
cells_ranking <- AUCell_buildRankings(sub_obj@assays$RNA@counts)
cells_AUC <- AUCell_calcAUC(selected_genes, cells_ranking, aucMaxRank=nrow(cells_ranking)*0.1)
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
sub_obj$AUC  <- AUCell_socre
sub_obj <- subset(sub_obj, donor %in% c('HDCM1','HDCM3','HDCM4','HDCM6','HDCM8'), invert = T)

sub_obj0 <- sub_obj
sub_obj$LVEF <- ''
sub_obj$LVEF[sub_obj$donor %in% c('TWCM-10-5','TWCM-13-17','TWCM-13-208','TWCM-LVAD3')] <- 20
sub_obj$LVEF[sub_obj$donor %in% c('TWCM-11-3')] <- 27
sub_obj$LVEF[sub_obj$donor %in% c('TWCM-11-93')] <- 13
sub_obj$LVEF[sub_obj$donor %in% c('TWCM-13-47')] <- 7
sub_obj$LVEF[sub_obj$donor %in% c('TWCM-13-84')] <- 12.5
sub_obj$LVEF[sub_obj$donor %in% c('TWCM-13-102')] <- 25
sub_obj$LVEF[sub_obj$donor %in% c('TWCM-13-181')] <- 28.1
sub_obj$LVEF[sub_obj$donor %in% c('H_ZC-LVAD')] <- 30
sub_obj$LVEF[sub_obj$donor %in% c('TWCM-13-285','TWCM-LVAD2')] <- 17.5
sub_obj$LVEF <- as.numeric(sub_obj$LVEF)

meta <- sub_obj@meta.data
data <- meta[,c('donor','sample','cell_type','age_group','classic1','disease','AUC','LVEF')]
adata <-
  aggregate(data[, c('AUC','LVEF')], list(data$sample), FUN = mean)
rownames(adata) <- adata$Group.1
adata <-adata[,-1]

correlation <- cor.test(adata$LVEF, adata$AUC, method = "pearson", alternative = "two.sided")
cor_label <- paste0("Correlation: ", round(correlation$estimate, 3), "\n", "p-value: ", format.pval(correlation$p.value, digits = 3))
dev.off()
dev.new()
ggplot(adata, aes(x = AUC, y = LVEF)) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.4) +
  geom_point(color = 'grey') +
  theme_classic() +
  labs(title = "",x = 'AUC score',y = "LVEF(%)") +
  annotate("text", x = min(adata$AUC, na.rm = TRUE), y = max(adata$LVEF, na.rm = TRUE),
           label = cor_label, hjust = 0, vjust = 1, size = 4)
ggsave(file.path(outdir2, paste0('5_AUC_LVEF_cor.pdf')),
       last_plot(),height=3.5, width=3.5)

#### 10.Fig2f ----
gene_s <- readRDS(paste0(outdir, "cutgene_cut6.rds"))
selected_genes <- c("MTOR","IKBKB","PPARG","KALRN","NR3C2","SMAD3","DAB2IP",
                    "AKT3","RASA1","TCF7L1","INPP5D","ATXN2","BAZ2B","FLT1")
##
seurat_obj0 <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
Idents(seurat_obj0) <- 'classic1'
seurat_obj <- subset(seurat_obj0, idents = c("LECs","EndoECs"), invert = T)
seurat_obj$classic1 <- as.character(seurat_obj$classic1)
table(seurat_obj$classic1)
seurat_obj$classic1 <- factor(seurat_obj$classic1,levels = c("ArtECs","CapECs","VenECs"))
datalist <- list(exp=seurat_obj@assays[["RNA"]]@counts,
                 annot=seurat_obj@meta.data)

annotLevels = list(level1class=datalist$annot$disease,level2class=datalist$annot$classic1)
ctd = generate_celltype_data(exp=datalist[[1]],annotLevels=annotLevels,groupName="filter")
print(ctd)
load(ctd)
human.hits <- gene_s
human.bg = unique(rownames(seurat_obj@assays[["RNA"]]@counts))
reps=14819
level=2

full_results = bootstrap_enrichment_test(sct_data=ctd, hits=human.hits, bg=human.bg,
                                         reps=reps, annotLevel=level,
                                         genelistSpecies = "human",
                                         sctSpecies = "human",
                                         sctSpecies_origin = "human",
                                         output_species = "human",)
saveRDS(full_results, paste0(outdir, 'EWCE_res_endoonly_classic_VEC.rds'))
full_results$results$CellType <- factor(full_results$results$CellType,levels = c("ArtECs","CapECs","VenECs"))
pdf(paste0(outdir2,'ewce_bar_endo_classic1_VEC.pdf'), width = 2.8, height = 4)
ewce_barplot(full_results$results,mtc_method="BH")
dev.off()

meta <- full_results[["gene_data"]] |> as.data.frame()
meta <- meta[,c('gene','CellType','p')]
pivot_df <- meta %>%
  pivot_wider(names_from = CellType, values_from = p, values_fill = list(p = 0)) %>%
  as.data.frame()
rownames(pivot_df) <- pivot_df$gene
pivot_df <- pivot_df[,-1]
pivot_df <- pivot_df[selected_genes,]
classic1 <- c("ArtECs","CapECs","VenECs")
names <- colnames(pivot_df)
colnames(pivot_df) <- names
pivot_df <- pivot_df[,classic1]
new_df <- pivot_df
new_df[new_df == 1] <- " "
new_df[new_df == 0] <- "*"

df_endo <- new_df

##
seurat_obj0 <- qread("~/PH/results/named/files/Pan_all.qs")
seurat_obj <- subset(seurat_obj0, idents = c("LEC.c12.FLT4","Endo.c13.NPR3"), invert = T)
seurat_obj$type <- seurat_obj$cell_type
celltypes <- c("Art.c01.DKK2","Art.c02.NEBL","Art.c03.VEGFA","Cap.c04.PAPSS2","Cap.c05.RGCC","Cap.c06.RDH10",
               "Cap.c07.CA2","Cap.c08.TMEM163","Cap.c09.RGS6","Ven.c10.VCAM1","Ven.c11.ACTA2")
seurat_obj$type[seurat_obj$type %in% celltypes] <- 'Endothelial'

datalist <- list(exp=seurat_obj@assays[["RNA"]]@counts,
                 annot=seurat_obj@meta.data)
annotLevels = list(level1class=datalist$annot$cell_type,level2class=datalist$annot$type)
ctd = generate_celltype_data(exp=datalist[[1]],annotLevels=annotLevels,groupName="filter")
print(ctd)
load(ctd)
human.hits <- gene_s
human.bg = unique(rownames(seurat_obj@assays[["RNA"]]@counts))
reps=14819
level=2

full_results = bootstrap_enrichment_test(sct_data=ctd, hits=human.hits, bg=human.bg,
                                         reps=reps, annotLevel=level,
                                         genelistSpecies = "human",
                                         sctSpecies = "human",
                                         sctSpecies_origin = "human",
                                         output_species = "human",)
saveRDS(full_results, paste0(outdir, 'EWCE_res_alltypes_classic1_VEC.rds'))

alltypes <- c('Cardiomyocyte' ,'Endothelial' ,'Fibroblast' ,'Lymphoid' ,'Myeloid')
full_results$results$CellType <- factor(full_results$results$CellType, 
                                        levels = alltypes)
pdf(paste0(outdir2,'ewce_bar_alltypes_classic1_VEC.pdf'), width = 3.5, height = 4)
ewce_barplot(full_results$results,mtc_method="BH")
dev.off()

meta <- full_results[["gene_data"]] |> as.data.frame()
meta <- meta[,c('gene','CellType','p')]
pivot_df <- meta %>%
  pivot_wider(names_from = CellType, values_from = p, values_fill = list(p = 0)) %>%
  as.data.frame()
rownames(pivot_df) <- pivot_df$gene
pivot_df <- pivot_df[,-1]
pivot_df <- pivot_df[selected_genes,]

new_df <- pivot_df
new_df[new_df == 1] <- " "
new_df[new_df == 0] <- "*"
df_type <- new_df[,c('Endothelial', 'Cardiomyocyte', 'Fibroblast', 'Lymphoid', 'Myeloid')]

##
df_all <- cbind(df_endo,df_type)
df_all <- df_all[,-4]
colnames(df_all)[1:3] <- c('ArtECs','CapECs','VenECs')

## heatmap
seurat_obj$cell_type[seurat_obj$cell_type %in% c("Art.c01.DKK2","Art.c02.NEBL","Art.c03.VEGFA")] <- 'ArtECs'
seurat_obj$cell_type[seurat_obj$cell_type %in% c("Cap.c04.PAPSS2","Cap.c05.RGCC","Cap.c06.RDH10",
                                                 "Cap.c07.CA2","Cap.c08.TMEM163","Cap.c09.RGS6")] <- 'CapECs'
seurat_obj$cell_type[seurat_obj$cell_type %in% c("Ven.c10.VCAM1","Ven.c11.ACTA2")] <- 'VenECs'
adata <- t(as.data.frame(seurat_obj@assays$RNA@data[selected_genes,]))
adata <- as.data.frame(adata)
adata$cell_type <- seurat_obj$cell_type

pdata <- aggregate(adata[, -length(colnames(adata))], list(adata$cell_type), FUN = mean)
rownames(pdata) <- pdata[,1]
pdata <- pdata[,-1]
pdata <- pdata[c('ArtECs','CapECs','VenECs',
                 'Cardiomyocyte', 'Fibroblast', 'Lymphoid', 'Myeloid'),]
pdata_f <- t(pdata)
data_scale <- round(t(apply(pdata_f, 1, scale)),2)
colnames(data_scale) <- colnames(pdata_f)

color_palette <- c(colorRampPalette(colors = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0"))(50),
                   colorRampPalette(colors = c("#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(50))
p <- pheatmap(data_scale,cluster_cols = F, 
              cluster_rows = F, scale="none",
              gaps_col = 3, gaps_row = 7,
              display_numbers = df_all, fontsize_number = 10,
              border="white",
              cellwidth = 12,cellheight =12,
              color = color_palette,
              fontsize_row = 10,fontsize_col = 10)
ggsave(file.path(outdir2, 'ewce_heatmap_classic1_VEC.pdf'),
       p,limitsize = FALSE, height=5, width=5)

#### 11.Fig2g ----
category <- read.csv('drug-target_dicts.csv')
s_gene <- c('ATXN2','HTT','ADRA1D','DAB2IP','BAZ2B','PPARD','PPARG','NR3C2','SLC22A1',
            'SMAD3','NFKB1','MTOR','MAPK1','FLT1','CACNB4','HDAC9','SMAD3','GABBR1')
category2 <- category |> filter(category$Gene %in% s_gene)
s_drug <- c('CHEMBL772|RESERPINE',
            'CHEMBL1484|NICARDIPINE',
            'CHEMBL254219|DIGITOXIN', 
            'CHEMBL1463345|CANRENONE', 
            'PIOGLITAZONE HYDROCHLORIDE',
            'CHEMBL1567|SUNITINIB MALATE',
            'TIVOZANIB HYDROCHLORIDE',
            'CHEMBL843|ROSIGLITAZONE MALEATE',
            'LENVATINIB MESYLATE' )
category3 <- category2 |> filter(category2$Drug %in% s_drug)
category3$D_G_names <- paste(category3$Drug, category3$Gene,sep = '|')
unique(category3$D_G_names)

##
library(AUCell)
seurat_obj <- qread('~/PH/results/files/Pan_all_CVD_D2.qs')
seurat_obj@meta.data <- seurat_obj@meta.data[ ,c("ID", "orig.ident", "disease", "cell_type")]

category <- read.csv('./files/drug-target_dicts.csv')
dataset <- category[,1:2]
colnames(dataset) <- c('gene_symbol', 'gs_name')
dataset <- unique(dataset)
geneSets <- lapply(unique(dataset$gs_name), 
                   function(x){dataset$gene_symbol[dataset$gs_name == x]})
names(geneSets) <- unique(dataset$gs_name)
cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@data) 
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
colnames(AUCell_socre) <- colnames(AUCell_socre)

seurat_obj$ID <- rownames(seurat_obj@meta.data)
meta <- AUCell_socre
meta <- as.data.frame(meta)
meta$cell_type <- seurat_obj$cell_type
object <- aggregate(meta[, -length(colnames(meta))], list(meta$cell_type), FUN = mean)
rownames(object) <- object$Group.1
object <-object[,-1]
saveRDS(object,file.path(outdir, "AUCscore_celltypes.rds"))

object2 <- object[, unique(category3$Drug)]
colnames(object2) <- unique(category3$D_G_names)[-1]
object2$`CHEMBL1484|NICARDIPINE|ATXN2` <- object2$`CHEMBL1484|NICARDIPINE|MTOR`
object2 <- object2[,c(7,1:6)]
object3 <- as.data.frame(t(object2))
object3 <- data[, c(12:15)]
object3$ArtECs <- rowMeans(data[, 1:3])
object3$CapECs <- rowMeans(data[, 4:9])
object3$VenECs <- rowMeans(data[, 10:11])
saveRDS(object3,file.path(outdir, "AUCscore_classic_sub.rds"))

object3 <- object3[,c(5,6,7,1,2,3,4)]
object3[object3 > 0.3] = 0.3
p1 <- pheatmap(object3, treeheight_row = 5,
               cluster_cols = F, cluster_rows = F, scale="none",
               fontsize_number = 20, border="grey", 
               shape = "circle",cellwidth = 10,cellheight =10,
               gaps_col = c(3),
               col = c(colorRampPalette(colors = c("white","#FFFFF3","#FDDBC7","#F4A582","#D6604D","#B2182B"))(50)),
               fontsize_row = 10,fontsize_col = 12 
) 
ggsave(file.path(outdir2, 'AUCscore_classic.pdf'),
       p1,limitsize = FALSE, height=2.5, width=5)
