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

#### 0.deg ----
seurat_obj0 <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
diseases <- c("ACM","DCM","HCM","ICM","AMI","CHD","CS","COV")
for(i in diseases){
  seurat_obj <- seurat_obj0[ , seurat_obj0$disease %in% c(i, "Normal")]
  Idents(seurat_obj) <- 'disease'

  outdir3 <- paste0('~/PH/results/DEG/', i, '/')
  dir.create(outdir3)
  
  ## DEG
  markers <- FindMarkers(seurat_obj,
                         ident.1 = i,
                         ident.2 = "Normal",
                         min.pct = 0,
                         logfc.threshold = 0)
  markers$gene <- rownames(markers)
  
  logFC.cutoff <- 0.25
  pvalue.cutoff <- 0.05
  
  markers$change <-
    as.factor(ifelse(
      markers$p_val_adj < pvalue.cutoff &
        abs(markers$avg_log2FC) > logFC.cutoff,
      ifelse(markers$avg_log2FC > logFC.cutoff , 'UP', 'DOWN'),
      'NOT'
    ))
  table(markers$change)
  
  ## save
  saveRDS(markers, file.path(outdir3, paste0(i,"_deg.rds")))
  write.table(markers, file.path(outdir3, paste0(i,"_gsea_input.txt")), quote = F, sep = ",", row.names = F)
}

DEG <- data.frame()
for (i in diseases) {
  path <- paste0('~/PH/results/DEG/',i,'/',i,'_deg.rds')
  deg <- readRDS(path)
  deg$cluster <- i
  DEG <- rbind(DEG,deg)
}
saveRDS(DEG, '~/PH/results/DEG/DEG_diseases.rds')

#### 1.Fig2a ----
## UP
DEG <- readRDS('~/PH/results/DEG/DEG_diseases.rds') |> filter(avg_log2FC > 0.25 & p_val_adj < 0.05)
DEG_split <- split(DEG, DEG$cluster)
for(i in 1:length(table(DEG$cluster))){
  DEG_split[[i]] <- DEG_split[[i]]$gene
}
upsetData=fromList(DEG_split)

color_ct <- c("#FCCDE5", "#377EB8", "#8DD3C7", "#FDB462",
              "#FB8072", "#80B1D3", "#1B9E77", "#CCEBC5")
pdf(file=file.path(outdir2, paste0(case, 'UP_upset.pdf')),onefile = FALSE,width=6,height=5)
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
        list(query=intersects, params=list("ACM","DCM","HCM","ICM", "AMI","CHD","CS","COV"), color= "#c82d31", active=T),
        list(query=intersects, params=list("ACM","DCM","HCM","ICM", "AMI","CS","COV"), color= "#c82d31", active=T),
        list(query=intersects, params=list("ACM","DCM","HCM","ICM", "AMI","CS","COV"), color= "#c82d31", active=T),
        list(query=intersects, params=list("AMI"), color="#eeb401", active=T),
        list(query=intersects, params=list('HCM'), color="#eeb401", active=T),
        list(query=intersects, params=list('CS'), color="#eeb401", active=T),
        list(query=intersects, params=list('COV'), color="#eeb401", active=T)
      ))  
dev.off()

## DOWN
DEG <- readRDS('~/PH/results/DEG/DEG_diseases.rds') |> filter(avg_log2FC < -0.25 & p_val_adj < 0.05)
DEG_split <- split(DEG, DEG$cluster)
for(i in 1:length(table(DEG$cluster))){
  DEG_split[[i]] <- DEG_split[[i]]$gene
}
upsetData=fromList(DEG_split)

color_ct <- c("#CCEBC5", "#1B9E77", "#FDB462", "#377EB8",
              "#8DD3C7", "#FCCDE5", "#FB8072", "#80B1D3")
pdf(file=file.path(outdir2, paste0(case, 'DOWN_upset.pdf')),onefile = FALSE,width=6,height=5)
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
        list(query=intersects, params=list("ACM","DCM","HCM","ICM", "AMI","CHD","CS","COV"), color= "#0e2c82", active=T),
        list(query=intersects, params=list("HCM","ICM", "AMI","CHD","CS","COV"), color= "#0e2c82", active=T),
        list(query=intersects, params=list("ACM","DCM","HCM","ICM","CHD","CS","COV"), color= "#0e2c82", active=T),
        list(query=intersects, params=list('CHD'), color="#eeb401", active=T)
      ))
dev.off()

## Venn
input_dir <- '~/PH/results/DEG/'

diseases <- c('ACM','DCM','HCM','ICM','AMI','CHD','CS','COV')
DEG <- data.frame()
for (i in diseases) {
  dir <- paste0(input_dir,i,'/',i,'_deg.rds')
  deg <- readRDS(dir)
  deg2 <- deg[,c('gene', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj', 'change')]
  deg2$cluster <- i
  DEG <- rbind(DEG,deg2)
}
DEG$cluster <- factor(DEG$cluster,levels = diseases)
## UP 
library(venn)
library(VennDiagram)

DEG_up <- DEG |> filter(DEG$change == 'UP')
DEG_split <- split(DEG_up, DEG_up$cluster)
for(i in 1:length(table(DEG_up$cluster))){
  DEG_split[[i]] <- DEG_split[[i]]$gene
}
upsetData=fromList(DEG_split)
##
pdf(paste0(outdir2,"up_venn.pdf"), height = 4, width = 5)
venn(DEG_split[1:4],
     zcolor = disease_colors[1:4],
     opacity = 0.8,
     box = F,
     ilcs = 0.8,
     sncs = 1
)
dev.off()

## DOWN
DEG_down <- DEG |> filter(DEG$change == 'DOWN')
DEG_split <- split(DEG_down, DEG_down$cluster)
for(i in 1:length(table(DEG_down$cluster))){
  DEG_split[[i]] <- DEG_split[[i]]$gene
}
upsetData=fromList(DEG_split)
##
pdf(paste0(outdir2,"down_venn.pdf"), height = 4, width = 5)
venn(DEG_split[1:4],
     zcolor = disease_colors[1:4],
     opacity = 0.8,
     box = F,
     ilcs = 0.8,
     sncs = 1
)
dev.off()

#### 2.FS5b-c ----
library(clusterProfiler)
library(org.Hs.eg.db)
library(patchwork)
library(enrichplot)
#### 2.1 GO ----
## UP
DEG <- readRDS("~/PH/results/DEG/DEG_diseases.rds") |> filter(change == "UP") 
DEG_split <- split(DEG, DEG$cluster)
for(i in 1:length(table(DEG$cluster))){
  DEG_split[[i]] <- DEG_split[[i]]$gene
}

diseases <- c("ACM","DCM","HCM","ICM","AMI","CHD","CS","COV")
pathway_up <- data.frame()
for (i in 1:8) {
  type <- diseases[i]
  genelist <- DEG_split[[type]]
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
  term$Cluster <- type
  
  pathway_up <- rbind(pathway_up, term)
}
pathway_up <- pathway_up[pathway_up$pvalue < 0.05,]
saveRDS(pathway_up, './files/pathway_up.rds')

## DOWN
DEG <- readRDS("~/PH/results/DEG/DEG_diseases.rds") |> filter(change == "DOWN") 
DEG_split <- split(DEG, DEG$cluster)
for(i in 1:length(table(DEG$cluster))){
  DEG_split[[i]] <- DEG_split[[i]]$gene
}

diseases <- c("ACM","DCM","HCM","ICM","AMI","CHD","CS","COV")
pathway_down <- data.frame()
for (i in 1:8) {
  type <- diseases[i]
  genelist <- DEG_split[[type]]
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
  term$Cluster <- type
  
  pathway_down <- rbind(pathway_down, term)
}
pathway_down <- pathway_down[pathway_down$pvalue < 0.05,]
saveRDS(pathway_down, './files/pathway_down.rds')

#### 2.2 plot ----
DEG <- readRDS("~/PH/results/DEG/DEG_diseases.rds") |> filter(avg_log2FC >= 0.25 & p_val < 0.05) 
DEG_split <- split(DEG, DEG$cluster)
for(i in 1:length(table(DEG$cluster))){
  DEG_split[[i]] <- DEG_split[[i]]$gene
}

data <- lapply(X = DEG_split, FUN = function(x) {
  x <- bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  x <- x[,-1]
})

GO_pathway <- compareCluster(data, 
                             fun='enrichGO', ont= 'BP',
                             OrgDb='org.Hs.eg.db' ,
                             pvalueCutoff = 0.05, qvalueCutoff = 1)
GO_results <- pairwise_termsim(GO_pathway, showCategory = dim(GO_pathway)[1])
pathway_my <- readRDS('./files/pathway_up.rds')
GO_pathway@compareClusterResult <- pathway_my
GO_results <- pairwise_termsim(GO_pathway, showCategory = dim(GO_pathway)[1])
GO_results@compareClusterResult$numb <- ave(GO_results@compareClusterResult$Description, 
                                            GO_results@compareClusterResult$Description, 
                                            FUN = length)
saveRDS(GO_results, file.path(outdir, 'GO_results_UP_my.rds'))

# plot
GO_results <- readRDS('./files/GO_results_UP_my.rds')
terms <- c(
  'glycoprotein biosynthetic process',
  'glycoprotein metabolic process',
  'glycosylation',
  'protein glycosylation',
  'protein O-linked glycosylation',
  'demethylation',
  'DNA dealkylation',
  'venous blood vessel development',
  'vascular transport',
  'transport across blood-brain barrier',
  'ventricular cardiac muscle cell action potential',
  'vascular process in circulatory system',
  'regulation of vascular associated smooth muscle contraction',
  'regulation of vascular permeability',
  'vasculogenesis',
  'response to BMP',
  'Wnt signaling pathway',
  'cell-cell signaling by wnt',
  'canonical Wnt signaling pathway',
  'regulation of heart rate',
  'regulation of heart rate by cardiac conduction',
  'heart contraction',
  'G protein-coupled receptor signaling pathway involved in heart process',
  'endocardium development',
  'cellular response to fluid shear stress',
  'cellular response to vitamin',
  'cellular response to vitamin D',
  'regulation of vitamin D receptor signaling pathway',
  'response to decreased oxygen levels',
  'response to oxygen levels',
  'response to hypoxia',
  'positive regulation of cytokine-mediated signaling pathway',
  'regulation of lipoprotein metabolic process',
  'regulation of chemokine-mediated signaling pathway',
  'inflammasome complex assembly',
  'cellular response to interleukin-6',
  'interleukin-6-mediated signaling pathway',
  'interleukin-1-mediated signaling pathway',
  'positive regulation of interleukin-2 production',
  'cell-matrix adhesion',
  'cell-substrate junction organization',
  'peptidyl-serine phosphorylation',
  'positive regulation of T cell migration',
  'regulation of dendritic cell apoptotic process',
  'regulation of lymphocyte chemotaxis',
  'B cell receptor signaling pathway',
  'positive regulation of lymphocyte chemotaxis',
  'small GTPase mediated signal transduction',
  'Ras protein signal transduction',
  'protein methylation',
  'histone modification',
  'protein acylation',
  'peptidyl-serine phosphorylation',
  'peptidyl-serine modification',
  'immune response-activating signaling pathway',
  'positive regulation of blood vessel endothelial cell migration',
  'activation of innate immune response',
  'immune response-activating cell surface receptor signaling pathway',
  'activation of immune response',
  'regulation of I-kappaB kinase/NF-kappaB signaling',
  'I-kappaB kinase/NF-kappaB signaling',
  'positive regulation of blood vessel endothelial cell migration',
  'blood vessel endothelial cell migration',
  'endothelial cell migration'
)
terms[!(terms %in% GO_results@compareClusterResult$Description)]
pdf(file = paste0(outdir2, "DEG_FUN_UP.pdf"), width = 5,height = 5)
emapplot(GO_results, showCategory = terms,
         layout = "kk",
         group_category =T,shadowtext = T,
         cex_category = 0.5, cex_line = 0.3, cex_label_category = 1,
         node_label = "category", repel = T,
         cex_pie2axis = 3, legend_n = 5) + 
  scale_fill_manual(values = disease_colors[2:9])
dev.off()

## DOWN
DEG <- readRDS("~/PH/results/DEG/DEG_diseases.rds") |> filter(avg_log2FC <= -0.25 & p_val < 0.05) 
DEG_split <- split(DEG, DEG$cluster)
for(i in 1:length(table(DEG$cluster))){
  DEG_split[[i]] <- DEG_split[[i]]$gene
}

data <- lapply(X = DEG_split, FUN = function(x) {
  x <- bitr(x, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  x <- x[,-1]
})
 
GO_pathway <- compareCluster(data, 
                             fun='enrichGO', ont= 'BP',
                             OrgDb='org.Hs.eg.db' ,
                             pvalueCutoff = 0.05, qvalueCutoff = 1)
GO_results <- pairwise_termsim(GO_pathway, showCategory = dim(GO_pathway)[1])

pathway_my <- readRDS('./files/pathway_down.rds')
GO_pathway@compareClusterResult <- pathway_my
GO_results <- pairwise_termsim(GO_pathway, showCategory = dim(GO_pathway)[1])
GO_results@compareClusterResult$numb <- ave(GO_results@compareClusterResult$Description, 
                                            GO_results@compareClusterResult$Description, 
                                            FUN = length)
saveRDS(GO_results, file.path(outdir, 'GO_results_DOWN_my.rds'))

# plot
GO_results <- readRDS(file.path(outdir, 'GO_results_DOWN_my.rds'))
terms <- c(
  'regulation of protein localization to Cajal body',
  'positive regulation of protein localization to Cajal body',
  'protein localization to Cajal body',
  'RNA localization to Cajal body',
  'telomerase RNA localization to Cajal body',
  'cell migration involved in sprouting angiogenesis',
  'regulation of angiogenesis',
  'positive regulation of angiogenesis',
  'regulation of vasculogenesis',
  'vascular associated smooth muscle cell development',
  'blood vessel endothelial cell migration',
  'positive regulation of vasculature development',
  'endothelial cell migration',
  'blood vessel endothelial cell migration',
  'regulation of endothelial cell migration',
  'regulation of DNA repair',
  'positive regulation of response to DNA damage stimulus',
  'regulation of response to DNA damage stimulus',
  'telomere maintenance in response to DNA damage',
  'histone H3 deacetylation',
  'histone H2A acetylation',
  'protein acetylation',
  'positive regulation of DNA repair',
  'stem cell population maintenance',
  'maintenance of cell polarity',
  'maintenance of apical/basal cell polarity',
  'stem cell population maintenance',
  'regulation of hematopoietic progenitor cell differentiation',
  'immunoglobulin mediated immune response',
  'regulation of type 2 immune response',
  'positive regulation of innate immune response',
  'activation of innate immune response',
  'regulation of innate immune response',
  'myeloid cell differentiation',
  'cellular response to interleukin-1',
  'response to interleukin-1',
  'interferon-beta production',
  'interleukin-1 production',
  'regulation of interleukin-1 production',
  'regulation of interferon-beta production',
  'negative regulation of NIK/NF-kappaB signaling',
  'T cell activation involved in immune response',
  'inflammatory cell apoptotic process',
  'aorta morphogenesis',
  'interleukin-10 production',
  'regulation of interleukin-10 production',
  'myeloid cell differentiation',
  'activation of innate immune response',
  'B cell mediated immunity',
  'immunological memory process',
  'regulation of innate immune response',
  'T cell activation involved in immune response',
  'inflammatory cell apoptotic process',
  'myeloid cell homeostasis',
  'antigen processing and presentation of endogenous peptide antigen via MHC class I',
  'antigen processing and presentation of peptide antigen via MHC class I',
  'antigen processing and presentation of peptide antigen via MHC class II',
  'antigen processing and presentation of peptide or polysaccharide antigen via MHC class II',
  'MHC class II protein complex assembly',
  'MHC protein complex assembly',
  'peptide antigen assembly with MHC class II protein complex',
  'T cell activation via T cell receptor contact with antigen bound to MHC molecule on antigen presenting cell',
  'blood vessel endothelial cell migration',
  'endothelial cell migration',
  'regulation of endothelial cell migration',
  'nucleoside triphosphate biosynthetic process',
  'purine ribonucleoside triphosphate metabolic process',
  'release of cytochrome c from mitochondria',
  'aerobic respiration',
  'electron transport chain',
  'oxidative phosphorylation',
  'mitochondrial ATP synthesis coupled electron transport',
  'NADH dehydrogenase complex assembly',
  'ATP synthesis coupled electron transport',
  'ATP metabolic process',
  'cytoplasmic translation'
)
terms[!(terms %in% GO_results@compareClusterResult$Description)]
pdf(file = paste0(outdir2, "DEG_FUN_DOWN.pdf"), width = 5,height = 5)
emapplot(GO_results, showCategory = terms,
         layout = "kk", # kk, mds, graphopt, fr
         group_category = T,shadowtext = T, # 是否分组聚类
         cex_category = 0.5, cex_line = 0.3, cex_label_category = 1, #圆圈线标签
         node_label = "category", repel = T, #nCluster =8,
         cex_pie2axis = 3, legend_n = 5) + 
  scale_fill_manual(values = disease_colors[2:9])
dev.off()

#### 3.Fig2b ----
## UP
GO_results <- readRDS(file.path(outdir, 'GO_results_UP_my.rds'))
GO_results <- readRDS('./files/GO_results_UP_my.rds')
terms <- c(
  'canonical Wnt signaling pathway',
  'activation of immune response',
  'I-kappaB kinase/NF-kappaB signaling',
  'interleukin-1-mediated signaling pathway',
  'inflammasome complex assembly',
  'regulation of leukocyte degranulation',
  'vascular process in circulatory system',
  'cellular response to vitamin D',
  'embryonic heart tube development',
  'cellular response to interleukin-6',
  'vascular transport',
  'interleukin-6-mediated signaling pathway',
  'AV node cell to bundle of His cell communication',
  'cell-cell signaling involved in cardiac conduction',
  'response to decreased oxygen levels',
  'response to hypoxia'
)
terms[!(terms %in% GO_results@compareClusterResult$Description)]

meta <- GO_results@compareClusterResult[GO_results@compareClusterResult$Description %in% terms,]
meta$logP <- -log10(meta$pvalue)
meta$Cluster <- factor(meta$Cluster, levels = c('Normal','ACM','DCM','HCM','ICM','AMI','CHD','CS','COV'))
meta$Description <- factor(meta$Description, levels = terms)
meta$Count[meta$Count < 20] <- 20
meta_up <- meta

## plot
ggplot(meta, aes(x = Description, y = Cluster)) +
  geom_point(aes(color = logP, size = Count), shape = 15) +
  scale_size_continuous(range = c(3,5))+
  scale_color_gradient(low = "#F7CE68", high = "#FF2525") +
  labs(color = "-log10(pvalue)", x = '', y = '') +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 10, angle = 30, hjust = 1),
        axis.text.y = element_text(color = "black", size = 10),
        plot.margin = margin(l = 80, t = 30))
ggsave(paste0(outdir2, 'heatmap_up_s.pdf'), last_plot(), height = 4.5, width = 12)

## DOWN
GO_results <- readRDS(file.path(outdir, 'GO_results_DOWN_my.rds'))
GO_results <- readRDS('./files/GO_results_DOWN_my.rds')
terms <- c(
  'ATP metabolic process',
  'antigen processing and presentation of peptide antigen via MHC class II',
  'antigen processing and presentation of peptide antigen via MHC class I',
  'glucose metabolic process',
  'fatty acid metabolic process',
  'tricarboxylic acid metabolic process',
  'acyl-CoA biosynthetic process',
  'purine ribonucleoside metabolic process',
  'regulation of leukocyte migration',
  'cell migration involved in sprouting angiogenesis',
  'regulation of monocyte chemotaxis',
  'cellular response to interleukin-1',
  'interferon-beta production',
  'regulation of DNA repair',
  'stress granule assembly',
  'maintenance of cell polarity',
  'positive regulation of epithelial to mesenchymal transition'
)

terms[!(terms %in% GO_results@compareClusterResult$Description)]

meta <- GO_results@compareClusterResult[GO_results@compareClusterResult$Description %in% terms,]
meta$logP <- -log10(meta$pvalue)
meta$Cluster <- factor(meta$Cluster, levels = c('Normal','ACM','DCM','HCM','ICM','AMI','CHD','CS','COV'))
meta$Description <- factor(meta$Description, levels = terms)
meta_down <- meta
meta_down$logP <- -meta_down$logP
meta_down$logP[meta_down$logP < -8] <- -8

meta_all <- rbind(meta_up, meta_down)
ggplot(meta_all, aes(x = Description, y = Cluster)) +
  geom_point(aes(color = logP, size = Count), shape = 15) +
  scale_size_continuous(range = c(2,5))+
  scale_color_gradientn(colors = c("#087DEA","#2c8fec","#74b4f0","#a7cef3","#F7CE68","orange","#FF2525"))+
  labs(color = "-log10(pvalue)", x = '', y = '') +
  theme_bw() +
  theme(axis.text.x = element_text(color = "black", size = 10, angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black", size = 10),
        plot.margin = margin(l = 120, t = 30))
ggsave(paste0(outdir2, 'heatmap_combined.pdf'), last_plot(), height = 6, width = 10)

#### 4.FigS5d ----
seurat_obj <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")

features <- c("HIF1A", "NDRG1", "PFKFB3", "PDK3", "SERPINE1", "EGFR", "BCL2",#'Hypoxia',
              "MGAT5", "B4GALT1", "B4GALT5", "ST6GAL1", # 'Glycosylation'
              "VEGFA", "JAG1", "JAG2", "THBD","TIMP1", "FSTL1",'PDGFRB','KDR','TGFBR2','GJC1','TBX20','TGFBR3',"PIK3R3", "HRAS",# "Angiogenesis"
              "NFKB1", "IL1R1", "IL18R1", "TNFAIP8", "STAT3", "IFNGR2", 'CD40','CD36','ANXA1','HLA-E', 'HLA-DRB1', #'Inflammatory response',
              'APH1A','DLL4','GATA2','HES1','HEY1','JAG2','PDCD10',# NOTCH
              'FUS', 'POLDIP2', 'CCNG1', 'PSMD14', 'ERCC1', 'UBC',# DNA damage
              "NDUFS2","NDUFS3",  "COX4I1", "UQCRC1", "ATP6AP1", "SDHA", "CYCS",# 'Oxidative phosphorylation'
              "GPX4", "SLC40A1", "FTH1", "FTL" #'Ferroptosis
)
features <- features[!duplicated(features)]

types <- c('ACM', 'DCM', 'HCM', 'ICM', 'AMI', 'CHD',  'CS', 'COV' )
DEG <- readRDS("~/PH/results/DEG/DEG_diseases.rds")
deg_s <- subset(DEG,gene %in% features)
##设置分面 
deg_s$types <- ""
deg_s$types[deg_s$gene%in% c("HIF1A", "NDRG1", "PFKFB3", "PDK3", "SERPINE1", "EGFR", "BCL2")] <- "Hypoxia"
deg_s$types[deg_s$gene%in% c("MGAT5", "B4GALT1", "B4GALT5", "ST6GAL1")] <- "Glycosylation"
deg_s$types[deg_s$gene%in% c("VEGFA", "JAG1", "JAG2", "THBD","TIMP1", "FSTL1",'PDGFRB','KDR','TGFBR2','GJC1','TBX20','TGFBR3',"PIK3R3", "HRAS")] <- "Angiogenesis"
deg_s$types[deg_s$gene%in% c("NFKB1", "IL1R1", "IL18R1", "TNFAIP8", "STAT3", "IFNGR2", 'CD40','CD36','ANXA1','HLA-E', 'HLA-DRB1')] <- "Inflammatory response"
deg_s$types[deg_s$gene%in% c('APH1A','DLL4','GATA2','HES1','HEY1','JAG2','PDCD10')] <- "Notch signaling"
deg_s$types[deg_s$gene%in% c('FUS', 'POLDIP2', 'CCNG1', 'PSMD14', 'ERCC1', 'UBC')] <- "Response to DNA damage"
deg_s$types[deg_s$gene%in% c("NDUFS2","NDUFS3", "COX4I1", "UQCRC1", "ATP6AP1", "SDHA", "CYCS")] <- "Oxidative phosphorylation"
deg_s$types[deg_s$gene%in% c("GPX4", "SLC40A1", "FTH1", "FTL")] <- "Ferroptosis"

deg_s$types <- factor(deg_s$types, levels = c("Hypoxia", "Glycosylation", "Angiogenesis",
                                              "Inflammatory response", "Notch signaling","Response to DNA damage",
                                              "Oxidative phosphorylation", "Ferroptosis"))

p <- DotPlot(seurat_obj,features=features[!duplicated(features)], group.by = 'disease')
data <- p[['data']]
colnames(data) <- c("avg.exp", "pct.exp", "gene", "cluster", "avg.exp.scaled")
deg_s <-merge(deg_s,data, by=c('gene','cluster'))
colnames(deg_s)
deg_s$avg_log2FC[deg_s$avg_log2FC < -1.4] = -1.4
deg_s <- deg_s |> filter(abs(avg_log2FC) > 0.25)|> filter(p_val<0.05)

p <- deg_s %>%
  catdotplot(
    x = factor(gene, level = features),
    y = factor(cluster, level = c('ACM','DCM','HCM','ICM','AMI','CHD','CS','COV' )),
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
ggsave(file.path(outdir2,'Dotplot_pathway_genes.pdf'), p, height=3, width=12)

#### 5.Fig2c ----
classic1 <- c('ArtECs','CapECs','VenECs','LECs','EndoECs')
types <- c("Normal","ACM","DCM","HCM","ICM","AMI","CHD","CS","COV")

df_ja=c()
for(j in classic1){
  seurat_obj<-subset(seurat_obj0, classic1==j)
  
  Idents(seurat_obj) <- seurat_obj$disease
  df_gene=FindAllMarkers(seurat_obj,only.pos = T,logfc.threshold =0.25)
  table(df_gene$cluster)
  
  ja=c()
  for (i in types[2:9]) {
    a=df_gene[df_gene$cluster=='Normal',]
    a=a$gene
    b=df_gene[df_gene$cluster==i,]
    b=b$gene
    jaccard=length(intersect(a,b))/length(union(a,b))
    
    ja=c(ja,jaccard)
  }
  df_ja <- cbind(df_ja,ja)
}
rownames(df_ja)=types[2:9]
colnames(df_ja)=classic1
saveRDS(df_ja,paste0(outdir,'df_jaccard.rds'))

##
df <- as.data.frame(t(df_ja))
df$max <- 0.04
df$min <- 0.01
df$mean <- rowMeans(df, na.rm= TRUE)
data <- as.data.frame(t(df))
data <- data[9:11,]

## 绘图 
library(fmsb)
radarchart(data)
pdf(file=file.path(outdir2, 'similarity.pdf'),onefile = FALSE,width=5.5,height=5)
radarchart(data,
           pcol = "#00AFBB",
           pfcol =  scales::alpha("#00AFBB", 0.5),
           plty = "solid",
           cglty = "solid",
           cglcol = "black",
           cglwd =0.5)
dev.off()

#### 6.FigS6a ----
seurat_obj <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
Idents(seurat_obj) <- seurat_obj$classic1

HVG <- FindAllMarkers(seurat_obj)
saveRDS(HVG, file.path(outdir, 'HVG_classic.rds'))
HVG <- readRDS(file.path(outdir, 'HVG_classic.rds'))

library(parallel)
set.seed(123)
for( i in c('ArtECs','CapECs','EndoECs','LECs','VenECs' )){
  case = i
  Genes <- HVG |> filter(cluster == case) |> pull(gene)
  Cells <- seurat_obj@meta.data |> filter(classic1 == case) |> pull(ID)
  sub_obj <- seurat_obj[Genes, Cells]
  expr_matrix <- sub_obj@assays$RNA@data
  disease_labels <- sub_obj@meta.data$disease
  expr_matrix_list <- split(seq_len(ncol(expr_matrix)), disease_labels)
  disease_expr_list <- lapply(expr_matrix_list, function(cells) {
    expr_matrix[, cells, drop = FALSE]
  })
  names(disease_expr_list) <- unique(disease_labels)
  qsave(disease_expr_list, file.path(outdir, paste0(case,  '_expr_list.qs')))
  
  ## calculate
  disease_expr_list <- qread(file.path(outdir, paste0(case,  '_expr_list.qs')))
  library(Matrix)

  calculate_mean_expression <- function(expr_list) {
    mean_expr_list <- lapply(expr_list, function(mat) {
      rowMeans(as.matrix(mat), na.rm = TRUE)
    })
    return(mean_expr_list)
  }
  
  calculate_similarity <- function(mean_expr_list) {
    expr_matrix <- do.call(cbind, mean_expr_list)
    similarity_matrix <- cor(expr_matrix, use = "pairwise.complete.obs")
    overall_similarity <- mean(similarity_matrix, na.rm = TRUE)
    return(overall_similarity)
  }
  
  permutations <- 1000
  permutation_test <- function(disease_expr_list, permutations = 1000) {
    observed_similarity <- calculate_similarity(calculate_mean_expression(disease_expr_list))
    
    permuted_similarities <- numeric(permutations)
    for (i in 1:permutations) {
      permuted_list <- lapply(disease_expr_list, function(mat) {
        permuted_mat <- as.matrix(mat)
        permuted_mat[] <- sample(permuted_mat)
        return(permuted_mat)
      })
      permuted_similarities[i] <- calculate_similarity(calculate_mean_expression(permuted_list))
    }
    
    p_value <- mean(permuted_similarities >= observed_similarity)
    return(list(observed_similarity = observed_similarity, p_value = p_value))
  }
  
  results <- permutation_test(disease_expr_list, permutations = 1000)
  print(results$observed_similarity)
  print(results$p_value)
  results$group <- case
  write.table(results, file.path(outdir, paste0(case,  '_similarity.txt')), row.names = F)
}

##
data <- read.table('./files/pearson.txt', header = T)
data <- data[1:5,-2]
rownames(data) <- data$group
data$max <- 1
data$min <- 0.6
data$observed_similarity <- as.numeric(data$observed_similarity)
data <- data[,-2]
pdata <- as.data.frame(t(data))
pdata <- pdata[c(2,3,1),]

## plot
library(fmsb)
radarchart(pdata)
pdf(file=file.path(outdir2, 'similarity_pearson.pdf'),onefile = FALSE,width=5.5,height=5)
radarchart(pdata,
           pcol = "#00AFBB",
           pfcol =  scales::alpha("#00AFBB", 0.5),
           plty = "solid",
           cglty = "solid",
           cglcol = "black",
           cglwd =0.5)
dev.off()

#### 7.FigS6b ----
prop.table(table(seurat_obj$classic1, seurat_obj$sample), margin = 2)
prop <- as.data.frame(prop.table(table(seurat_obj$classic1, seurat_obj$sample), margin = 2))
colnames(prop) <- c('classic1','sample','prop')
meta <- data.frame(sample=seurat_obj$sample,
                   group=seurat_obj$group)
meta$group <- as.character(meta$group)
meta <- meta %>% distinct(sample,group, .keep_all = T)

prop$group <- ''
i <- 1
for (i in 1:length(rownames(prop))) {
  prop$group[i] <- meta[meta$sample == prop$sample[i],]$group
}

## 
prop$group <- factor(prop$group, levels = c('Health', 'CVD'))
pdata <- prop
fdata <- pdata[pdata$prop > 0,]
data <- fdata
Health <- data %>% filter(group == "Health")
CVD <- data %>% filter(group == "CVD")

ggplot()+
  geom_half_violin(data = Health, aes(x = classic1, y = prop),
                   position = position_dodge(width = 1),
                   scale = 'width',
                   color = 'black', fill = "#B3DE69",
                   side = 'l', size = 0.3)+
  geom_half_violin(data = CVD, aes(x = classic1, y = prop),
                   position = position_dodge(width = 1),
                   scale = 'width',
                   color = 'black', fill = "#b55f31",
                   side = 'r', size = 0.3)+
  stat_compare_means(data = data, aes(x = classic1, y = prop, group = group),
                     method = 'wilcox.test',
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "-")),
                     label = "p.signif")+
  stat_summary(data = prop, aes(x = classic1,y = prop, group = group),
               fun = mean,
               fun.min = function(x){quantile(x)[2]},
               fun.max = function(x){quantile(x)[4]},
               geom = "pointrange",
               size=0.2,
               position = position_dodge(width = 0.2))+
  theme_bw()+
  theme(axis.text.x = element_text(colour = 'black', hjust = 1, angle = 45, size = 10), # hjust = 1, angle = 45,
        axis.text.y = element_text(colour = 'black'),
        legend.position = "top",
        legend.justification = "right")+
  xlab("") + ylab("Proporation")
ggsave(paste0('./plots/classic1_prop.pdf'),last_plot(),width = 4,height = 3)

#### 8.FigS6c ----
#### 8.0 deg ----
outdir3 <- './deg/'
for(i in classic1){
  seurat_obj <- seurat_obj0[ , seurat_obj0$classic1 %in% i]
  Idents(seurat_obj) <- 'group'
  markers <- FindMarkers(seurat_obj,
                         ident.1 = 'CVD',
                         ident.2 = "Health",
                         min.pct = 0,
                         logfc.threshold = 0)
  markers$celltype <- i
  markers$gene <- rownames(markers)
  
  logFC.cutoff <- 0.25
  pvalue.cutoff <- 0.05
  
  markers$change <-
    as.factor(ifelse(
      markers$p_val_adj < pvalue.cutoff &
        abs(markers$avg_log2FC) > logFC.cutoff,
      ifelse(markers$avg_log2FC > logFC.cutoff , 'UP', 'DOWN'),
      'NOT'
    ))
  saveRDS(markers, file.path(outdi3r, paste0(i,"_deg.rds")))
  write.table(markers, file.path(outdir3, paste0(i,"_gsea_input.txt")), quote = F, sep = ",", row.names = F)
}

#### 8.1 FigS6c ----
inputdir <- './deg/'
for (i in classic1) {
  assign(paste0('deg_',i),readRDS(paste0(inputdir,i,'_deg.rds')) |>
           filter(change=='UP') |> pull(gene))
}

classic1_list <- list(`deg_ArtECs`,`deg_CapECs`,`deg_VenECs`,`deg_LECs`,`deg_EndoECs`)
names(classic1_list) <- classic1
upsetData=fromList(classic1_list)

classic1_colors <- c("#a1d5b9","#f2ccac","#edd064","#57C3F3",'#6778AE')
pdf(file=file.path(outdir2, 'classic1_upset.pdf'),onefile = FALSE,width=5,height=4)
upset(upsetData,
      nsets = length(classic1), 
      nintersects = 15,
      order.by = "freq",
      show.numbers = "yes",
      number.angles = 0,
      point.size = 3,
      matrix.color="#b0b9b8",                  
      line.size = 1,                   
      mainbar.y.label = "Gene Intersections", 
      sets.x.label = "Set Size",           
      sets.bar.color = classic1_colors, 
      main.bar.color ="black",
      mb.ratio = c(0.55,0.45), 
      queries = list(
        list(query=intersects, params=list(classic1), color= "#1fab89", active=T),
        list(query=intersects, params=list('ArtECs','CapECs','VenECs'), color= "#1fab89", active=T),
        list(query=intersects, params=list('ArtECs','CapECs','VenECs','LECs'), color= "#1fab89", active=T)
      )
)
dev.off()

## inter_genes
DEG_3 <- classic1_list[!names(classic1_list) %in% c("EndoECs", "LECs")]
inter_3 <- Reduce(intersect, DEG_3)
inter_3_filter <- inter_3[!inter_3 %in% c(`deg_EndoECs`,`deg_LECs`)]
DEG_4 <- classic1_list[names(classic1_list) != "EndoECs"]
inter_4 <- Reduce(intersect, DEG_4)
inter_4_filter <- inter_4[!inter_4 %in% `deg_EndoECs`]
DEG_5 <- classic1_list
inter_5_filter <- Reduce(intersect, DEG_5)
inter_45_filter <- unique(union(inter_4_filter,inter_5_filter))
DEG_celltype <- unique(union(inter_3_filter,inter_45_filter))
saveRDS(DEG_celltype,paste0(outdir,"DEG_celltype.rds"))

#### 9.Fig2d ----
deg_s <- readRDS("~//PH/results/DEG/DEG_diseases.rds")
genes <- unique(deg_s$gene)
DEG <- readRDS("~/PH/results/DEG/DEG_diseases.rds") |> dplyr::filter(change == "UP")
data <- as.data.frame(table(DEG$gene))
count <- data.frame(gene=genes,
                    freq='')
for (i in genes) {
  count[count$gene==i,]$freq <- ifelse(i %in% data$Var1, data[data$Var1 == i,]$Freq, 0)
}
count <- count[order(count$gene),]

df_deg <- data.frame()
for (i in types) {
  deg <- deg_s[deg_s$cluster == i,]
  deg <- deg[order(deg$avg_log2FC),]
  deg$rank <- 1:nrow(deg)
  df_deg <- rbind(df_deg,deg)
}

genes <- unique(deg_s$gene)
df_rank <- data.frame()
for (i in genes) {
  df <- df_deg[df_deg$gene == i,]
  a <- sum(df$rank)/length(genes)
  df_rank <- rbind(df_rank,data.frame(gene=i,
                                      avg_rank=a))
}
df_rank <- df_rank[order(df_rank$gene),]

library(ggplot2)
library(gcookbook)
adata <- cbind(count,df_rank)
adata <- adata[,c(3,2,4)]
DEG_disease <- adata[adata$freq>=7 & adata$avg_rank>=7,]$gene
saveRDS(DEG_disease, paste0(outdir, "DEG_disease.rds"))
# SEEP
selected_genes <- intersect(DEG_celltype,DEG_disease) #1369
saveRDS(selected_genes,paste0(outdir,'selected_genes_345adj.rds'))

##
genes <- c('DOCK1','PPARG','VAV3','ELMO1','CBLB','SPRED2','IGF1R', 'NR3C2', 'ZBTB46')
adata$cutoff <- 'none'
adata[adata$gene %in% DEG_disease,]$cutoff <- 'label'
adata$label <- ''
adata[adata$gene %in% genes,]$label <- adata[adata$gene %in% genes,]$gene

ggplot(adata,aes(x=freq,y=avg_rank)) +
  geom_point(aes(colour = factor(cutoff)))+
  scale_color_manual(values = c('red','grey'))+
  geom_vline(xintercept = 8,
             lty = 2, col = "black", lwd = 0.5) +
  geom_hline(yintercept = 7,
             lty = 2, col = "black", lwd = 0.5) +
  ggrepel::geom_text_repel(aes(label = label),max.overlaps = Inf,
                           color = "black",size = 3,fontface = "italic")+
  xlab('Number of disease types') + ylab('Average rank score')+
  scale_y_continuous(limits = c(0,11))+
  theme_classic()+NoLegend()+
  theme(axis.text.x = element_text(colour = 'black'),
        axis.text.y = element_text(colour = 'black'),
        aspect.ratio = 1.2)
ggsave(file.path(outdir2, 'rank_point.pdf'),
       last_plot(), height=5, width=6)

#### 10.Fig2e ----
## Venn
library(ggvenn)
venn_list <- list(type_deg = DEG_celltype, 
                  rank_genes = DEG_disease)
ggvenn(venn_list,show_percentage = F,
       stroke_color = "white",
       fill_color = c("#156077","#f46f20"),
       fill_alpha = 0.8,
       set_name_color = "black",
       text_size = 6,
       text_color = "black")
ggsave(paste0(outdir2,"venn_selected_genes.pdf"), last_plot(), height = 3.5, width = 4)
## function
selected_genes <- readRDS(paste0(outdir, "selected_genes_345adj.rds"))
genelist <- selected_genes
# GO
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
saveRDS(term, file.path(outdir, 'SEEP_GO.rds'))
# KEGG
convert <- bitr(genelist, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
genelist <- convert$ENTREZID
library(R.utils)
R.utils::setOption("clusterProfiler.download.method",'auto')
ego <- enrichKEGG(
  gene = genelist,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
)
term <- ego@result
saveRDS(term, file.path(outdir, 'selected_genes_KEGG.rds'))

#### 
termgo <- readRDS(paste0(outdir,'SEEP_GO.rds'))
termkegg <- readRDS(paste0(outdir,'selected_genes_KEGG.rds'))
term <- rbind(termgo,termkegg)
term <- term[order(-term$pvalue),]
terms <- c("cell-matrix adhesion", "stress-activated MAPK cascade", "response to transforming growth factor beta", 
           "Wnt signaling pathway", "hippo signaling",
           "Cellular senescence", "Inositol phosphate metabolism", "VEGF signaling pathway","regulation of DNA repair",
           "cellular response to angiotensin")
data <- term[term$Description %in% terms,]
data$Description <- factor(data$Description, levels = terms)
data <- data[c(-5),]

data$formatted_geneID <- sapply(data$geneID, function(gene_ids) {
  paste0('"', unlist(strsplit(gene_ids, "/")), '"', collapse = ", ")
})
# ID transform
data$formatted_geneID[c(2,5,6)] <- sapply(data$geneID[c(2,5,6)], function(gene_ids) {
  entrez_ids <- unlist(strsplit(gene_ids, "/"))
  converted <- bitr(entrez_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
  paste0('"', converted$SYMBOL, '"', collapse = ", ")
})
cat(data[1, 'formatted_geneID'])

for (i in seq_along(data$formatted_geneID)) {
  cat(data$formatted_geneID[i], "\n")
}

data$genes <- c(
  # "stress-activated MAPK cascade"
  "MAPK1, MAPK8, MAPK8IP3, MAP2K4",
  #"cell-matrix adhesion"
  "ABL1, ACER2, ADAMTS9, BCR, CDH13",
  # "response to transforming growth factor beta"
  "SMAD2, SMAD3, SKI",
  # "hippo signaling"
  "YAP1, TEAD4, LATS2",
  # "Inositol phosphate metabolism"
  "PIK3CB, PLCG1, PIK3C2A",
  # "Cellular senescence"
  "RB1, FOXO1, FOXO3",
  # "Wnt signaling pathway"
  "TCF7L2, GSK3B, LRP6, BTRC",
  # "regulation of DNA repair"
  "ATR, ARID1A, SMARCA2, TRIP12",
  # "VEGF signaling pathway"
  "AKT3, RAF1, PIK3CB",
  # "cellular response to angiotensin"
  "NFKB1, ROCK1, ROCK2, ARID1B"
)
data$logp <- -log10(data$pvalue)
data$labelx=rep(0,nrow(data))
data$labely=seq(nrow(data),1)
ggplot(data = data, aes(logp, reorder(Description,logp))) +
  geom_bar(stat="identity",fill='#a3daff',alpha=1,width = 0.6) +
  geom_text(data = data,
            aes(x = 0.1, y = rev(Description), label = genes,),
            size = 4, hjust = 0, vjust = 2.6,
            fontface = 'italic',color = '#2166AC') +
  geom_text(aes(x = 0.1, y = Description, label = Description),
            size=4, hjust =0)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,size=14),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(colour = 'black',size = 10,vjust = 3),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  scale_y_discrete(expand = c(0.1,0))+
  scale_x_continuous(expand = c(0,0))+
  NoLegend()+
  xlab("-log10(pvalue)")
ggsave(file.path(outdir2, 'func_Barplot.pdf'),
       last_plot(), height=5, width=4.5)

#### 11.FigS6d ----
bp <-
  enrichGO(
    inter_3_filter,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2)
term <- bp@result
term <- term[order(-term$pvalue),]
terms <- c(
  'cell-matrix adhesion',
  'response to transforming growth factor beta',
  'activation of innate immune response',
  'regulation of I-kappaB kinase/NF-kappaB signaling',
  'I-kappaB kinase/NF-kappaB signaling',
  'Wnt signaling pathway',
  'coronary vasculature development',
  'heart contraction',
  'cardiac cell development'
)
data <- term[term$Description %in% terms,]
data$Description <- factor(data$Description, levels = terms)
data$logp <- -log10(data$pvalue)
data$labelx=rep(0,nrow(data))
data$labely=seq(nrow(data),1)
ggplot(data = data,aes(logp, reorder(Description,logp))) +
  geom_bar(stat="identity", fill="#f1bbba", alpha=1, width = 0.8) +
  geom_text(aes(x = 0.1, y = Description,label = Description),
            size=4, hjust =0)+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,size=14),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_text(colour = 'black',size = 10,vjust = 3),
        axis.line.x = element_line(colour = 'black', linewidth = 1),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black', linewidth = 1),
        axis.title.x = element_text(colour = 'black', size = 12))+
  scale_y_discrete(expand = c(0.1,0))+
  scale_x_continuous(expand = c(0,0))+
  NoLegend()+
  xlab("-log10(pvalue)")
ggsave(file.path(outdir2, 'inter3_func_Barplot.pdf'), last_plot(), height=4, width=4.5)

#### 12.FigS6e
## VenECs
genes <- setdiff(classic1_list[[3]],
                 union(classic1_list[[1]],
                       c(classic1_list[[2]],classic1_list[[4]],classic1_list[[5]])))
bp <-
  enrichGO(
    genes,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2)
term <- bp@result
terms <- c('cell-substrate junction organization',
           'regulation of GTPase activity',
           'regulation of chemotaxis',
           'regulation of cell-substrate adhesion',
           'regulation of leukocyte migration',
           'cellular response to oxidative stress',
           'regulation of Wnt signaling pathway',
           'leukocyte chemotaxis',
           'extracellular matrix organization'
)

data <- term[term$Description %in% terms,]
data$labelx=rep(0,nrow(data))
data$labely=seq(nrow(data),1)
ggplot(data, aes(x = -log10(pvalue), y = reorder(Description,-log10(pvalue)))) +
  geom_bar(stat="identity", alpha=1, fill= "#a1d5b9", width = 0.8) +
  geom_text(aes(x=labelx, y=labely, label = data$Description),
            size=4.5, hjust =0)+
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
ggsave(paste0(outdir2, 'Ven_Barplot.pdf'),
       last_plot(), height=4.5, width=5.5)

## LECs
genes <- setdiff(classic1_list[[4]],
                 union(classic1_list[[1]],
                       c(classic1_list[[2]],classic1_list[[3]],classic1_list[[5]])))
bp <-
  enrichGO(
    genes,
    OrgDb = org.Hs.eg.db,
    keyType = 'SYMBOL',
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2)
term <- bp@result
terms <- c('lymph vessel development',
           'positive regulation of MAPK cascade',
           'phosphatidylinositol-mediated signaling',
           'small GTPase mediated signal transduction',
           'Notch signaling pathway',
           'regulation of cell-cell adhesion',
           'lymphangiogenesis',
           'blood vessel endothelial cell differentiation',
           'cell-matrix adhesion'
)

data <- term[term$Description %in% terms,]
data$labelx=rep(0,nrow(data))
data$labely=seq(nrow(data),1)
ggplot(data, aes(x = -log10(pvalue), y = reorder(Description,-log10(pvalue)))) +
  geom_bar(stat="identity", alpha=1, fill= "#57C3F3", width = 0.8) +
  geom_text(aes(x=labelx, y=labely, label = data$Description), size=4.5, hjust =0)+
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
ggsave(paste0(outdir2, 'LEC_Barplot.pdf'),
       last_plot(), height=4.5, width=5.5)

#### 12.FigS6f ----
disease_colors <- c("#B3DE69","#FB8072", "#80B1D3", "#FDB462", "#8DD3C7",
                    "#FCCDE5","#CCEBC5", "#377EB8", "#1B9E77") 
age_colors <- c('#FB8072', '#1B9E77', '#80B1D3', '#BC80BD', '#FDB462', '#377EB8')
sex_colors <- c("#ee6470", "#00a6e1")

cells_ranking <- AUCell_buildRankings(seurat_obj0@assays$RNA@counts)
cells_AUC <- AUCell_calcAUC(selected_genes, cells_ranking, aucMaxRank=nrow(cells_ranking)*0.1)
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
seurat_obj0$AUC  <- AUCell_socre
data <- seurat_obj0@meta.data

## 12.1 disease
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
       ggplot2::last_plot(),width = 6, height = 2.5)

## 12.2 age
data_D <- data[data$group == 'CVD',]
sub_data <- data_D %>%
  dplyr::filter(age_group %in% c('20-30','30-40','40-50','50-60','60-70','>70'))
sub_data$age_group <- factor(sub_data$age_group, levels = c('20-30','30-40','40-50','50-60','60-70','>70'))

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
ggsave(file.path(outdir2, 'Score_vlnplot_age_D.pdf'),
       ggplot2::last_plot(),width = 4, height = 2.5)

## 12.3 sex
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

#### 13.FigS6g ----
sub_obj <- subset(seurat_obj0, Article %in% c('Data5_HF_AndrewLKoenig_2022'))
cells_ranking <- AUCell_buildRankings(sub_obj@assays$RNA@counts)
cells_AUC <- AUCell_calcAUC(selected_genes, cells_ranking, aucMaxRank=nrow(cells_ranking)*0.1)
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
sub_obj$AUC  <- AUCell_socre

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

####
meta <- sub_obj@meta.data
data <- meta[,c('donor','sample','cell_type','age_group','classic1','disease','AUC','LVEF')]
adata <-aggregate(data[, c('AUC','LVEF')], list(data$sample), FUN = mean)
rownames(adata) <- adata$Group.1
adata <- adata[,-1]

correlation <- cor.test(adata$AUC, adata$LVEF, method = "spearman", alternative = "less")
cor_coef <- round(correlation$estimate, 3)
p_value <- round(correlation$p.value, 4)
ggplot(adata, aes(x = AUC, y = LVEF)) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.4) +
  geom_point(color = 'grey') +
  theme_classic() +
  labs(title = "", x = 'SEEP score', y = "LVEF(%)") +
  annotate("text", x = min(adata$AUC, na.rm = TRUE), y = max(adata$LVEF, na.rm = TRUE),
           label = paste("Correlation =", cor_coef, "\np-value =", p_value), 
           hjust = 0, vjust = 1, size = 4)
ggsave(file.path(outdir2, paste0('AUC_LVEF_cor.pdf')), ggplot2::last_plot(), height=3.5, width=3.5)

#### 14.Fig2f ----
library(readxl)
folder_path <- './database/'
file_list <- list.files(folder_path)

#### 14.1 load database ----
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

#### 14.2 intersect genes ----
disgenes <- unique(result$gene)
saveRDS(disgenes, paste0(outdir, 'disgenet_disgenes.rds'))

genes_inter <- intersect(disgenes, selected_genes)
DG_list <- genelist[genelist$gene %in% genes_inter,]
write.csv(DG_list,"./files/DG_list.csv")

#### 15.Fig2g ----

#### 15.1 EC ----
## 载入数据
seurat_obj <- qs::qread("~/PH/results/named/files/sample_15_1800_named_type.qs")
datalist <- list(exp=seurat_obj@assays[["RNA"]]@counts,
                 annot=seurat_obj@meta.data)
annotLevels = list(level1class=datalist$annot$disease,level2class=datalist$annot$cell_type)
ctd = generate_celltype_data(exp=datalist[[1]],annotLevels=annotLevels,groupName="filter")
print(ctd)
load(ctd)

#### Application to genetic data
human.hits <- selected_genes
human.bg = unique(rownames(seurat_obj@assays[["RNA"]]@counts))
reps=14819
level=2 

# Bootstrap significance testing, without controlling for transcript length and GC content
full_results = bootstrap_enrichment_test(sct_data=ctd, hits=human.hits, bg=human.bg,
                                         reps=reps, annotLevel=level,
                                         genelistSpecies = "human",
                                         sctSpecies = "human",
                                         sctSpecies_origin = "human",
                                         output_species = "human",)

#### 15.2 barplot_EC ----
full_results$results$CellType <- factor(full_results$results$CellType, 
                                        levels = c("Art_c01_DKK2","Art_c02_NEBL","Art_c03_VEGFA",
                                                   "Cap_c04_PAPSS2","Cap_c05_RGCC","Cap_c06_RDH10",
                                                   "Cap_c07_CA2","Cap_c08_TMEM163","Cap_c09_RGS6",
                                                   "Ven_c10_VCAM1","Ven_c11_ACTA2","LEC_c12_FLT4","Endo_c13_NPR3"))
pdf(paste0(outdir2,'ewce_bar_endo.pdf'), width = 6, height = 4)
ewce_barplot(full_results$results,mtc_method="BH")
dev.off()

## 整理结果
meta <- full_results[["gene_data"]] |> as.data.frame()
meta <- meta[,c('gene','CellType','p')]
pivot_df <- meta %>%
  pivot_wider(names_from = CellType, values_from = p, values_fill = list(p = 0)) %>%
  as.data.frame()
rownames(pivot_df) <- pivot_df$gene
pivot_df <- pivot_df[,-1]
pivot_df <- pivot_df[features,]

celltypes <- c("Art.c01.DKK2","Art.c02.NEBL","Art.c03.VEGFA","Cap.c04.PAPSS2","Cap.c05.RGCC","Cap.c06.RDH10",
               "Cap.c07.CA2","Cap.c08.TMEM163","Cap.c09.RGS6","Ven.c10.VCAM1","Ven.c11.ACTA2","LEC.c12.FLT4","Endo.c13.NPR3")
names <- colnames(pivot_df)
names <- gsub("_", "\\.", names)
colnames(pivot_df) <- names
pivot_df <- pivot_df[,celltypes]

new_df <- pivot_df
new_df[new_df == 1] <- " "
new_df[new_df == 0] <- "*"

df_endo <- new_df

#### 15.3 alltypes ----
seurat_obj <- qread("~/PH/results/named/files/Pan_all.qs")
seurat_obj$type <- seurat_obj$cell_type
seurat_obj$type[seurat_obj$type %in% celltypes] <- 'Endothelial'

datalist <- list(exp=seurat_obj@assays[["RNA"]]@counts,
                 annot=seurat_obj@meta.data)
annotLevels = list(level1class=datalist$annot$cell_type,level2class=datalist$annot$type)
ctd = generate_celltype_data(exp=datalist[[1]],annotLevels=annotLevels,groupName="filter")
print(ctd)
load(ctd)

#### Application to genetic data
human.hits <- selected_genes
human.bg = unique(rownames(seurat_obj@assays[["RNA"]]@counts))
reps=14819
level=2

# Bootstrap significance testing, without controlling for transcript length and GC content
full_results = bootstrap_enrichment_test(sct_data=ctd, hits=human.hits, bg=human.bg,
                                         reps=reps, annotLevel=level,
                                         genelistSpecies = "human",
                                         sctSpecies = "human",
                                         sctSpecies_origin = "human",
                                         output_species = "human",)

#### 15.4 barplot_alltypes ----
alltypes <- c('Cardiomyocyte' ,'Endothelial' ,'Fibroblast' ,'Lymphoid' ,'Myeloid')
full_results$results$CellType <- factor(full_results$results$CellType, 
                                        levels = alltypes)
pdf(paste0(outdir2,'ewce_bar_alltypes.pdf'), width = 3.5, height = 4)
ewce_barplot(full_results$results,mtc_method="BH")
dev.off()

##
meta <- full_results[["gene_data"]] |> as.data.frame()
meta <- meta[,c('gene','CellType','p')]
pivot_df <- meta %>%
  pivot_wider(names_from = CellType, values_from = p, values_fill = list(p = 0)) %>%
  as.data.frame()
rownames(pivot_df) <- pivot_df$gene
pivot_df <- pivot_df[,-1]
pivot_df <- pivot_df[features,]

new_df <- pivot_df
new_df[new_df == 1] <- " "
new_df[new_df == 0] <- "*"
df_type <- new_df[,c('Endothelial', 'Cardiomyocyte', 'Fibroblast', 'Lymphoid', 'Myeloid')]

##
df_all <- cbind(df_endo,df_type)
df_all <- df_all[,-14]

#### 15.5 heatmap ----
adata <- t(as.data.frame(seurat_obj@assays$RNA@data[features,]))
adata <- as.data.frame(adata)
adata$cell_type <- seurat_obj$cell_type

pdata <- aggregate(adata[, -length(colnames(adata))], list(adata$cell_type), FUN = mean)
rownames(pdata) <- pdata[,1]
pdata <- pdata[,-1]
pdata <- pdata[c('Cardiomyocyte', 'Fibroblast', 'Lymphoid', 'Myeloid',
                 "Art.c01.DKK2","Art.c02.NEBL","Art.c03.VEGFA","Cap.c04.PAPSS2","Cap.c05.RGCC","Cap.c06.RDH10",
                 "Cap.c07.CA2","Cap.c08.TMEM163","Cap.c09.RGS6","Ven.c10.VCAM1","Ven.c11.ACTA2","LEC.c12.FLT4","Endo.c13.NPR3"
),]
#
pdata_f <- t(pdata)
data_scale <- round(t(apply(pdata_f, 1, scale)),2)
colnames(data_scale) <- colnames(pdata_f)
##
color_palette <- c(colorRampPalette(colors = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0"))(50),
                   colorRampPalette(colors = c("#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(50))
p <- pheatmap(data_scale, scale="none",
              cluster_cols = F,  cluster_rows = F,
              gaps_col = c(4,7,13,15,16), gaps_row = 14,
              display_numbers = df_all, fontsize_number = 10,
              border="white", cellwidth = 10,cellheight =13,
              color = color_palette,
              fontsize_row = 10,fontsize_col = 10)
ggsave(file.path(outdir2, 'ewce_heatmap.pdf'),
       p, limitsize = FALSE, height=8, width=6)

#### 16.Fig2h ----
library(sceasy)
library(Seurat)
library(reticulate)

#### 16.0 transfer ----
use_condaenv("sceasy", required = TRUE)
loompy <- reticulate::import('loompy')
seurat_obj <- qread('~/PH/results/files/Pan_all_CVD_D2.qs')
seurat_obj@meta.data <- seurat_obj@meta.data[ ,c("ID", "orig.ident", "disease", "cell_type")]

sceasy::convertFormat(seurat_obj, from="seurat", to="anndata",
                      outFile='/home/chencx/figures/PH_data.h5ad')

seurat_obj@assays$RNA$data <- seurat_obj@assays$RNA$counts
sceasy::convertFormat(seurat_obj, from="seurat", to="anndata",
                      outFile='/home/chencx/figures/PH_raw.h5ad')

#### 16.1 drug2cell ----
#### 16.2 heatmap ----
seurat_obj@meta.data <- seurat_obj@meta.data[ ,c("ID", "orig.ident", "disease", "cell_type")]
library(tidyr)
library(dplyr)
category <- read.csv('./files/drug-target_dicts.csv')
scores <- read.csv('./files/drug_scores.csv')
pvadj <- read.csv('./files/drug_pvals_adj.csv')
logFC <- read.csv('./files/drug_logFC.csv')
names <- read.csv('./files/drug_names.csv')

####
data_list <- list(scores = scores, pvadj = pvadj, logFC = logFC)
data_names <- names(data_list)
for (i in seq_along(data_list)) {
  data <- data_list[[i]]
  data_name <- data_names[i]
  
  scores_long <- data %>%
    mutate(Row = row_number()) %>%
    pivot_longer(-Row, names_to = "CellType", values_to = data_name)
  names_long <- names %>%
    mutate(Row = row_number()) %>%
    pivot_longer(-Row, names_to = "CellType", values_to = "Drug")
  merged_long <- inner_join(names_long, scores_long, by = c("Row", "CellType"))
  write.csv(merged_long, file.path(outdir, paste0(data_name, "_Drag_long.csv")), row.names = FALSE)
}

#### drug-gene ----
library(AUCell)
category <- read.csv('./files/drug-target_dicts.csv')
dataset <- category[,1:2]
colnames(dataset) <- c('gene_symbol', 'gs_name')
dataset <- unique(dataset)

##
geneSets <- lapply(unique(dataset$gs_name), 
                   function(x){dataset$gene_symbol[dataset$gs_name == x]})
names(geneSets) <- unique(dataset$gs_name)
cells_rankings <- AUCell_buildRankings(seurat_obj@assays$RNA@data) 
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
colnames(AUCell_socre) <- colnames(AUCell_socre)
saveRDS(AUCell_socre,file.path(outdir, "AUCscore.rds"))

##
seurat_obj$ID <- rownames(seurat_obj@meta.data)
meta <- AUCell_socre
meta <- as.data.frame(meta)
meta$cell_type <- seurat_obj$cell_type
object <- aggregate(meta[, -length(colnames(meta))], list(meta$cell_type), FUN = mean)
rownames(object) <- object$Group.1
object <-object[,-1]
saveRDS(object,file.path(outdir, "AUCscore_celltypes.rds"))

##
object <- readRDS(file.path(outdir, "AUCscore_celltypes.rds"))
s_drug <- read.table(file.path(outdir, 's_drugs.txt'), sep = ',') |> pull()
s_Genes <- readRDS(paste0(outdir, "selected_genes_345adj.rds"))

s_drug2 <- unique(intersect(category2$Drug, s_drug))
category3 <- category |> filter(Drug %in% s_drug2) |> filter(Category == 'C')
category4 <- category |> filter(Drug %in% s_drug2) |> filter(Category == 'C') |> filter(Gene %in% s_Genes) 
unique(category4$Gene)
unique(category4$Drug)
category5 <- category |> filter(Drug %in% s_drug2)
unique(category5$Gene) 
unique(category5$Drug)
s_drug3 <- c('CHEMBL595|PIOGLITAZONE',
             'CHEMBL843|ROSIGLITAZONE MALEATE',
             'CHEMBL2108675|DUPILUMAB',
             'CHEMBL815|DINOPROST',
             'CHEMBL1484|NICARDIPINE',
             'CHEMBL6966|VERAPAMIL',
             'CHEMBL1294|QUINIDINE'
)
category_s <- category |> filter(Drug %in% s_drug3)|> filter(Gene %in% s_Genes)

##
object2 <- object[, colnames(object) %in% s_drug3]
object2 <- as.data.frame(t(object2))
saveRDS(object2,file.path(outdir, "AUCscore_celltypes_sub.rds"))

object2[object2>0.5] = 0.5
p1 <- pheatmap(object2, treeheight_row = 5,
               cluster_cols = F, cluster_rows = T, scale="none",
               fontsize_number = 20, border="grey",
               shape = "circle",cellwidth = 10,cellheight =10,
               col = c(colorRampPalette(colors = c("white","#FFFFF3","#FDDBC7","#F4A582","#D6604D","#B2182B"))(50)),
               fontsize_row = 10,fontsize_col = 12 
) 
ggsave(file.path(outdir2, 'AUCscpre_sub2.pdf'),
       p1,limitsize = FALSE, height=3, width=8)

