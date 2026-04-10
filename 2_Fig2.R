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

terms <- c('dephosphorylation',
           'Wnt signaling pathway',
           'small GTPase mediated signal transduction',
           'transforming growth factor beta receptor signaling pathway',
           'miRNA metabolic process',
           'I-KappaB kinase/NF-kappaB signaling',
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

features <- c("NDRG1", "PFKFB3",
              'PTPRM', 'PTPN4', 'INPP5D', 'MTMR3','PPP6R3',
              'NFKB1', 'IKBKB', 'MAPK14', 'MAPK8',
              'CD44','ITGA4','CD276','EMILIN2',
              'MYH6','TNNT2','DSP','MYBPC3','SCN5A',
              'VEGFA', 'KDR', 'ANGPT2', 'DLL4', 'CXCL12', 'GATA4',
              'NOTCH3', 'HES1', 'HEY1', 'JAG2',
              'NDUFS7', 'SDHB', 'CYCS', 'COX4I1', 'UQCRFS1',
              'ACSL4',"GPX4",'AIFM2','DHODH','SLC40A1', 'FTH1', 'FTL')
features <- features[!duplicated(features)]

types <- c('ACM', 'DCM', 'HCM', 'ICM', 'AMI', 'CHD', 'CS', 'COV' )
DEG <- readRDS("~/PH/results/DEG/files/DEGs_disease_number.rds")
DEG$gene <- DEG$Gene_Symbol
DEG$avg_log2FC <- DEG$logFC
deg_s <- subset(DEG, gene %in% features)

deg_s$types <- ""
deg_s$types[deg_s$gene%in% c("NDRG1", "PFKFB3")] <- "Hypoxia"
deg_s$types[deg_s$gene%in% c('NFKB1', 'IKBKB', 'MAPK14', 'MAPK8')] <- "Inflammation"
deg_s$types[deg_s$gene%in% c('VEGFA', 'KDR', 'ANGPT2', 'DLL4', 'CXCL12', 'GATA4')] <- "Angiogenesis"
deg_s$types[deg_s$gene%in% c('PTPRM', 'PTPN4', 'INPP5D', 'MTMR3','PPP6R3')] <- "Dephosphorylation"
deg_s$types[deg_s$gene%in% c('MYH6','TNNT2','DSP','MYBPC3','SCN5A')] <- "Cardiac contraction"
deg_s$types[deg_s$gene%in% c('CD44','ITGA4','CD276','EMILIN2')] <- "Cell adhesion"
deg_s$types[deg_s$gene%in% c('NOTCH3','HES1','HEY1','JAG2')] <- "NOTCH signaling"
deg_s$types[deg_s$gene%in% c('NDUFS7', 'SDHB', 'CYCS', 'COX4I1', 'UQCRFS1')] <- "Energy metabolism"
deg_s$types[deg_s$gene%in% c('ACSL4',"GPX4",'AIFM2','DHODH','SLC40A1', 'FTH1', 'FTL')] <- "Ferroptosis"

deg_s$types <- factor(deg_s$types, levels = c("Dephosphorylation", "Hypoxia", "Inflammation", 
                                              'Cardiac contraction',"Cell adhesion","Energy metabolism", "Angiogenesis",
                                              'NOTCH signaling',"Ferroptosis"))

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




#### inter genes ----
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

#### 4.Fig2f AUCell ----
library(AUCell)
seurat_obj0 <- readRDS("~/project/PH/src_bulk/data/seurat_obj.rds")
selected_genes <- readRDS(paste0(outdir, "cutgene_cut6.rds")) 
cells_ranking <- AUCell_buildRankings(seurat_obj0@assays$RNA@counts)
cells_AUC <- AUCell_calcAUC(selected_genes, cells_ranking, aucMaxRank=nrow(cells_ranking)*0.1)
AUCell_socre <- as.data.frame(t(cells_AUC@assays@data@listData[["AUC"]]))
seurat_obj0$AUC  <- AUCell_socre
data <- seurat_obj0@meta.data

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

library(ggsignif)
data_D <- data[data$group == 'CVD',]
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

#### 5.Fig2g ----
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
  labs(title = "",x = 'SEEP score',y = "LVEF(%)") +
  annotate("text", x = min(adata$AUC, na.rm = TRUE), y = max(adata$LVEF, na.rm = TRUE),
           label = cor_label, hjust = 0, vjust = 1, size = 4)
ggsave(file.path(outdir2, paste0('5_AUC_LVEF_cor.pdf')),
       last_plot(),height=3.5, width=3.5)



#### 6.Fig2i ----------------------------
library(Seurat)
library(ggplot2)
library(tidyverse)
library(qs)
library(AUCell)
library(pheatmap)
setwd('~/project/PH/drug2cell/results/2603/')
outdir <- './'
outdir2 <- './'

seurat_obj <- qread('~/project/PH/alltype/merge/results/files/T6_ori_merged_CVD_1000.qs')
seurat_obj$type <- ifelse(
  seurat_obj$celltype == "Endothelial",
  as.character(seurat_obj$cell_type),  
  as.character(seurat_obj$celltype)   
)
table(seurat_obj$type)
seurat_obj@meta.data <- seurat_obj@meta.data[ ,c("ID", "orig.ident", "disease", 
                                                 "celltype", "cell_type", 'type')]

category <- read.csv('~/figures/files/drug-target_dicts.csv')
category2 <- category |> filter(Category == 'C')
category2$Drug <- paste(category2$Drug, category2$Category, sep = '_')
category3 <- category |> filter(Category == 'No-category')
dataset <- rbind(category2, category3)
dataset <- dataset[,1:2]
colnames(dataset) <- c('gene_symbol', 'gs_name')
dataset <- unique(dataset) 

AUCell_socre <- readRDS(file.path(outdir, "AUCscore_CVD_1000.rds"))
seurat_obj$ID <- rownames(seurat_obj@meta.data)
meta <- AUCell_socre
meta <- as.data.frame(meta)
meta$cell_type <- seurat_obj$type
object <-
  aggregate(meta[, -length(colnames(meta))], list(meta$cell_type), FUN = mean)
rownames(object) <- object$Group.1
object <- object[,-1]
saveRDS(object,file.path(outdir, "AUCscore_C_celltypes_CVD.rds"))

zero_var_cols <- sapply(object, function(col) var(col) == 0)

zero_var_rows <- apply(object, 1, function(row) var(row) == 0)

sum(zero_var_cols)
sum(zero_var_rows)
object2 <- as.data.frame(t(object[, !zero_var_cols]))
rownames(object2) <- gsub('_C', '', rownames(object2))
object2 <- object2[,c("Art.c01.DKK2", "Art.c02.NEBL", "Art.c03.VEGFA", "Cap.c04.PAPSS2", "Cap.c05.RGCC", 
                      "Cap.c06.RDH10", "Cap.c07.CA2", "Cap.c08.TMEM163", "Cap.c09.RGS6", 
                      "Ven.c10.VCAM1", "Ven.c11.TMSB4X", "LEC.c12.FLT4", "Endo.c13.NPR3", 
                      "Cardiomyocyte", "Fibroblast", "Mural", "Myeloid", "Lymphoid"         
)]
saveRDS(object2, file.path(outdir, "AUCscore_C_celltypes_CVD2.rds"))

object3 <- as.data.frame(object2)
object3$ArtECs <- rowMeans(object3[, 1:3])
object3$CapECs <- rowMeans(object3[, 4:9])
object3$VenECs <- rowMeans(object3[, 10:11])
object3 <- object3[, c('ArtECs', 'CapECs', 'VenECs', 
                       "Cardiomyocyte", "Fibroblast", 
                       "Mural", "Myeloid", "Lymphoid")]
saveRDS(object3, file.path(outdir, "AUCscore_C_class_CVD.rds"))

## 4.heatmap
p1 <- pheatmap(object2,cluster_cols = F,
               cluster_rows = T, scale="row",
               fontsize_number = 20, border="white",
               shape = "circle",cellwidth = 10,cellheight =8,
               color=colorRampPalette(c("darkblue",'white', "darkred"))(10000),
               fontsize_row = 10,fontsize_col = 12 ) 

ggsave(file.path(outdir2,'AUCscore_C_all_CVD_celltype.pdf'),
       p1,limitsize = FALSE,
       height=120, width=10)

mat <- mat + 1e-8
p2 <- pheatmap(object3,cluster_cols = F,
               cluster_rows = T, scale="row",
               fontsize_number = 20, border="white", 
               shape = "circle",cellwidth = 10,cellheight =8,
               color=colorRampPalette(c("darkblue",'white', "darkred"))(10000),
               fontsize_row = 10,fontsize_col = 12 ) 

ggsave(file.path(outdir2,'AUCscore_C_all_CVD_class.pdf'),
       p2,limitsize = FALSE, height=120, width= 8)

s_drug <- c('CHEMBL1484|NICARDIPINE',
            'CHEMBL772|RESERPINE',
            'CHEMBL254219|DIGITOXIN', 
            'CHEMBL1463345|CANRENONE', 
            'CHEMBL843|ROSIGLITAZONE MALEATE',
            'CHEMBL1567|SUNITINIB MALATE')
object4 <- object2[intersect(s_drug, rownames(object2)), ]
saveRDS(object4, file.path(outdir, "AUCscore_C_celltype_CVD_s.rds"))

p3 <- pheatmap(object4,cluster_cols = F, cluster_rows = F, scale="row",
               fontsize_number = 20, border="white", 
               shape = "circle",cellwidth = 10,cellheight =8,
               color=colorRampPalette(c("darkblue",'white', "darkred"))(10000),
               fontsize_row = 10,fontsize_col = 12 ) 
ggsave(file.path(outdir2,'AUCscore_C_s_CVD_celltype.pdf'),
       p3,limitsize = FALSE, height=3, width= 6)



