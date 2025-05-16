#### cat_construct_metacells_sub_cell_type ----
cat_construct_metacells_sub_cell_type <- function(seurat_obj, k, name, min_cells = 50, max_shared = 10, assay = "RNA") {
  if ("sub_cell_type" %in% colnames(seurat_obj@meta.data) &
      "sample" %in% colnames(seurat_obj@meta.data)) {
    set.seed(717)
    #### Set up Seurat object for WGCNA ----
    seurat_obj <- hdWGCNA::SetupForWGCNA(
      seurat_obj,
      gene_select = "fraction",
      fraction = 0.05,
      wgcna_name = name
    )
    
    #### Construct metacells ----
    # construct metacells  in each group
    seurat_obj <- hdWGCNA::MetacellsByGroups(
      seurat_obj = seurat_obj,
      group.by = c("sub_cell_type", "sample"),
      k = k,
      min_cells = min_cells,
      max_shared = max_shared,
      mode = "sum",
      assay = assay,
      ident.group = "sub_cell_type"
    )
    seurat_obj <- hdWGCNA::NormalizeMetacells(seurat_obj)
    return(seurat_obj)
  }
  stop("The column name of sub_cell_type or sample does not exist in the meta.data of your Seurat object, please add it!")
}

#### do.tissueDist ----
do.tissueDist <- function(cellInfo.tb = cellInfo.tb,
                          meta.cluster = cellInfo.tb$meta.cluster,
                          colname.patient = "patient",
                          loc = cellInfo.tb$loc,
                          out.prefix,
                          pdf.width=3,
                          pdf.height=5,
                          verbose=0){
  ##input data 
  library(data.table)
  dir.create(dirname(out.prefix),F,T)
  
  cellInfo.tb = data.table(cellInfo.tb)
  cellInfo.tb$meta.cluster = as.character(meta.cluster)
  
  if(is.factor(loc)){
    cellInfo.tb$loc = loc
  }else{cellInfo.tb$loc = as.factor(loc)}
  
  loc.avai.vec <- levels(cellInfo.tb[["loc"]])
  count.dist <- unclass(cellInfo.tb[,table(meta.cluster,loc)])[,loc.avai.vec]
  freq.dist <- sweep(count.dist,1,rowSums(count.dist),"/")
  freq.dist.bin <- floor(freq.dist * 100 / 10)
  print(freq.dist.bin)
  
  {
    count.dist.melt.ext.tb <- test.dist.table(count.dist)
    p.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="p.value")
    OR.dist.tb <- dcast(count.dist.melt.ext.tb,rid~cid,value.var="OR")
    OR.dist.mtx <- as.matrix(OR.dist.tb[,-1])
    rownames(OR.dist.mtx) <- OR.dist.tb[[1]]
  }
  
  sscVis::plotMatrix.simple(OR.dist.mtx,
                            out.prefix=sprintf("%s.OR.dist",out.prefix),
                            show.number=F,
                            waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                            exp.name=expression(italic(OR)),
                            z.hi=4,
                            palatte=viridis::viridis(7),
                            pdf.width = 4, pdf.height = pdf.height)
  if(verbose==1){
    return(list("count.dist.melt.ext.tb"=count.dist.melt.ext.tb,
                "p.dist.tb"=p.dist.tb,
                "OR.dist.tb"=OR.dist.tb,
                "OR.dist.mtx"=OR.dist.mtx))
  }else{
    return(OR.dist.mtx)
  }
}

#### hp_run_pyscenic ----
hp_run_pyscenic <-
  function(x,
           outdir = "./pyscenic",
           pyscenic = NULL,
           species = "human",
           features = NULL,
           tfs_fname = NULL,
           database_fname = NULL,
           annotations_fname = NULL,
           mode = "custom_multiprocessing",
           ncores = 10,
           seed = 717) {
    if (is.null(pyscenic)) {
      pyscenic <- "/opt/pyscenic/0.12.1/bin/pyscenic"
      hp_check_file(pyscenic)
    }
    if (species == "human") {
      tfs_fname <- "/DATA/public/cistarget/tf_lists/allTFs_hg38.txt"
      database_fname <-
        c(
          "/DATA/public/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
          "/DATA/public/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
        )
      annotations_fname <-
        "/DATA/public/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
    } else if (species == "mouse") {
      tfs_fname <- "/DATA/public/cistarget/tf_lists/allTFs_mm.txt"
      database_fname <-
        c(
          "/DATA/public/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
          "/DATA/public/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
        )
      annotations_fname <-
        "/DATA/public/cistarget/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl"
    } else {
      stop("Please check if the species is human or mouse!")
    }
    hp_check_file(x = c(tfs_fname, database_fname, annotations_fname))
    outdir <- hp_set_path(outdir)
    expression_mtx_fname <- file.path(outdir, "expression_mtx.csv")
    module_fname <-
      file.path(outdir, "expression_mtx.adjacencies.tsv")
    signatures_fname <- file.path(outdir, "regulons.gmt")
    auc_mtx_fname <- file.path(outdir, "auc_mtx.csv")
    tfs_target_fname <- file.path(outdir, "tfs_targer.tsv")
    auc_g_mtx_fname <- file.path(outdir, "auc_g_mtx.csv")
    # pre-processing
    if (!is.null(features)) {
      # x <- x[features,]
      tfs <- read.table(tfs_fname)[, 1]
      filtered_tfs <- intersect(tfs, features)
      tfs_fname <- file.path(outdir, "tfs.txt")
      write.table(
        x = filtered_tfs,
        file = tfs_fname,
        row.names = FALSE,
        quote = FALSE,
        col.names = FALSE
      )
      counts <- x@assays$RNA@counts[features,] |>
        as.data.frame()
    } else {
      counts <- x@assays$RNA@counts |>
        as.data.frame()
    }
    if (!file.exists(expression_mtx_fname)) {
      log_info("Write counts to a csv file...")
      tic("Write csv")
      fwrite(
        x = counts,
        file = expression_mtx_fname,
        quote = FALSE,
        row.names = TRUE
      )
      toc()
    }
    # grn
    if (!file.exists(module_fname)) {
      log_info("Derive co-expression modules from expression matrix...")
      tic("Run grn")
      system2(
        command = pyscenic,
        args = c(
          "grn",
          "--transpose",
          "--num_workers",
          ncores,
          "--seed",
          seed,
          "--output",
          module_fname,
          expression_mtx_fname,
          tfs_fname
        )
      )
      toc()
    }
    # ctx
    if (!file.exists(signatures_fname)) {
      log_info(
        "Find enriched motifs for a gene signature and optionally prune targets from this signature based on cis-regulatory cues..."
      )
      tic("Run ctx")
      system2(
        command = pyscenic,
        args = c(
          "ctx",
          "--transpose",
          "--output",
          signatures_fname,
          "--mode",
          mode,
          "--annotations_fname",
          annotations_fname,
          "--num_workers",
          ncores,
          "--expression_mtx_fname",
          expression_mtx_fname,
          module_fname,
          database_fname
        )
      )
      toc()
    }
    # aucell
    if (!file.exists(auc_mtx_fname)) {
      log_info("Quantify activity of gene signatures across single cells...")
      tic("Run aucell")
      system2(
        command = pyscenic,
        args = c(
          "aucell",
          "--transpose",
          "--output",
          auc_mtx_fname,
          "--num_workers",
          ncores,
          expression_mtx_fname,
          signatures_fname
        )
      )
      toc()
    }

    # tfs_target
    if (!file.exists(tfs_target_fname)) {
      signatures <-
        clusterProfiler::read.gmt(signatures_fname)
      signatures_count <- signatures |>
        count(term)
      signatures <-
        left_join(x = signatures, y = signatures_count,
                  by = "term")
      signatures <- signatures |>
        mutate(term = str_split_fixed(
          string = term,
          pattern = "\\(",
          n = 2
        )[, 1]) |>
        mutate(symbol = paste0(term, " (", n, "g)")) |>
        dplyr::select(symbol, term, gene) |>
        rename(target_gene = gene, tf = term)
      write.csv(signatures,
                tfs_target_fname,
                row.names = FALSE)
    } else {
      signatures <- read.csv(tfs_target_fname)
    }
    # auc_g_mtx
    if (!file.exists(auc_g_mtx_fname)) {
      auc_mtx <-
        read.csv(auc_mtx_fname,
                 row.names = 1,
                 check.names = FALSE)
      rownames(auc_mtx) <- unique(signatures$tf)
      fwrite(
        x = auc_mtx,
        file = auc_g_mtx_fname,
        quote = FALSE,
        row.names = TRUE
      )
    } else {
      auc_mtx <- read.csv(auc_g_mtx_fname,
                          row.names = 1,
                          check.names = FALSE)
    }
    colnames(auc_mtx) <-
      gsub(pattern = "\\.",
           replacement = "-",
           x = colnames(auc_mtx))
    x[["scenic"]] <- CreateAssayObject(counts = auc_mtx)

    return(x)
  }

