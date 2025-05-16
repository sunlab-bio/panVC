#### EWCE ----
check_mtc_method <- function (mtc_method) 
{
  err_msg <- paste0("ERROR: Invalid mtc_method argument. Please see", 
                    " '?p.adjust' for valid methods.")
  if (!mtc_method %in% c(stats::p.adjust.methods)) {
    stop(err_msg)
  }
}

ewce_barplot <- function (total_res, mtc_method = "bonferroni", q_threshold = 0.05, 
                          ctd = NULL, annotLevel = 1, heights = c(0.3, 1), make_dendro = FALSE, 
                          verbose = TRUE) 
{
  requireNamespace("ggplot2")
  requireNamespace("patchwork")
  check_mtc_method(mtc_method = mtc_method)
  multiList <- TRUE
  if (is.null(total_res$list)) 
    multiList <- FALSE
  if (isTRUE(make_dendro)) {
    if (is.null(ctd)) {
      messager("Warning: Can only add the dendrogram when ctd is provided.", 
               "Setting make_dendro=FALSE.", v = verbose)
      make_dendro <- FALSE
    }
    else {
      if (length(ctd[[annotLevel]]$plotting) > 0) {
        annotLevel <- which(unlist(lapply(ctd, FUN = cells_in_ctd, 
                                          cells = as.character(total_res$CellType))) == 
                              1)
        err_msg2 <- paste0("All of the cells within total_res should come", 
                           " from a single annotation layer of the CTD")
        if (length(annotLevel) == 0) {
          stop(err_msg2)
        }
      }
      if (length(ctd[[annotLevel]]$plotting) > 0) {
        total_res$CellType <- factor(x = fix_celltype_names(total_res$CellType), 
                                     levels = fix_celltype_names(ctd[[annotLevel]]$plotting$cell_ordering), 
                                     ordered = TRUE)
      }
    }
  }
  if (!"q" %in% colnames(total_res)) {
    total_res$q <- stats::p.adjust(total_res$p, method = mtc_method)
  }
  ast_q <- rep("", dim(total_res)[1])
  ast_q[total_res$q < q_threshold] <- "*"
  total_res$ast_q <- ast_q
  total_res$sd_from_mean[total_res$sd_from_mean < 0] <- 0
  graph_theme <- ggplot2::theme_bw(base_size = 12, base_family = "Helvetica") + 
    ggplot2::theme(text = ggplot2::element_text(size = 14), 
                   axis.title.y = ggplot2::element_text(vjust = 0.6), 
                   strip.background = ggplot2::element_rect(fill = "white"), 
                   strip.text = ggplot2::element_text(color = "black"))
  upperLim <- max(abs(total_res$sd_from_mean), na.rm = TRUE)
  total_res$y_ast <- total_res$sd_from_mean * 1.05
  total_res$abs_sd <- abs(total_res$sd_from_mean)
  if ("Direction" %in% colnames(total_res)) {
    the_plot <- ggplot2::ggplot(total_res) + ggplot2::geom_bar(ggplot2::aes_string(x = "CellType", 
                                                                                   y = "abs_sd", fill = "Direction"), position = "dodge", 
                                                               stat = "identity") + graph_theme
  }
  else {
    the_plot <- ggplot2::ggplot(total_res) + ggplot2::geom_bar(ggplot2::aes_string(x = "CellType", 
                                                                                   y = "abs_sd", fill = "abs_sd"), stat = "identity") + 
      ggplot2::scale_fill_gradientn(colors = c(colorRampPalette(colors = c("#2166AC","#4393C3","#92C5DE","#D1E5F0"))(50),
                                               colorRampPalette(colors = c("#FDDBC7","#F4A582","#D6604D","#B2182B"))(50))) + 
      graph_theme + ggplot2::theme(legend.position = "none")
  }
  the_plot <- the_plot + ggplot2::theme(plot.margin = ggplot2::unit(c(0.5, 
                                                                      0, 0, 0), "mm"), axis.text.x = ggplot2::element_text(angle = 55, 
                                                                                                                           hjust = 1)) + ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", 
                                                                                                                                                                                             fill = NA, linewidth = 1)) + ggplot2::xlab("Cell type") + 
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0)) + 
    ggplot2::ylab("Std.Devs. from the mean")
  the_plot <- the_plot + ggplot2::scale_y_continuous(breaks = c(0, 
                                                                ceiling(upperLim * 0.66)), expand = c(0, 1.1)) + ggplot2::geom_text(ggplot2::aes_string(label = "ast_q", 
                                                                                                                                                        x = "CellType", y = "y_ast"), size = 10)
  if (isTRUE(multiList)) {
    the_plot <- the_plot + ggplot2::facet_grid("list ~ .", 
                                               scales = "free", space = "free_x")
  }
  output <- list()
  output$plain <- the_plot
  if (isTRUE(make_dendro)) {
    ctdIN <- prep_dendro(ctdIN = ctd[[annotLevel]], expand = c(0, 
                                                               0.66))
    output$withDendro <- patchwork::wrap_plots(ctdIN$plotting$ggdendro_horizontal, 
                                               the_plot, heights = heights, ncol = 1)
  }
  return(output)
}

#### catdotplot ----
catdotplot <- function(adata,
                       x,
                       y,
                       color,
                       size,
                       title = NULL,
                       annotation_row = NULL,
                       annotation_col = NULL,
                       cluster_row = FALSE,
                       cluster_col = FALSE,
                       dot_scale = 6) {
  x <- enquo(x)
  y <- enquo(y)
  size <- enquo(size)
  color <- enquo(color)
  
  if (!is.null(annotation_col)) {
    x_order <- annotation_col %>% pull(1)
    adata <- adata %>%
      mutate(x = as_factor(x)) %>%
      mutate(x = fct_relevel(x, x_order))
  }
  if (!is.null(annotation_row)) {
    y_order <- annotation_row %>% pull(1)
    adata <- adata %>%
      mutate(y = as_factor(y)) %>%
      mutate(y = fct_relevel(y, y_order))
  }
  p <- adata %>%
    ggplot(aes(
      x = {{ x }},
      y = {{ y }},
      size = {{ size }},
      color = {{ color }}
    )) +
    geom_point() +
    theme_cat() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        colour = "black",
        size = 6,
        # face = "italic"
      ),
      axis.text.y = element_text(
        colour = "black",
        # face = "italic",
        size = 6
      ),
      axis.title = element_blank(),
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      legend.key.size = unit(6, "pt"),
      legend.key = element_blank(),
      legend.margin = margin(r = 0, unit = "pt"),
      plot.title = element_text(
        size = 6,
        colour = "black",
        face = "plain",
        hjust = 0.5
      )
    ) +
    scale_color_gradient2(low = "#00a6e1", mid = "white", high = "#ee6470") +
    scale_radius(range = c(0, dot_scale))
  if (!is.null(title)) {
    p <- p + labs(title = title)
  }
  if (!is.null(annotation_col)) {
    pc <- catlabel(annotation_col) +
      scale_fill_manual(values = paletteer::paletteer_d("RColorBrewer::Set3"))
    
    p <- p %>%
      aplot::insert_top(pc, height = 0.02)
  }
  if (!is.null(annotation_row)) {
    pr <- catlabel(annotation_row, flip = T)
    # scale_fill_manual(values = paletteer::paletteer_d("RColorBrewer::Paired"))
    p <- p %>%
      aplot::insert_right(pr, width = 0.02)
  }
  if (cluster_row) {
    cluster_row <- adata %>%
      ungroup() %>%
      select(x, y, {{ color }}) %>%
      pivot_wider(names_from = x, values_from = {{ color }}) %>%
      replace(is.na(.), 0) %>%
      column_to_rownames("y") %>%
      dist() %>%
      hclust() %>%
      ggtree::ggtree()
    p <- p %>%
      aplot::insert_left(cluster_row, width = 0)
  }
  if (cluster_col) {
    cluster_col <- adata %>%
      ungroup() %>%
      select(x, y, {{ color }}) %>%
      pivot_wider(names_from = y, values_from = {{ color }}) %>%
      replace(is.na(.), 0) %>%
      column_to_rownames("x") %>%
      dist() %>%
      hclust() %>%
      ggtree::ggtree() +
      ggtree::layout_dendrogram()
    p <- p %>%
      aplot::insert_top(cluster_col, height = 0)
  }
  return(p)
}

#### theme_cat ----
theme_cat <- function(...) {
  theme(
    # Axis
    axis.title = element_text(
      size = 8,
      face = "plain",
      colour = "black"
    ),
    axis.text = element_text(
      size = 8,
      face = "plain",
      colour = "black"
    ),
    axis.line = element_blank(),
    axis.ticks = element_line(
      size = 0.25,
      lineend = "square",
      colour = "black"
    ),
    # Panel
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(
      size = 0.25,
      fill = NA,
      colour = "black"
    ),
    # Legend
    legend.background = element_blank(),
    legend.title = element_text(
      size = 8,
      face = "plain",
      colour = "black"
    ),
    legend.text = element_text(
      size = 8,
      face = "plain",
      colour = "black"
    ),
    legend.key = element_blank(),
    legend.key.height = unit(8, "pt"),
    legend.key.width = unit(8, "pt"),
    # Plot
    plot.background = element_blank(),
    plot.title = element_text(
      size = 8,
      hjust = 0.5,
      face = "plain",
      colour = "black"
    ),
    # plot.margin = margin(
    #   t = 0,
    #   r = 0,
    #   b = 0,
    #   l = 0
    # ),
    # Facetting
    strip.background = element_blank(),
    strip.text = element_text(
      size = 8,
      face = "plain",
      colour = "black"
    ),
    ...
  )
}

#### pheatmap_rowname_left ----
library(pheatmap)
library(grid)
library(gtable)

# Modified pheatmap:::heatmap_motor
heatmap_motor <- function (matrix, border_color, cellwidth, cellheight, tree_col, 
                           tree_row, treeheight_col, treeheight_row, filename, width, 
                           height, breaks, color, legend, annotation_row, annotation_col, 
                           annotation_colors, annotation_legend, annotation_names_row, 
                           annotation_names_col, main, fontsize, fontsize_row, fontsize_col, 
                           hjust_col, vjust_col, angle_col, fmat, fontsize_number, number_color, 
                           gaps_col, gaps_row, labels_row, labels_col, ...) 
{
  lo = pheatmap:::lo(coln = labels_col, rown = labels_row, nrow = nrow(matrix), 
                     ncol = ncol(matrix), cellwidth = cellwidth, cellheight = cellheight, 
                     treeheight_col = treeheight_col, treeheight_row = treeheight_row, 
                     legend = legend, annotation_col = annotation_col, annotation_row = annotation_row, 
                     annotation_colors = annotation_colors, annotation_legend = annotation_legend, 
                     annotation_names_row = annotation_names_row, annotation_names_col = annotation_names_col, 
                     main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
                     fontsize_col = fontsize_col, angle_col = angle_col, gaps_row = gaps_row, 
                     gaps_col = gaps_col, ...)
  res = lo$gt
  mindim = lo$mindim
  if (!is.na(filename)) {
    if (is.na(height)) {
      height = convertHeight(gtable_height(res), "inches", valueOnly = T)
    }
    if (is.na(width)) {
      width = convertWidth(gtable_width(res), "inches", valueOnly = T)
    }
    r = regexpr("\\.[a-zA-Z]*$", filename)
    if (r == -1) 
      stop("Improper filename")
    ending = substr(filename, r + 1, r + attr(r, "match.length"))
    f = switch(ending, pdf = function(x, ...) pdf(x, ...), 
               png = function(x, ...) png(x, units = "in", res = 300, 
                                          ...), jpeg = function(x, ...) jpeg(x, units = "in", 
                                                                             res = 300, ...), jpg = function(x, ...) jpeg(x, 
                                                                                                                          units = "in", res = 300, ...), tiff = function(x, 
                                                                                                                                                                         ...) tiff(x, units = "in", res = 300, compression = "lzw", 
                                                                                                                                                                                   ...), bmp = function(x, ...) bmp(x, units = "in", 
                                                                                                                                                                                                                    res = 300, ...), stop("File type should be: pdf, png, bmp, jpg, tiff"))
    f(filename, height = height, width = width)
    gt = heatmap_motor(matrix, cellwidth = cellwidth, cellheight = cellheight, 
                       border_color = border_color, tree_col = tree_col, 
                       tree_row = tree_row, treeheight_col = treeheight_col, 
                       treeheight_row = treeheight_row, breaks = breaks, 
                       color = color, legend = legend, annotation_col = annotation_col, 
                       annotation_row = annotation_row, annotation_colors = annotation_colors, 
                       annotation_legend = annotation_legend, annotation_names_row = annotation_names_row, 
                       annotation_names_col = annotation_names_col, filename = NA, 
                       main = main, fontsize = fontsize, fontsize_row = fontsize_row, 
                       fontsize_col = fontsize_col, hjust_col = hjust_col, 
                       vjust_col = vjust_col, angle_col = angle_col, fmat = fmat, 
                       fontsize_number = fontsize_number, number_color = number_color, 
                       labels_row = labels_row, labels_col = labels_col, 
                       gaps_col = gaps_col, gaps_row = gaps_row, ...)
    grid.draw(gt)
    dev.off()
    return(gt)
  }
  if (mindim < 3) 
    border_color = NA
  if (!is.na(main)) {
    elem = pheatmap:::draw_main(main, fontsize = 1.3 * fontsize, ...)
    res = gtable_add_grob(res, elem, t = 1, l = 3, name = "main", 
                          clip = "off")
  }
  if (!pheatmap:::is.na2(tree_col) & treeheight_col != 0) {
    elem = pheatmap:::draw_dendrogram(tree_col, gaps_col, horizontal = T)
    res = gtable_add_grob(res, elem, t = 2, l = 3, name = "col_tree")
  }
  if (!pheatmap:::is.na2(tree_row) & treeheight_row != 0) {
    elem = pheatmap:::draw_dendrogram(tree_row, gaps_row, horizontal = F)
    res = gtable_add_grob(res, elem, t = 4, l = 1, name = "row_tree")
  }
  elem = pheatmap:::draw_matrix(matrix, border_color, gaps_row, gaps_col, 
                                fmat, fontsize_number, number_color)
  res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", 
                        name = "matrix")
  if (length(labels_col) != 0) {
    pars = list(labels_col, gaps = gaps_col, fontsize = fontsize_col, 
                hjust_col = hjust_col, vjust_col = vjust_col, angle_col = angle_col, 
                ...)
    elem = do.call(pheatmap:::draw_colnames, pars)
    res = gtable_add_grob(res, elem, t = 5, l = 3, clip = "off", 
                          name = "col_names")
  }
  if (length(labels_row) != 0) {
    pars = list(labels_row, gaps = gaps_row, fontsize = fontsize_row, 
                ...)
    elem = do.call(pheatmap:::draw_rownames, pars)
    res = gtable_add_grob(res, elem, t = 4, l = 3, clip = "off", 
                          name = "row_names")
  }
  if (!pheatmap:::is.na2(annotation_col)) {
    converted_annotation = convert_annotations(annotation_col, 
                                               annotation_colors)
    elem = pheatmap:::draw_annotations(converted_annotation, border_color, 
                                       gaps_col, fontsize, horizontal = T)
    res = gtable_add_grob(res, elem, t = 3, l = 3, clip = "off", 
                          name = "col_annotation")
    if (annotation_names_col) {
      elem = pheatmap:::draw_annotation_names(annotation_col, fontsize, 
                                              horizontal = T)
      res = gtable_add_grob(res, elem, t = 3, l = 4, clip = "off", 
                            name = "col_annotation_names")
    }
  }
  if (!pheatmap:::is.na2(annotation_row)) {
    converted_annotation = convert_annotations(annotation_row, 
                                               annotation_colors)
    elem = pheatmap:::draw_annotations(converted_annotation, border_color, 
                                       gaps_row, fontsize, horizontal = F)
    res = gtable_add_grob(res, elem, t = 4, l = 2, clip = "off", 
                          name = "row_annotation")
    if (annotation_names_row) {
      elem = pheatmap:::draw_annotation_names(annotation_row, fontsize, 
                                              horizontal = F, hjust_col = hjust_col, vjust_col = vjust_col, 
                                              angle_col = angle_col)
      res = gtable_add_grob(res, elem, t = 5, l = 2, clip = "off", 
                            name = "row_annotation_names")
    }
  }
  annotation = c(annotation_col[length(annotation_col):1], 
                 annotation_row[length(annotation_row):1])
  annotation = annotation[unlist(lapply(annotation, function(x) !pheatmap:::is.na2(x)))]
  if (length(annotation) > 0 & annotation_legend) {
    elem = pheatmap:::draw_annotation_legend(annotation, annotation_colors, 
                                             border_color, fontsize = fontsize, ...)
    t = ifelse(is.null(labels_row), 4, 3)
    res = gtable_add_grob(res, elem, t = t, l = 6, b = 5, 
                          clip = "off", name = "annotation_legend")
  }
  if (!pheatmap:::is.na2(legend)) {
    elem = pheatmap:::draw_legend(color, breaks, legend, fontsize = fontsize, 
                                  ...)
    t = ifelse(is.null(labels_row), 4, 3)
    res = gtable_add_grob(res, elem, t = t, l = 5, b = 5, 
                          clip = "off", name = "legend")
  }
  return(res)
}

# Modified pheatmap:::lo    
lo <- function (rown, coln, nrow, ncol, cellheight = NA, cellwidth = NA, 
                treeheight_col, treeheight_row, legend, annotation_row, annotation_col, 
                annotation_colors, annotation_legend, annotation_names_row, 
                annotation_names_col, main, fontsize, fontsize_row, fontsize_col, 
                angle_col, gaps_row, gaps_col, ...) 
{
  if (!is.null(coln[1]) | (!pheatmap:::is.na2(annotation_row) & annotation_names_row)) {
    if (!is.null(coln[1])) {
      t = coln
    }
    else {
      t = ""
    }
    tw = strwidth(t, units = "in", cex = fontsize_col/fontsize)
    if (annotation_names_row) {
      t = c(t, colnames(annotation_row))
      tw = c(tw, strwidth(colnames(annotation_row), units = "in"))
    }
    longest_coln = which.max(tw)
    gp = list(fontsize = ifelse(longest_coln <= length(coln), 
                                fontsize_col, fontsize), ...)
    coln_height = unit(1, "grobheight", textGrob(t[longest_coln], 
                                                 rot = angle_col, gp = do.call(gpar, gp))) + unit(10, 
                                                                                                  "bigpts")
  }
  else {
    coln_height = unit(5, "bigpts")
  }
  if (!is.null(rown[1])) {
    t = rown
    tw = strwidth(t, units = "in", cex = fontsize_row/fontsize)
    if (annotation_names_col) {
      t = c(t, colnames(annotation_col))
      tw = c(tw, strwidth(colnames(annotation_col), units = "in"))
    }
    longest_rown = which.max(tw)
    gp = list(fontsize = ifelse(longest_rown <= length(rown), 
                                fontsize_row, fontsize), ...)
    rown_width = unit(1, "grobwidth", textGrob(t[longest_rown], 
                                               rot = 0, gp = do.call(gpar, gp))) + unit(10, "bigpts")
  }
  else {
    rown_width = unit(5, "bigpts")
  }
  gp = list(fontsize = fontsize, ...)
  if (!pheatmap:::is.na2(legend)) {
    longest_break = which.max(nchar(names(legend)))
    longest_break = unit(1.1, "grobwidth", 
                         textGrob(as.character(names(legend))[longest_break], 
                                  gp = do.call(gpar, gp)))
    title_length = unit(1.1, "grobwidth", textGrob("Scale", 
                                                   gp = gpar(fontface = "bold", ...)))
    legend_width = unit(12, "bigpts") + longest_break * 1.2
    legend_width = max(title_length, legend_width)
  }
  else {
    legend_width = unit(0, "bigpts")
  }
  if (is.na(main)) {
    main_height = unit(0, "npc")
  }
  else {
    main_height = unit(1.5, "grobheight", textGrob(main, 
                                                   gp = gpar(fontsize = 1.3 * fontsize, ...)))
  }
  textheight = unit(fontsize, "bigpts")
  if (!pheatmap:::is.na2(annotation_col)) {
    annot_col_height = ncol(annotation_col) * (textheight + 
                                                 unit(2, "bigpts")) + unit(2, "bigpts")
    t = c(as.vector(as.matrix(annotation_col)), colnames(annotation_col))
    annot_col_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], 
                                                             gp = gpar(...))) + unit(12, "bigpts")
    if (!annotation_legend) {
      annot_col_legend_width = unit(0, "npc")
    }
  }
  else {
    annot_col_height = unit(0, "bigpts")
    annot_col_legend_width = unit(0, "bigpts")
  }
  if (!pheatmap:::is.na2(annotation_row)) {
    annot_row_width = ncol(annotation_row) * (textheight + 
                                                unit(2, "bigpts")) + unit(2, "bigpts")
    t = c(as.vector(as.matrix(annotation_row)), colnames(annotation_row))
    annot_row_legend_width = unit(1.2, "grobwidth", textGrob(t[which.max(nchar(t))], 
                                                             gp = gpar(...))) + unit(12, "bigpts")
    if (!annotation_legend) {
      annot_row_legend_width = unit(0, "npc")
    }
  }
  else {
    annot_row_width = unit(0, "bigpts")
    annot_row_legend_width = unit(0, "bigpts")
  }
  annot_legend_width = max(annot_row_legend_width, annot_col_legend_width)
  treeheight_col = unit(treeheight_col, "bigpts") + unit(5, 
                                                         "bigpts")
  treeheight_row = unit(treeheight_row, "bigpts") + unit(5, 
                                                         "bigpts")
  if (is.na(cellwidth)) {
    mat_width = unit(1, "npc") - rown_width - legend_width - 
      treeheight_row - annot_row_width - annot_legend_width
  }
  else {
    mat_width = unit(cellwidth * ncol, "bigpts") + length(gaps_col) * 
      unit(4, "bigpts")
  }
  if (is.na(cellheight)) {
    mat_height = unit(1, "npc") - main_height - coln_height - 
      treeheight_col - annot_col_height
  }
  else {
    mat_height = unit(cellheight * nrow, "bigpts") + length(gaps_row) * 
      unit(4, "bigpts")
  }
  gt = gtable(widths = unit.c(treeheight_row, rown_width,  
                              mat_width, treeheight_row, legend_width, annot_legend_width), 
              heights = unit.c(main_height, treeheight_col, annot_col_height, 
                               mat_height, coln_height), vp = viewport(gp = do.call(gpar, 
                                                                                    gp)))
  cw = convertWidth(mat_width - (length(gaps_col) * unit(4, 
                                                         "bigpts")), "bigpts", valueOnly = T)/ncol
  ch = convertHeight(mat_height - (length(gaps_row) * unit(4, 
                                                           "bigpts")), "bigpts", valueOnly = T)/nrow
  mindim = min(cw, ch)
  res = list(gt = gt, mindim = mindim)
  return(res)
}

# Modified pheatmap:::draw_rownames      
draw_rownames <- function (rown, gaps, ...) 
{
  coord = pheatmap:::find_coordinates(length(rown), gaps)
  y = unit(1, "npc") - (coord$coord - 0.5 * coord$size)
  res = textGrob(rown, x = unit(-3, "bigpts"), y = y, vjust = 0.5, 
                 hjust = 1, gp = gpar(...))
  return(res)
}

assignInNamespace(x="draw_rownames", value=draw_rownames, ns="pheatmap")
assignInNamespace(x="lo", value=lo, ns="pheatmap")
assignInNamespace(x="heatmap_motor", value=heatmap_motor, ns="pheatmap") 

#### catvolcano ----
library(tidyverse)

catvolcano <-
  function(adata,
           x,
           y,
           label,
           text = NULL,
           p_value = 0.05,
           log2FC = 0.25) {
    adata <- adata %>%
      dplyr::rename(x = {{ x }}, y = {{ y }}, label = {{ label }})
    adata <- adata %>%
      mutate(change = ifelse(
        y <= p_value & abs(x) >= log2FC,
        ifelse(x >= log2FC, "Upregulated", "Downregulated"),
        "Stable"
      ))
    max_x <- ceiling(max(abs(adata$x)))
    if (!is.null(text)) {
      adata <- adata %>% mutate(label = if_else(label %in% text, label, ""))
    }
    print(adata %>% count(change))
    # print(adata[adata$label != "",])
    p <- adata %>%
      ggplot(aes(
        x = x,
        y = -log10(y),
        colour = change
      )) +
      ggrastr::geom_point_rast(raster.dpi = 600, size = 0.5) +
      geom_vline(
        xintercept = c(-log2FC, log2FC),
        lty = 2,
        col = "black",
        lwd = 0.5
      ) +
      geom_hline(
        yintercept = -log10(p_value),
        lty = 2,
        col = "black",
        lwd = 0.5
      ) +
      labs(
        x = "log2(Fold change)",
        y = "-log10 (P value)"
      ) +
      scale_color_manual(
        values = c("#00a6e1", "lightgrey", "#ee6470"),
        labels = c(
          str_c(
            "Down (",
            adata %>% count(change) %>% filter(change == "Downregulated") %>% pull(n),
            ")"
          ),
          "Stable",
          str_c(
            "Up (",
            adata %>% count(change) %>% filter(change == "Upregulated") %>% pull(n),
            ")"
          )
        )
      ) +
      theme_cat() +
      theme(
        aspect.ratio = 1,
        legend.position = "top",
        legend.title = element_blank(),
        legend.margin = margin(b = -8)
      ) +
      xlim(c(-max_x, max_x))
    # theme(
    #   plot.title = element_blank(),
    #   legend.position = "top",
    #   legend.title = element_blank(),
    #   legend.text = element_text(size = 6, colour = "black"),
    #   legend.key = element_blank(),
    #   aspect.ratio = 1,
    #   panel.background = element_blank(),
    #   panel.grid = element_blank(),
    #   axis.text = element_text(size = 6, colour = "black"),
    #   axis.ticks = element_line(size = 0.5, colour = "black"),
    #   axis.line = element_blank(),
    #   axis.title = element_text(size = 6, colour = "black"),
    #   panel.border = element_rect(
    #     size = 0.5,
    #     colour = "black",
    #     fill = F
    #   )
    # )
    
    if (!is.null(text)) {
      p <- p + ggrepel::geom_text_repel(
        aes(label = label),
        # max.overlaps = 1000000,
        max.overlaps = Inf,
        color = "black", 
        size = 3,
        fontface = "italic"
      )
    }
    
    return(p)
  }

#### GSEA ----
tableGrob2 <- function(d, p = NULL) {
  d <- d[order(rownames(d)), ]
  tp <- gridExtra::tableGrob(d,
                             theme = gridExtra::ttheme_minimal(base_size = 6)
  )
  if (is.null(p)) {
    return(tp)
  }
  p_data <- ggplot_build(p)$data[[1]]
  p_data <- p_data[order(p_data[["group"]]), ]
  pcol <- unique(p_data[["colour"]])
  j <- which(tp$layout$name == "rowhead-fg")
  for (i in seq_along(pcol)) {
    tp$grobs[j][[i + 1]][["gp"]] <- grid::gpar(col = pcol[i])
  }
  return(tp)
}
cat_gseaplot <-
  function(x,
           geneSetID,
           title = "",
           color = "#8ac45f",
           base_size = 6,
           rel_heights = c(1.6, 0.4, 1),
           subplots = 1:3,
           pvalue_table = FALSE,
           ES_geom = "line") {
    ES_geom <- match.arg(ES_geom, c("line", "dot"))
    geneList <- position <- NULL
    if (length(geneSetID) == 1) {
      gsdata <- enrichplot:::gsInfo(x, geneSetID)
    } else {
      gsdata <-
        do.call(rbind, lapply(geneSetID, enrichplot:::gsInfo, object = x))
    }
    p <- ggplot(gsdata, aes_(x = ~x)) +
      xlab(NULL) +
      theme_cat() +
      theme(axis.line.x.bottom = element_blank()) +
      scale_x_continuous(expand = c(0, 0))
    if (ES_geom == "line") {
      es_layer <- geom_line(aes_(y = ~runningScore, color = ~Description),
                            size = 0.5
      )
    } else {
      es_layer <-
        geom_point(
          aes_(y = ~runningScore, color = ~Description),
          size = 0.5,
          data = subset(gsdata, position == 1)
        )
    }
    p.res <- p + es_layer + theme(
      legend.position = c(0.8, 0.8),
      legend.title = element_blank(),
      legend.background = element_rect(fill = "transparent")
    )
    p.res <- p.res + ylab("Enrichment Score") + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      plot.margin = margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0,
        unit = "cm"
      )
    )
    i <- 0
    for (term in unique(gsdata$Description)) {
      idx <- which(gsdata$ymin != 0 & gsdata$Description ==
                     term)
      gsdata[idx, "ymin"] <- i
      gsdata[idx, "ymax"] <- i + 1
      i <- i + 1
    }
    p2 <- ggplot(gsdata, aes_(x = ~x)) +
      geom_linerange(aes_(
        ymin = ~ymin,
        ymax = ~ymax,
        color = ~Description
      )) +
      theme_classic(base_size) +
      theme(
        legend.position = "none",
        plot.margin = margin(t = 0, b = 0),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line.x = element_blank(),
        axis.line.x.bottom = element_blank(),
        axis.line.y.left = element_blank(),
        panel.border = element_rect(
          size = 0.25,
          fill = NA,
          colour = "black"
        ),
        plot.title = element_blank(),
        axis.title = element_blank()
      ) +
      scale_x_continuous(expand = c(
        0,
        0
      )) +
      scale_y_continuous(expand = c(0, 0))
    if (length(geneSetID) == 1) {
      v <- seq(1, sum(gsdata$position), length.out = 9)
      inv <- findInterval(rev(cumsum(gsdata$position)), v)
      if (min(inv) == 0) {
        inv <- inv + 1
      }
      col <-
        c(
          rev(RColorBrewer::brewer.pal(5, "Blues")),
          RColorBrewer::brewer.pal(
            5,
            "Reds"
          )
        )
      ymin <- min(p2$data$ymin)
      yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
      xmin <- which(!duplicated(inv))
      xmax <-
        xmin + as.numeric(table(inv)[as.character(unique(inv))])
      d <- data.frame(
        ymin = ymin,
        ymax = yy,
        xmin = xmin,
        xmax = xmax,
        col = col[unique(inv)]
      )
      p2 <- p2 + geom_rect(
        aes_(
          xmin = ~xmin,
          xmax = ~xmax,
          ymin = ~ymin,
          ymax = ~ymax,
          fill = ~ I(col)
        ),
        data = d,
        alpha = 0.9,
        inherit.aes = FALSE
      )
    }
    df2 <- p$data
    df2$y <- p$data$geneList[df2$x]
    p.pos <- p + geom_segment(
      data = df2,
      aes_(
        x = ~x,
        xend = ~x,
        y = ~y,
        yend = 0
      ),
      color = "grey"
    )
    p.pos <-
      p.pos + ylab("Ranked List Metric") + xlab("Rank in Ordered Dataset") +
      theme(plot.margin = margin(
        t = -0.1,
        r = 0.2,
        b = 0.2,
        l = 0.2,
        unit = "cm"
      ))
    if (!is.null(title) && !is.na(title) && title != "") {
      p.res <- p.res + ggtitle(title)
    }
    if (length(color) == length(geneSetID)) {
      p.res <- p.res + scale_color_manual(values = color)
      if (length(color) == 1) {
        p.res <- p.res + theme(legend.position = "none")
        p2 <- p2 + scale_color_manual(values = "black")
      } else {
        p2 <- p2 + scale_color_manual(values = color)
      }
    }
    if (pvalue_table) {
      # pd <- x[geneSetID, c("Description", "NES", "qvalues")]
      pd <- x[geneSetID, c("Description", "NES", "p.adjust")]
      colnames(pd)[3] <- c("FDR")
      pd <- pd[, -1]
      # pd <- round(pd, 2)
      # print(pd$FDR)
      pd$NES <- round(pd$NES, 2)
      pd$FDR <- format(pd$FDR, scientific = T, digits = 2)
      if (length(geneSetID) == 1) {
        rownames(pd) <- ""
        # pd <- as.data.frame(t(pd))
        # pd$V2 <- c("NES", "FDR")
        # pd <- pd[,c("V2", "V1")]
      }
      print(pd)
      tp <- tableGrob2(pd, p.res)
      p.res <-
        p.res + theme(legend.position = "none") + annotation_custom(
          tp,
          xmin = quantile(p.res$data$x, 0.75),
          xmax = quantile(
            p.res$data$x,
            0.8
          ),
          ymin = quantile(
            p.res$data$runningScore,
            0.6
          ),
          ymax = quantile(
            p.res$data$runningScore,
            0.8
          )
        )
    }
    plotlist <- list(p.res, p2, p.pos)[subplots]
    n <- length(plotlist)
    plotlist[[n]] <- plotlist[[n]] + theme(
      axis.line.x = element_line(),
      axis.ticks.x = element_line(),
      axis.text.x = element_text()
    )
    if (length(subplots) == 1) {
      # return(plotlist[[1]] + theme(plot.margin = margin(
      #   t = 0.2,
      #   r = 0.2, b = 0.2, l = 0.2, unit = "cm"
      # )))
      return(plotlist[[1]])
    }
    if (length(rel_heights) > length(subplots)) {
      rel_heights <- rel_heights[subplots]
    }
    aplot::plot_list(
      gglist = plotlist,
      ncol = 1,
      heights = rel_heights
    )
  }