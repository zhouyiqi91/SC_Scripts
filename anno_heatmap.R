DoHeatmap<- function (object, data.use = NULL, use.scaled = TRUE, cells.use = NULL, 
          genes.use = NULL, disp.min = -2.5, disp.max = 2.5, group.by = "ident", 
          group.order = NULL, draw.line = TRUE, col.low = "#FF00FF", 
          col.mid = "#000000", col.high = "#FFFF00", slim.col.label = FALSE, 
          remove.key = FALSE, rotate.key = FALSE, title = NULL, cex.col = 10, 
          cex.row = 10, group.label.loc = "bottom", group.label.rot = FALSE, 
          group.cex = 15, group.spacing = 0.15, assay.type = "RNA", 
          do.plot = TRUE,row_annotation=NULL) 
{
  if (is.null(x = data.use)) {
    if (use.scaled) {
      data.use <- GetAssayData(object, assay.type = assay.type, 
                               slot = "scale.data")
    }
    else {
      data.use <- GetAssayData(object, assay.type = assay.type, 
                               slot = "data")
    }
  }
  cells.use <- SetIfNull(x = cells.use, default = object@cell.names)
  cells.use <- intersect(x = cells.use, y = colnames(x = data.use))
  if (length(x = cells.use) == 0) {
    stop("No cells given to cells.use present in object")
  }
  genes.use <- SetIfNull(x = genes.use, default = rownames(x = data.use))
  genes.use <- intersect(x = genes.use, y = rownames(x = data.use))
  if (length(x = genes.use) == 0) {
    stop("No genes given to genes.use present in object")
  }
  if (is.null(x = group.by) || group.by == "ident") {
    cells.ident <- object@ident[cells.use]
  }
  else {
    cells.ident <- factor(x = FetchData(object = object, 
                                        cells.use = cells.use, vars.all = group.by)[, 1])
    names(x = cells.ident) <- cells.use
  }
  cells.ident <- factor(x = cells.ident, labels = intersect(x = levels(x = cells.ident), 
                                                            y = cells.ident))
  data.use <- data.use[genes.use, cells.use, drop = FALSE]
  if ((!use.scaled)) {
    data.use = as.matrix(x = data.use)
    if (disp.max == 2.5) 
      disp.max = 10
  }
  data.use <- MinMax(data = data.use, min = disp.min, max = disp.max)
  data.use <- as.data.frame(x = t(x = data.use))
  data.use$cell <- rownames(x = data.use)
  colnames(x = data.use) <- make.unique(names = colnames(x = data.use))
  data.use <- data.use %>% melt(id.vars = "cell")
  names(x = data.use)[names(x = data.use) == "variable"] <- "gene"
  names(x = data.use)[names(x = data.use) == "value"] <- "expression"
  data.use$ident <- cells.ident[data.use$cell]
  data.use$row_anno <- row_annotation[data.use$gene]
  if (!is.null(group.order)) {
    if (length(group.order) == length(levels(data.use$ident)) && 
        all(group.order %in% levels(data.use$ident))) {
      data.use$ident <- factor(data.use$ident, levels = group.order)
    }
    else {
      stop("Invalid group.order")
    }
  }
  breaks <- seq(from = min(data.use$expression), to = max(data.use$expression), 
                length = length(x = PurpleAndYellow()) + 1)
  data.use$gene <- with(data = data.use, expr = factor(x = gene, 
                                                       levels = rev(x = unique(x = data.use$gene))))
  data.use$cell <- with(data = data.use, expr = factor(x = cell, 
                                                       levels = cells.use))
  if (rotate.key) {
    key.direction <- "horizontal"
    key.title.pos <- "top"
  }
  else {
    key.direction <- "vertical"
    key.title.pos <- "left"
  }
  heatmap <- ggplot(data = data.use, mapping = aes(x = cell, 
                                                   y = gene, fill = expression)) + geom_tile() + scale_fill_gradient2(low = col.low, 
                                                                                                                      mid = col.mid, high = col.high, name = "Expression", 
                                                                                                                      guide = guide_colorbar(direction = key.direction, title.position = key.title.pos)) + 
    scale_y_discrete(position = "right") + 
    theme(axis.line = element_blank(), axis.title.y = element_blank(), 
          axis.ticks.y = element_blank(), strip.text.x = element_text(size = group.cex), 
          axis.text.y = element_text(size = cex.row), axis.text.x = element_text(size = cex.col), 
          axis.title.x = element_blank())
  if (slim.col.label) {
    heatmap <- heatmap + theme(axis.title.x = element_blank(), 
                               axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
                               axis.line = element_blank(), axis.title.y = element_blank(), 
                               axis.ticks.y = element_blank())
  }
  else {
    heatmap <- heatmap + theme(axis.text.x = element_text(angle = 90))
  }
  if (!is.null(x = group.by)) {
    if (group.label.loc == "top") {
      switch <- NULL
    }
    else {
      switch <- "x"
    }
    heatmap <- heatmap + facet_grid(facets = row_anno~ident, drop = TRUE, 
                                    space = "free", scales = "free", switch = switch, 
    ) + scale_x_discrete(expand = c(0, 0), drop = TRUE) +theme(
      strip.background = element_blank(),
      strip.text.y = element_blank()
    )
    if (draw.line) {
      panel.spacing <- unit(x = group.spacing, units = "lines")
    }
    else {
      panel.spacing <- unit(x = 0, units = "lines")
    }
    heatmap <- heatmap + theme(strip.background = element_blank(), 
                               panel.spacing = panel.spacing)
    if (group.label.rot) {
      heatmap <- heatmap + theme(strip.text.x = element_text(angle = 90))
    }
  }
  if (remove.key) {
    heatmap <- heatmap + theme(legend.position = "none")
  }
  if (!is.null(x = title)) {
    heatmap <- heatmap + labs(title = title)
  }
  if (do.plot) {
    heatmap
  }
  dd <<- data.use
  return(heatmap)
}

SetIfNull <- function (x, default) 
{
  if (is.null(x = x)) {
    return(default)
  }
  else {
    return(x)
  }
}

main <- function(){
	rds <- readRDS("/SGRNJ/Database/test/tests/auto_assign_test/rds/he_30.PRO.rds")
	marker_test <- FindAllMarkers(rds,genes.use = rds@var.genes,max.cells.per.ident = 300)
	top_gene <- marker_test %>% group_by(cluster) %>% top_n(3, avg_logFC)
	genes <- top_gene$gene
	no_dup <- !duplicated(genes)
	no_dup_genes <- genes[no_dup]

	rds <- ScaleData(rds,genes.use = genes)
	anno <- top_gene$cluster[no_dup]
	names(anno) <- top_3$gene[no_dup]

	print(DoHeatmap(object=rds,genes.use=no_dup_genes,
                slim.col.label=TRUE,
                group.label.rot=FALSE, 
                title = "cluster top markers",group.spacing = 0.2,rotate.key = TRUE,
                
                col.low="blue",col.mid="white",col.high="red",row_annotation = anno))
}
