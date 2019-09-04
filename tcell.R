### pheatmap
library(Seurat)
library(tidyverse)
library(pheatmap)
library(argparser)

argv <- arg_parser('')
argv <- add_argument(argv,"--rds", help="T cell Seurat rds")
argv <- parse_args(argv)

all_data <- readRDS(argv$rds)
markers <- read_tsv("/SGRNJ/Database/sc_marker/tcell/t_cell_marker.txt")
genes <- markers$gene

## z-score
z <- all_data@data[genes,]
x <- z
for (i in 1:nrow(z)){
  x[i,] <- scale(z[i,])
}

data.all <- data.frame(row.names = genes)

for (j in levels(x = all_data@ident)) {
  temp.cells <- WhichCells(object = all_data, ident = j)
  data.use <- x[,temp.cells,drop=FALSE]
  data.temp <- apply(X = data.use, MARGIN = 1, FUN = mean)
  data.all <- cbind(data.all, data.temp)
  colnames(x = data.all)[ncol(x = data.all)] <- j
}




anno <- data.frame(type=markers$type)
rownames(anno) <- genes
anno$type <- factor(anno$type,levels<-unique(anno$type))

# exp
#exp <- AverageExpression(all_data)
#exp1 <- exp[markers$gene,]

pdf("tcell_markers.pdf",width=15)
print (pheatmap(data.all,cluster_row = FALSE, cluster_col = T,display_numbers =T，
	annotation_row = anno,anotation_names_row=FALSE,gaps_row = head(as.numeric(cumsum(table(anno$type))), -1)))
dev.off()
