#!/usr/bin/env Rscript
library(argparser)
library(Seurat)
library(tidyverse)

argv <- arg_parser('')
argv <- add_argument(argv,"--outdir", help="the output dir.",default=getwd())
argv <- add_argument(argv,"--rds",help="Seurat rds")
argv <- add_argument(argv,"--level",help="calculate correlation on single cell level or cluster level",default="cluster")
argv <- add_argument(argv,"--resolution", help="tSNE resolution",default=0.8)
argv <- parse_args(argv)

#read args

outdir <- argv$outdir
rds <- argv$rds
level <- argv$level
resolution <- argv$resolution
origin.cluster <- paste("res.",resolution,sep="")

all_data <- readRDS(rds)
lm22 <- read_tsv("/SGRNJ/Database/download/bulk/LM22/LM22.txt")

color2 <- c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple",
"DarkSlateGray","Gold","DarkGreen","DeepPink2","Red4",
"#4682B4","#FFDAB9","#708090","#836FFF","#CDC673",
"#CD9B1D","#FF6EB4","#CDB5CD","#008B8B","#43CD80",
"#6A5ACD","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62"
                      ,"#68228B","#CDB7B5","#CD853F","#6B8E23","#E6E6FA","#FFDAB9")

setwd(outdir)

mye <- all_data@data
genes <- lm22$`Gene symbol`
genes1 <- rownames(mye)
genes2 <- intersect(genes,genes1)

l <- filter(lm22,`Gene symbol` %in% genes2)

if (level=="sc"){
  m <- mye[genes2,]
  barcodes <- colnames(m)
  index <- 0
  result <- tibble(barcode=character(),ident=character(),cor=numeric())
  for (barcode in barcodes){
    index <- index + 1
    print (index)
    exp <- m[,barcode]
    #str(exp)
    cors <- c()
    for (cell in colnames(l)[-1]){
      exp1 <- l[,cell]
      cors[cell] <- cor(exp,exp1,method="spearman")
    }
    max_cor <- max(cors)
    max_ident <- names(cors[cors==max_cor])
    result <- add_row(result,barcode=barcode,ident=max_ident,cor=max_cor)
  }
  result1 <- filter(result,duplicated(barcode)==F)

  x <- all_data@meta.data
  x1 <- rownames_to_column(x,var="barcode")
  x2 <- left_join(x1,result1,by="barcode")
  x3 <- column_to_rownames(x2,var="barcode")
  all_data@meta.data <- x3

  all_data <- SetAllIdent(all_data,id="ident")

  pdf("./sc_level_22ident.pdf",width=12)
  print (TSNEPlot(all_data,colors.use = color2,do.return=T) )
  dev.off()

  y <- all_data@meta.data
  y$first_ident <- map_chr(y$ident,
        function (x) unlist(strsplit(x,split=" "))[1])
  all_data@meta.data <- y
  all_data <- SetAllIdent(all_data,id="first_ident")
  pdf("./sc_level_10ident.pdf")
  TSNEPlot(all_data,colors.use = color2)
  dev.off()

  saveRDS(all_data,"./sc_level.rds")

}

if (level=="cluster"){
  avg_exp <- AverageExpression(all_data)  
  genes1 <- rownames(avg_exp)
  genes2 <- intersect(genes,genes1)

  m <- avg_exp[genes2,]
  l <- filter(lm22,`Gene symbol` %in% genes2)


  clusters <- colnames(m)
  cells <- colnames(l)[-1]

  result <- tibble(cluster=character(),ident=character(),cor=numeric())

  for (cluster in clusters){
    exp <- m[,cluster]
    cors <- c()
    for (cell in cells){
      exp1 <- l[,cell]
      cors[cell] <- cor(exp,exp1,method="spearman")
    }
    max_cor <- max(cors)
    max_ident <- names(cors[cors==max_cor])
    result <- add_row(result,cluster=cluster,ident=max_ident,cor=max_cor)
  }

  trans <- result$ident
  names(trans) <- result$cluster
  all_data@meta.data$cor_ident <- trans[all_data@meta.data$res.0.8]

  all_data <- SetAllIdent(all_data,id="cor_ident")
  pdf("./cluster_level_22ident.pdf")
  TSNEPlot(all_data,colors.use = color2)
  dev.off()

  saveRDS(all_data,"./cluster_level.rds")

}



