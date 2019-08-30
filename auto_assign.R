#!/usr/bin/env Rscript
library(argparser)
library(Seurat)
library(tidyverse)

argv <- arg_parser('')
argv <- add_argument(argv,"--outdir", help="the output dir.",default=getwd())
argv <- add_argument(argv,"--rds",help="Seurat rds")
argv <- add_argument(argv,"--type_marker_tsv",help="cell type marker tsv")
#argv <- add_argument(argv,"--resolution", help="tSNE resolution",default=0.8)
argv <- parse_args(argv)

#read args

outdir <- argv$outdir
rds <- argv$rds
type_marker_tsv <- argv$type_marker_tsv
#resolution <- argv$resolution
rds <- argv$rds
#origin.cluster <- paste("res.",resolution,sep="")

all_data <- readRDS(rds)
marker_file <- read.table(type_marker_tsv,header = TRUE,sep="\t",row.names = 1,stringsAsFactors=FALSE)

setwd(outdir)

#
cell_name <- rownames(marker_file)
n_cell_name <- length(cell_name)

#reset
#all_data <- SetAllIdent(object = all_data, id = origin.cluster)
clusters <- sort(unique(all_data@ident))

#create dir
auto_dir <- paste0(outdir,"/auto_assign/")
auto_png <- paste0(auto_dir,"png/")
dir.create(auto_dir)
dir.create(auto_png)

# type marker
setwd(auto_dir)
c = 0 
for (cluster in clusters){
  index = 0
  for (cell in cell_name){
    index = index + 1
    features = unlist(strsplit(marker_file[index,1],","))
    for (feature in features){
      tryCatch({
        dat <- FindMarkers(all_data,genes.use=feature,ident.1=cluster,min.pct = 0,logfc.threshold = -Inf)
        dat$cell_type <- cell
        dat$cluster <- cluster
        dat <- rownames_to_column(dat,var="gene")
        if (c==0){
          all_dat <- dat
          c = c + 1
        } else {
          all_dat <- rbind(all_dat,dat)
          }
        }
        ,error=function(e){print(paste0(feature," not found in cluster ",cluster)) })
    }
  }
}

all_dat <- mutate(all_dat,pct.diff=pct.1-pct.2)
write_tsv(all_dat,"type_marker_exp.tsv")

# plot
color2 <- c("OrangeRed","SlateBlue3","DarkOrange","GreenYellow","Purple",
"DarkSlateGray","Gold","DarkGreen","DeepPink2","Red4",
"#4682B4","#FFDAB9","#708090","#836FFF","#CDC673",
"#CD9B1D","#FF6EB4","#CDB5CD","#008B8B","#43CD80",
"#6A5ACD","#66CD00","#CDC673","#CDAD00","#CD9B9B","#FF8247","#8B7355","#8B3A62",
"#68228B","#CDB7B5","#CD853F","#6B8E23","#E6E6FA","#FFDAB9","Green")

exp <- all_dat
a <- group_by(exp,cluster,cell_type)
setwd(auto_png)
for (cluster in clusters){
  c = a[a$cluster==cluster,]

  png(paste0(cluster,"_pctdiff.png"),width=1200,height=1000)
  p1 <- ggplot(c,aes(x=interaction(gene,cell_type),y=pct.diff,fill=cell_type)) +geom_bar(stat="identity")+ coord_flip() + scale_fill_manual(values=color2)
  print (p1)
  dev.off()

  png(paste0(cluster,"_logfc.png"),width=1200,height=1000)
  p2 <- ggplot(c,aes(x=interaction(gene,cell_type),avg_logFC,fill=cell_type)) +geom_bar(stat="identity")+ coord_flip() + scale_fill_manual(values=color2)
  print (p2)
  dev.off()
}

# auto assign
setwd(auto_dir)
as <- summarize(a,avg_pct.diff=mean(pct.diff),avg_logfc=mean(avg_logFC),avg_p_val_adj=mean(p_val_adj))    
as1 <- group_by(ungroup(as),cluster)
as1 <- mutate(as1,pct_rank = rank(avg_pct.diff),
              logfc_rank= rank(avg_logfc),total_rank=pct_rank+logfc_rank)
as2 <- as1 %>% ungroup %>% group_by(cluster) %>% 
  filter(total_rank==max(total_rank)) %>% arrange(as.numeric(cluster))
as3 <- select(as2,cluster,cell_type,avg_pct.diff,avg_logfc,avg_p_val_adj)
write_tsv(as3,"auto_cluster_type.tsv")