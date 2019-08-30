#!/usr/bin/env Rscript
library(argparser)
library(tidyverse)
library(Seurat)

argv <- arg_parser('')
argv <- add_argument(argv,"--rds",help="Seurat rds")
argv <- add_argument(argv,"--detection_rate",help="detection_rate",default=0.5)
argv <- add_argument(argv,"--save_rds",help="save no doublet rds",default="Y")
argv <- add_argument(argv,"--maxcell_perident",help="detection_rate",default=5000)
argv <- parse_args(argv)


#read args
print ("read rds.")
rds <- readRDS(argv$rds)
print ("read done.")
detection_rate <- argv$detection_rate
save_rds <- argv$save_rds

meta <- rds@meta.data
idents <- unique(meta$new_ident)

top1500 <- head(rownames(rds@hvg.info),1500)
print ("find markers")
markers <- FindAllMarkers(rds,genes.use = top1500,max.cells.per.ident=argv$maxcell_perident,min.diff.pct=0.5,only.pos=T)
print ("find markers done.")
markers$pct.diff <- markers$pct.1 - markers$pct.2

doublets <- list()
for (ident1 in idents)
{
  rest_ident <- idents[!idents %in% ident1]
  for (ident2 in rest_ident)
  {
    list_name <- paste(ident1,ident2,sep="_")
    print (list_name)
    ident2_marker <- markers[markers$cluster==ident2,]
    ident2_top_genes <- ident2_marker$gene
    print (ident2_top_genes)
    dropout_percent <- mean(ident2_marker$pct.1)
    find_number <- round(length(ident2_top_genes)*dropout_percent*detection_rate) 
    print (find_number)
    if (is.nan(find_number)) {next}
    if (find_number==0) {find_number <- 1}
    ident1_barcode <- rownames(meta[meta$new_ident==ident1,])
    ident_data <- rds@raw.data[ident2_top_genes,ident1_barcode]
    #percent10 <- round(length(ident1_barcode)/10)
    
    gene_list <- list()
    barcode_count <- c()
    
    for (gene in ident2_top_genes)
    {
      gene_data <- ident_data[gene,]
      mean_umi <- mean(gene_data)
      gene_data <- gene_data[gene_data>1 & gene_data>2*mean_umi]
      #gene_data <- gene_data[order(gene_data,decreasing = T)]
      #gene_barcode <- head(names(gene_data),percent10)
      gene_barcode <- (names(gene_data))
      for (barcode in gene_barcode)
      {
        if (!barcode %in% names(barcode_count))
        {
          barcode_count[barcode] <- 1
        } else 
        {
          barcode_count[barcode] <- barcode_count[barcode] + 1
        }
      }
    }
    doublet_barcode <- names(barcode_count[barcode_count>=find_number])
    doublets[[list_name]] <- doublet_barcode
  }
}

doublets_vector <- c()
for (item in doublets)
{
	doublets_vector <- c(doublets_vector,item)
}

named_vector <- c()
for (item in names(doublets))
{
  for (barcode in doublets[[item]])
  {
    
    if ( barcode %in% names(named_vector) ){
      value <- paste(named_vector[barcode],item,sep=" , ")
      named_vector[barcode] <- value
    } else{
      named_vector[barcode] <- item
    }
  }
}
named_dat <- data.frame(barcode=names(named_vector),identity=named_vector)
write_tsv(named_dat,"doublets_identity.txt")

saveRDS(doublets,"doublets.rds")

meta$is_doublets <- "normal"
meta[doublets_vector,]$is_doublets <- "doublets"
sample_dat <- group_by(meta,orig.ident) %>% 
	summarize(doublet_prop=length(is_doublets[which(is_doublets=="doublets")])/n()*100)
total_doublets <- sum(meta$is_doublets=="doublets")/length(meta$is_doublets)*100
sample_dat <- add_row(sample_dat,orig.ident="total",doublet_prop=total_doublets)
write_tsv(sample_dat,"sample_doublet.tsv")
pdf("doublet_sample.pdf",width=15)
ggplot(sample_dat,aes(x=orig.ident,y=doublet_prop)) + geom_bar(stat="identity") +
 theme(axis.text.x = element_text(angle = 90))
dev.off()


pdf("origin_TSNE.pdf",width=15)
TSNEPlot(rds,label.size = 8,do.label=T,pt.size=0.2)
dev.off()

all_barcode <- rownames(meta)
cells.use <- all_barcode[!all_barcode %in% doublets_vector]
nodoublet <- SubsetData(rds,cells.use = cells.use)
pdf("remove_TSNE.pdf",width=15)
TSNEPlot(nodoublet,label.size = 8,do.label=T,pt.size=0.2)
dev.off()
if (save_rds=="Y") saveRDS(nodoublet,"nodoublet.rds")

rds@meta.data <- meta
rds <- SetAllIdent(rds,"is_doublets")
pdf("doublets_TSNE.pdf",width=15)
TSNEPlot(rds,label.size = 8,pt.size=0.2)
dev.off()

