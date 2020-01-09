### transform
tranform_meta <- function(x){
  
  x[x$samples=="S31T",]$samples <- "S31TA"
  x[x$samples %in% c("S31T_1","S31T_2"),]$samples <- "S31TB"
  x$patients <- str_replace(x$samples,pattern = "_[1,2]$",replacement = "")
  x$patients <- str_replace(x$patients,pattern = "^S",replacement = "")
  x$patients <- str_replace(x$patients,pattern = "T",replacement = "")
  x$patients <- str_replace(x$patients,pattern = "_$",replacement = "")
  x <- x[!x$patients %in% c("37","48"),]
  return (x)
}
meta <- rds@meta.data
meta <- tranform_meta(meta)
meta1 <- rds@meta.data
rds@meta.data <- meta

# 细胞大类
assign_big_type <- function(x){
  if (grepl("Myeloid",x)) return("Myeloid")
  if (grepl("T_cell",x)) return("T_cell")
  if (grepl("B_cell",x)) return("B_cell")
  if (grepl("LUAD",x)) return("LUAD")
  if (grepl("LUSC",x)) return("LUSC")
  if (grepl("Fibroblast",x)) return("Fibroblast")
  if (grepl("Langerhans",x)) return("Myeloid")
  if (grepl("Mast",x)) return("Myeloid")
  return(x)
}
x$big_type <- map_chr(x$new_ident,.f=assign_big_type)


### boxplot
dat <- as.tibble(freq_table_4)
colnames(dat) <- c("cell_type","patient","percent")

dic <- info5$mutation
names(dic) <- rownames(info5)
dat$mutation <- dic[dat$patient]

library(ggpubr)
pdf("./all/all_immune_MUTvsnoMUT_boxplot_3type.pdf",width=15)
dat1 <- filter(dat,mutation %in% c("EGFR 19-del","0","EGFR L858R"))
ggplot(group_by(dat1,cell_type),aes(x=cell_type,y=percent,fill=mutation)) +
  geom_boxplot(position = "dodge") + coord_flip() +
  guides(fill = guide_legend(reverse = TRUE)) + 
  stat_compare_means( method="anova",paired = F,
                      aes(label = paste0("p = ", ..p.format..))
  )
dev.off()


###feature boxplot
BoxPlot <- function(rds,features,nrow=1){
  library(ggplot2)
  library(ggpubr)
  library(Seurat)
  allDat <- tibble(barcode=character(),
                   exp=double(),
                   feature=character(),
                   ident = character())
  idents <- unique(rds@ident)
  for (ident in idents){
    cells <- WhichCells(rds,ident=ident)
    for (feature in features){
      exp <- data.frame(FetchData(rds,cells.use=cells,vars.all=feature))
      exp <- rownames_to_column(exp,var="barcode")
      colnames(exp)[[2]] <- "exp"
      exp$feature <- feature
      exp$ident <- ident
      allDat <- rbind(allDat,exp)
    }
  }
  allDat$feature <- factor(allDat$feature,levels=features)
  ggplot(allDat,aes(x=ident,y=exp)) + geom_boxplot() +
    stat_compare_means(method="wilcox.test",paired = F, aes(label = paste0("p = ", ..p.format..)))+
    facet_wrap(~feature,nrow=nrow)
}


### 基因表达
library(tidyverse)
mtx <- read_tsv("/SGRNJ01/RandD3/RD2019016_VDJ/20191120_sc/PBMC_S_P1/05.count/PBMC_S_P1_matrix.xls")
mtx
a <- rowSums(mtx)
a <- mtx %>% mutate(rowsum=rowSums(.[2:length(colnames(mtx))]))
b <- select(a,c("geneID","rowsum"))
c <- arrange(b,desc(rowsum))

gn <- read_tsv("/SGRNJ/Public/Database/genome/homo_sapiens/ensembl_92/Homo_sapiens.GRCh38.genenamefile.txt")
d <- left_join(c,gn,by=c("geneID"="Gene stable ID"))

### TSNEPlot
TSNEPlot(rds,group.by="samples")


### marker heatmap

##heatmap
meta.type = "res.0.6"
#source("/SGRNJ/Public/Script/shouhou/SC_scripts_zyq/anno_heatmap.R")
rds <- readRDS("/SGRNJ01/Aftersales/P2018009_Lung/20190624_jiace/paper/rerun/fibro/rds/all_TSNE_samples.rds")
#marker_test <- FindAllMarkers(rds,genes.use = rds@var.genes,max.cells.per.ident = 300)
marker_test <- read_csv("/SGRNJ01/Aftersales/P2018009_Lung/20190624_jiace/paper/rerun/fibro/csv/all_markers.csv")
top_gene <- marker_test %>% group_by(cluster) %>% top_n(10, avg_logFC)
top_gene <- top_gene[!duplicated(genes),]
genes <- top_gene$gene


rds <- ScaleData(rds,genes.use = genes)
data <- rds@scale.data[genes,]

cell_type <- factor(sort(unique(rds@ident)))
colors <- color1[c(1:length(cell_type))]
names(colors) <- cell_type

anno_col <- rds@meta.data[,meta.type,drop=F]
colnames(anno_col) <- "cell_type"
anno_col$cell_type <- factor(anno_col$cell_type,levels=unique(anno_col$cell_type))

anno_row <- as.data.frame(top_gene[,c("cluster")])
rownames(anno_row) <- top_gene$gene
colnames(anno_row) <- "cell_marker"
annotation_colors = list(cell_type=colors,cell_marker=colors)

barcode <- c()
genes <- c()
for (type in cell_type){
  print(type)
  cell <- WhichCells(rds,ident=type)
  barcode <- c(barcode,cell)
  gene <- filter(top_gene,cluster==type)$gene
  genes <- c(genes,gene)
}
data1 <- data[genes,barcode]

pheatmap(data1,cluster_row = FALSE, cluster_col = F,display_numbers =F,
         annotation_col=anno_col,anotation_names_row=FALSE,annotation_row=anno_row,
         show_colnames=F,annotation_colors=annotation_colors,annotation_names_row = F)


###gsva
Rscript /SGRNJ01/RD_dir/pipeline_test/yaodan/gsva1/gsva.R\
 --geneset /SGRNJ/Database/download/GSEA/h.all.v6.2.symbols.gmt --input ../avg.tsv --outdir ./ --pvalue 1 --lfc 0 --ftsize 5 --type 0


###cellphonedb
cell200 <- SubsetData(rds,max.cells.per.ident = 200,subset.raw = T)
data <- cell200@data
rownames(data)
library(tidyverse)
gn <- read_tsv("/SGRNJ01/Aftersales/P2019008_SCOPE.endometrium/P1.heatmap/cellphone/Homo_sapiens_Ensemble_92_gene.txt")
dic <- gn$gene_id
names(dic) <- gn$gene_name
gene.symbol <- rownames(data)
gene.en <- dic[gene.symbol]
table(is.na(gene.en))

gene.in <- gene.symbol[gene.symbol %in% names(dic)]
data1 <- data[gene.in,]
rownames(data1) <- dic[rownames(data1)]
data2 <- as.matrix(data1)
write.table(data2, "/Public/Script/shouhou/cellphonedb/test_lung/cell200/cell200_count.tsv", sep='\t', quote=F)

meta_data <- cbind(rownames(cell200@meta.data), cell200@meta.data[,'big_type', drop=F])   #####  cluster is the user’s specific cluster column
write.table(meta_data, "/Public/Script/shouhou/cellphonedb/test_lung/cell200/cell200_meta.tsv", sep='\t', quote=F, row.names=F)


## ITH

library(Seurat)
library(tidyverse)
library(reshape2)


#argv
patient_col = "patients"
ident = "LUAD"
rds_name = "/SGRNJ01/Aftersales/P2018016_Lung/paper1/rm_doublets/luad/LUAD_all.rds"

#run
rds <- readRDS(rds_name)
rds <- SubsetData(rds,ident.use = ident)
patients = unique(luad@meta.data[,patient_col])
dat = tibble(patient=character(),iqr=numeric())

for (patient in patients){
  pcell = WhichCells(rds,subset.name = patient_col,
                     accept.value=patient)
  data = rds@data[,pcell]
  corr = cor(as.matrix(data))
  melt_cor <- melt(corr)
  IQR = IQR(melt_cor$value)
  dat = add_row(dat,patient=patient,iqr=IQR)
}

write_tsv(dat,"./dat.tsv")