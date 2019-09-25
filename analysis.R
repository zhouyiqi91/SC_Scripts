## transform
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

## 细胞大类
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


## boxplot
dat <- as.tibble(freq_table_4)
colnames(dat) <- c("cell_type","patient","percent")

dic <- info5$mutation
names(dic) <- rownames(info5)
dat$mutation <- dic[dat$patient]

pdf("./all/all_immune_MUTvsnoMUT_boxplot_3type.pdf",width=15)
dat1 <- filter(dat,mutation %in% c("EGFR 19-del","0","EGFR L858R"))
ggplot(group_by(dat1,cell_type),aes(x=cell_type,y=percent,fill=mutation)) +
  geom_boxplot(position = "dodge") + coord_flip() +
  guides(fill = guide_legend(reverse = TRUE)) + 
  stat_compare_means( method="anova",paired = F,
                      aes(label = paste0("p = ", ..p.format..))
  )
dev.off()