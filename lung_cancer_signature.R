luad <- readRDS("/SGRNJ01/Aftersales/P2018009_Lung/20190624_jiace/paper/rerun/luad/LUAD.rds")

genes <- readRDS("/SGRNJ/Database/sc_marker/alveolar/alveolar_gene.rds")
a_gene <- genes$normal
l_gene <- genes$cancer
TSNEPlot(luad,do.label = T)

gene <- c(a_gene,l_gene)
y <- luad@data[gene,]
a_mean <- colMeans(as.matrix(y[a_gene,]), na.rm = FALSE, dims = 1)
l_mean <- colMeans(as.matrix(y[l_gene,]), na.rm = FALSE, dims = 1)
a_mean[a_mean==0] <- min(a_mean[a_mean!=0])
l_mean[l_mean==0] <- min(l_mean[l_mean!=0])
res <- (log(a_mean,base = 2) - log(l_mean,base=2))
luad@meta.data$signature <- res
FeaturePlot(luad,"signature",no.legend = F,cols.use = c("red","green"))
TSNEPlot(luad,group.by="samples",no.legend = F,colors.use = color2,do.label=T)
