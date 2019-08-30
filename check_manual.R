#!/usr/bin/env Rscript
library(argparser)
library(tidyverse)

argv <- arg_parser('')
argv <- add_argument(argv,"--tsv",help="tsv name")
argv <- parse_args(argv)

dat <- read_tsv(argv$tsv)

if (TRUE %in% duplicated(dat$X1))
{
	stop ("duplicated cluster names")
}

dat$c <- paste0("C",dat$cluster)
dat$new_type <- paste(dat$c,dat$cell_type,sep="_")
dat1 <- dplyr::select(dat,cluster,new_type)
write_tsv(dat1,"./check_cluster_type.tsv")