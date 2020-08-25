#! /usr/bin/env Rscript
# Title     : heatmap.R
# Objective :
# Created by: logan
# Created on: 2020/8/25

args<-commandArgs(T)

library('pheatmap')

mt <- read.csv(args[1],sep='\t',row.names=1)
info <- read.csv(args[2],sep='\t',row.names=1)

pdf(args[3])

pheatmap(mt[order(info)],annotation_col=info)
dev.off()
