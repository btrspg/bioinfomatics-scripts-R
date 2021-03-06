#! /usr/bin/env Rscript
# @License : Copyright(c) Aegicare License
# @Time    : 2019-08-10 10:33
# @Author  : YUELONG.CHEN
# @Mail    : yuelong.chen.btr@gmail.com
# @File    : deseq2
# @Software: PyCharm

library('DESeq2')
library('ggplot2')
library("RColorBrewer")
library('pheatmap')

args <- commandArgs(T)

featurecounts = args[1]
sampleinfo=args[2]
tag=args[3]
outdir=args[4]


data = read.table(featurecounts,sep='\t',header=TRUE)
count_data=as.matrix(data[,2:dim(data)[2]])
rownames(count_data)=data$Geneid

#sample.info
# AG18120501	Control	1812	c
# AG18120502	Control	1812	c
# AG18120503	Control	1812	c
# AG19060501	Control	1906	c
# AG19052202	Control	1905	c
# AG18120505	D	1812	d1
# AG18120506	D	1812	d1
# AG19060502	D	1906	d1
# AG19060503	D	1906	d1
# AG18120504	D	1812	d2
# AG19052205	D	1905	d2

sample.info=read.table(sampleinfo,sep='\t',header=FALSE,row.names=NULL,col.names=c('sample','condition','batch','type'))
sample.info=sample.info[match(colnames(count_data),sample.info$sample),]
dds <- DESeqDataSetFromMatrix(count_data, colData=sample.info, design= ~ condition + batch)
dds$condition <- relevel(dds$condition, ref = "Control")

dds <- dds[ rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds,test='LRT',reduced = ~ batch)
res <- results(dds)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
resdata <- data.frame(resdata[,-1], row.names=resdata[,1])
write.csv(resdata, paste(outdir,paste(tag,'des.all.output.csv',sep='-'),sep='/'), row.names=TRUE,quote=FALSE)
write.csv(res,paste(outdir,paste(tag,'des.output.csv',sep='-'),sep='/') , row.names=TRUE,quote=FALSE)

rld <- rlog(dds)
vsd <- vst(dds)
ntd <- normTransform(dds)



sv = rowVars(assay(rld))/rowMeans(assay(rld))
rv = rowVars(assay(rld))

select = order(rv, decreasing=TRUE)[seq_len(min(100, length(rv)))]

print(dim(rld[select,]))
print(dim(rld))

pcaData <- plotPCA(rld[select,], intgroup=c("condition" ,"type"), returnData=T)
percentVar <- round(100*attr(pcaData, "percentVar"))
pca <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
	geom_point(size=3) +
	ggtitle("DESeq2 PCA") +geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black") +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
	ylab(paste0("PC2: ", percentVar[2], "% variance"))
ggsave(paste(outdir,paste(tag,'pca.png',sep='-'),sep='/'),pca)

pcaData2 <- plotPCA(rld, intgroup=c("condition", "type"), returnData=T)
percentVar2 <- round(100*attr(pcaData2, "percentVar"))
pca2 <- ggplot(pcaData2, aes(PC1, PC2, color=condition, shape=type)) +
	geom_point(size=3) +
	ggtitle("DESeq2 PCA") +geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black") +
	xlab(paste0("PC1: ", percentVar2[1], "% variance")) +
	ylab(paste0("PC2: ", percentVar2[2], "% variance"))
ggsave(paste(outdir,paste(tag,'pca-all.png',sep='-'),sep='/'),pca2)


sv = rowVars(assay(vsd))/rowMeans(assay(vsd))
rv = rowVars(assay(vsd))

select = order(rv, decreasing=TRUE)[seq_len(min(100, length(rv)))]

print(dim(vsd[select,]))
print(dim(vsd))

pcaData <- plotPCA(vsd[select,], intgroup=c("condition" ,"type"), returnData=T)
percentVar <- round(100*attr(pcaData, "percentVar"))
pca <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
	geom_point(size=3) +
	ggtitle("DESeq2 PCA") +geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black") +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
	ylab(paste0("PC2: ", percentVar[2], "% variance"))
ggsave(paste(outdir,paste(tag,'vsd.pca.png',sep='-'),sep='/'),pca)

pcaData2 <- plotPCA(vsd, intgroup=c("condition", "type"), returnData=T)
percentVar2 <- round(100*attr(pcaData2, "percentVar"))
pca2 <- ggplot(pcaData2, aes(PC1, PC2, color=condition, shape=type)) +
	geom_point(size=3) +
	ggtitle("DESeq2 PCA") +geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black") +
	xlab(paste0("PC1: ", percentVar2[1], "% variance")) +
	ylab(paste0("PC2: ", percentVar2[2], "% variance"))
ggsave(paste(outdir,paste(tag,'vsd.pca-all.png',sep='-'),sep='/'),pca2)

output = paste(outdir,paste(tag,'rlog.pheatmap.pdf',sep='-'),sep='/')
output2 = paste(outdir,paste(tag,'rlog.pheatmap-test.pdf',sep='-'),sep='/')
output3 = paste(outdir,paste(tag,'vst.pheatmap.pdf',sep='-'),sep='/')
output4 = paste(outdir,paste(tag,'vst.pheatmap-test.pdf',sep='-'),sep='/')
sample_to_sample_output = paste(outdir,paste(tag,'sample2sample.pdf',sep='-'),sep='/')


newres.padj1 =na.omit(res)[na.omit(res)$padj<0.01,]
newres.padj5 =na.omit(res)[na.omit(res)$padj<0.05,]
newres.p1 =na.omit(res)[na.omit(res)$pvalue<0.01,]
newres.p5 =na.omit(res)[na.omit(res)$pvalue<0.05,]

threshold = ''
if(dim(newres.padj1)[1]>300){
    newres=newres.padj1
    threshold = 'padj(0.01)'
}else if(dim(newres.padj5)[1]>300){
    newres=newres.padj5
    threshold = 'padj(0.05)'
}else if(dim(newres.p1)[1]>300){
    newres=newres.p1
    threshold = 'pvalue(0.01)'
}else if(dim(newres.p5)[1]>300){
    newres=newres.p5
    threshold = 'pvalue(0.05)'
}else{
    newres=na.omit(res)
    threshold = 'ALL'
}

write.csv(resdata[rownames(newres),],paste(outdir,paste(tag,'des.output.sig.csv',sep='-'),sep='/') , row.names=TRUE,quote=FALSE)

colData=as.data.frame(colData(dds)[,colnames(colData(dds))])[,c('condition','type')]
pheatmap(assay(rld)[rownames(newres),],
        show_rownames=TRUE,
        cellwidth = 22,
        fontsize_row = 0.5,
        cluster_cols=TRUE,
        cluster_rows=TRUE,
        annotation_col=colData,
        scale = 'row',
        filename=output,
        main=paste('Heatmap of DEGs in ',threshold,sep=''))

pheatmap(assay(rld),
        show_rownames=TRUE,
        cellwidth = 22,
        fontsize_row = 0.5,
        cluster_cols=TRUE,
        cluster_rows=TRUE,
        annotation_col=colData,
        scale = 'row',
        filename=output2,
        main=paste('Test Heatmap of DEGs in ',threshold,sep=''))


pheatmap(assay(vsd)[rownames(newres),],
        show_rownames=TRUE,
        cellwidth = 22,
        fontsize_row = 0.5,
        cluster_cols=TRUE,
        cluster_rows=TRUE,
        annotation_col=colData,
        scale = 'row',
        filename=output3,
        main=paste('Heatmap of DEGs in ',threshold,sep=''))

pheatmap(assay(vsd),
        show_rownames=TRUE,
        cellwidth = 22,
        fontsize_row = 0.5,
        cluster_cols=TRUE,
        cluster_rows=TRUE,
        annotation_col=colData,
        scale = 'row',
        filename=output4,
        main=paste('Test Heatmap of DEGs in ',threshold,sep=''))

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, vsd$sample,sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         filename=sample_to_sample_output,
         main='Sample to Sample Heatmap')