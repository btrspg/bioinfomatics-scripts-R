#! /usr/bin/env Rscript
# @License : Copyright(c) Aegicare License
# @Time    : 2019-08-10 17:14
# @Author  : YUELONG.CHEN
# @Mail    : yuelong.chen.btr@gmail.com
# @File    : enrichment.R
# @Software: PyCharm



args<-commandArgs(T)

deg_file = args[1]
prefix = args[2]
sps=args[3]


library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)
## methods

species='hsa'



load_data = function(csv_file){
    data = read.csv(csv_file,header=TRUE)
    return(data)
}

filter_data = function(deg_data,threshold_padj=0.05,threshold_fc=1,both=TRUE){
    if(both){
        deseq2.sig <- subset(deg_data, padj < threshold_padj & abs(log2FoldChange) > threshold_fc)


    }else if(threshold_fc>0){
        deseq2.sig <- subset(deg_data, padj < threshold_padj & log2FoldChange > threshold_fc)

    }else{
        deseq2.sig <- subset(deg_data, padj < threshold_padj & log2FoldChange < threshold_fc)

    }
    genelist = deseq2.sig$log2FoldChange
    names(genelist) <- deseq2.sig$GENEID
    genelist <- sort(genelist, decreasing = TRUE)
    return(genelist)

}

annotation_go = function(genelist,anno_db,ont="CC"){
    ego <- enrichGO(gene=names(genelist),
                OrgDb         = anno_db,
                keyType       = 'ENSEMBL',
                ont           = ont,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable=TRUE)

    # ego <- setReadable(ego, OrgDb = anno_db)
    return(ego)
}

annotation_kegg = function(gene_df,anno_db,org='hsa'){
    gene = bitr(names(gene_df), fromType = "ENSEMBL", toType=c("ENTREZID","SYMBOL"), anno_db, drop = TRUE)
    ekegg<-enrichKEGG(gene$ENTREZID, organism = org, keyType = "kegg", pvalueCutoff = 0.05,
            pAdjustMethod = "BH",  minGSSize = 10, maxGSSize = 500,
            qvalueCutoff = 0.2, use_internal_data = FALSE)
    return(ekegg)
}

plot_anno = function(ego,title){




        dotplot(ego, showCategory=30,title=title)



}




    # enrichMap(ego, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai,title=paste(deg_type,ont,sep='-'))
    # cnetplot(ego, foldChange=geneList,title=paste(deg_type,ont,sep='-'))
    # plotGOgraph(ego,title=paste(deg_type,ont,sep='-'))


gsea_anno = function(gene_list,anno_db,org="hsa"){
    gsego <- gseGO(geneList=gene_list, keyType = "ENTREZID",ont="ALL", OrgDb=anno_db, verbose=F)
    gsekegg<-gseKEGG(geneList=gene_list, organism = org, keyType = "kegg", exponent = 1,
        nPerm = 1000, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05,
        pAdjustMethod = "BH", verbose = TRUE, use_internal_data = FALSE,
        seed = FALSE, by = "fgsea")
}

plot_kegg=function(ekegg,foldChange=1){

    viewKEGG(ekegg, pathwayID, foldChange, color.low = "green",
        color.high = "red", kegg.native = TRUE, out.suffix = "clusterProfiler")

}

writer_annotation=function(anno,output){
    write.csv(as.data.frame(anno),output,quote=FALSE)
}


## pipeline

anno_db = load_db(sps)
deg_df = load_data(deg_file)
kegg_sps=which_species(sps)

deg_up_df=filter_data(deg_df,threshold_padj=0.05,threshold_fc=1.5,both=FALSE)
deg_down_df=filter_data(deg_df,threshold_padj=0.05,threshold_fc=-1.5,both=FALSE)
deg_all_df=filter_data(deg_df,threshold_padj=0.05,threshold_fc=1.5,both=TRUE)


## 上调基因GO富集
print("GO UP enrichment")
up_go_cc = annotation_go(deg_up_df,anno_db,ont="CC")
up_go_bp = annotation_go(deg_up_df,anno_db,ont="BP")
up_go_mf = annotation_go(deg_up_df,anno_db,ont="MF")
up_go_all = annotation_go(deg_up_df,anno_db,ont="ALL")

## 下调基因GO富集
print("GO DOWN enrichment")
down_go_cc = annotation_go(deg_down_df,anno_db,ont="CC")
down_go_bp = annotation_go(deg_down_df,anno_db,ont="BP")
down_go_mf = annotation_go(deg_down_df,anno_db,ont="MF")
down_go_all = annotation_go(deg_down_df,anno_db,ont="ALL")

## 差异基因（所有）GO富集
print("GO ALL enrichment")
all_go_cc = annotation_go(deg_all_df,anno_db,ont="CC")
all_go_bp = annotation_go(deg_all_df,anno_db,ont="BP")
all_go_mf = annotation_go(deg_all_df,anno_db,ont="MF")
all_go_all = annotation_go(deg_all_df,anno_db,ont="ALL")


## 差异基因KEGG富集

print("KEGG ALL enrichment")
up_kegg=annotation_kegg(deg_up_df,anno_db,org=kegg_sps)
down_kegg=annotation_kegg(deg_down_df,anno_db,org=kegg_sps)
all_kegg=annotation_kegg(deg_all_df,anno_db,org=kegg_sps)


print("WRITER enrichment")
writer_annotation(up_go_cc,paste(prefix,"GO_CC_UP.csv",sep='-'))
writer_annotation(up_go_bp,paste(prefix,"GO_BP_UP.csv",sep='-'))
writer_annotation(up_go_mf,paste(prefix,"GO_MF_UP.csv",sep='-'))
writer_annotation(up_go_all,paste(prefix,"GO_ALL_UP.csv",sep='-'))

writer_annotation(down_go_cc,paste(prefix,"GO_CC_DOWN.csv",sep='-'))
writer_annotation(down_go_bp,paste(prefix,"GO_BP_DOWN.csv",sep='-'))
writer_annotation(down_go_mf,paste(prefix,"GO_MF_DOWN.csv",sep='-'))
writer_annotation(down_go_all,paste(prefix,"GO_ALL_DOWN.csv",sep='-'))

writer_annotation(all_go_cc,paste(prefix,"GO_CC_ALL.csv",sep='-'))
writer_annotation(all_go_bp,paste(prefix,"GO_BP_ALL.csv",sep='-'))
writer_annotation(all_go_mf,paste(prefix,"GO_MF_ALL.csv",sep='-'))
writer_annotation(all_go_all,paste(prefix,"GO_ALL_ALL.csv",sep='-'))


writer_annotation(setReadable(up_kegg, anno_db, keytype="ENTREZID"),paste(prefix,"KEGG_UP.csv",sep='-'))
writer_annotation(setReadable(down_kegg, anno_db, keytype="ENTREZID"),paste(prefix,"KEGG_DOWN.csv",sep='-'))
writer_annotation(setReadable(all_kegg, anno_db, keytype="ENTREZID"),paste(prefix,"KEGG_all.csv",sep='-'))

print("PLOT enrichment")




# go_plot=c(up_go_cc,up_go_bp,up_go_mf,up_go_all,
#             down_go_cc,down_go_bp,down_go_mf,down_go_all,
#             all_go_cc,all_go_bp,all_go_mf,all_go_all)
# plot_name=c('up_go_cc','up_go_bp','up_go_mf','up_go_all',
#             'down_go_cc','down_go_bp','down_go_mf','down_go_all',
#             'all_go_cc','all_go_bp','all_go_mf','all_go_all')
#
# kegg_plot=c(up_kegg,down_kegg,all_kegg)
# plot_kegg_name=c('up_kegg','down_kegg','all_kegg')


pdf(paste(prefix,"annotation.pdf",sep="-"),width=15,height=9)
if(dim(as.data.frame(up_go_cc))[1]>1){
    dotplot(up_go_cc, showCategory=30,title='up_go_cc')
}
if(dim(as.data.frame(up_go_bp))[1]>1){
    dotplot(up_go_bp, showCategory=30,title='up_go_bp')
}
if(dim(as.data.frame(up_go_mf))[1]>1){
    dotplot(up_go_mf, showCategory=30,title='up_go_mf')
}
if(dim(as.data.frame(up_go_all))[1]>1){
    dotplot(up_go_all, showCategory=30,title='up_go_all')
}
if(dim(as.data.frame(down_go_bp))[1]>1){
    dotplot(down_go_bp, showCategory=30,title='down_go_bp')
}
if(dim(as.data.frame(down_go_mf))[1]>1){
    dotplot(down_go_mf, showCategory=30,title='down_go_mf')
}
if(dim(as.data.frame(down_go_cc))[1]>1){
    dotplot(down_go_cc, showCategory=30,title='down_go_cc')
}
if(dim(as.data.frame(down_go_all))[1]>1){
    dotplot(down_go_all, showCategory=30,title='down_go_all')
}
if(dim(as.data.frame(all_go_cc))[1]>1){
    dotplot(all_go_cc, showCategory=30,title='all_go_cc')
}
if(dim(as.data.frame(all_go_bp))[1]>1){
    dotplot(all_go_bp, showCategory=30,title='all_go_bp')
}
if(dim(as.data.frame(all_go_mf))[1]>1){
    dotplot(all_go_mf, showCategory=30,title='all_go_mf')
}
if(dim(as.data.frame(all_go_all))[1]>1){
    dotplot(all_go_all, showCategory=30,title='all_go_all')
}
if(dim(as.data.frame(up_kegg))[1]>1){
    dotplot(up_kegg, showCategory=30,title='up_kegg')
}
if(dim(as.data.frame(down_kegg))[1]>1){
    dotplot(down_kegg, showCategory=30,title='down_kegg')
}
if(dim(as.data.frame(all_kegg))[1]>1){
    dotplot(all_kegg, showCategory=30,title='all_kegg')
}





# for(pl in 1:length(plot_name)){
#     if(dim(as.data.frame(go_plot[[pl]]))[1]> 1){
#         print(plot_name[pl])
#         pdf(paste(prefix,plot_name[pl],"annotation.pdf",sep="-"),width=20,height=18)
#         dotplot(go_plot[[pl]], showCategory=30,title=plot_name[pl])
#         dev.off()
        # emapplot(go_plot[[pl]])
        # cnetplot(go_plot[[pl]],categorySize="pvalue",title=plot_name[pl])
        # goplot(go_plot[[pl]])
        # plotGOgraph(go_plot[[pl]])
    # }
# }

# for(pl in 1:length(plot_kegg_name)){
#     if(dim(as.data.frame(kegg_plot[[pl]]))[1]>1){
#         print(plot_kegg_name[pl])
#         print(kegg_plot[[pl]])
#         pdf(paste(prefix,plot_kegg_name[pl],"annotation.pdf",sep="-"),width=20,height=18)
#         dotplot(kegg_plot[[pl]], showCategory=30,title=plot_kegg_name[pl])
#         dev.off()
#     }
    # emapplot(kegg_plot[pl])
    # cnetplot(kegg_plot[pl],categorySize="pvalue",title=plot_name[pl])
    # goplot(kegg_plot[pl])
    # plotGOgraph(kegg_plot[pl])
# }


dev.off()


