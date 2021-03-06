#! /usr/bin/env Rscript
# @License : Copyright(c) Aegicare License
# @Time    : 2019-09-08 12:48
# @Author  : YUELONG.CHEN
# @Mail    : yuelong.chen.btr@gmail.com
# @File    : ${NAME}
# @Software: PyCharm




read_htseq <- function(featurecounts_file,sample_id) {
    df = read.table(featurecounts_file,sep='\t',header=FALSE,comment.char = "_")
    colnames(df) = c('Geneid',sample_id)
    return(df)
}

args <- commandArgs(T)

my_list_file = args[1]
my_output_file = args[2]


mylist=read.table(my_list_file,sep='\t',header=FALSE)

merge_data = read_htseq(toString(mylist[1,2]),toString(mylist[1,1]))

for(i in 2:dim(mylist)[1]){
    temp = read_htseq(toString(mylist[i,2]),toString(mylist[i,1]))
    merge_data = merge(merge_data,temp,by.x='Geneid',by.y='Geneid',all.x=TRUE,all.y=TRUE)
}

write.table(merge_data,my_output_file,sep='\t',row.names=FALSE,col.names=TRUE,quote=FALSE)

