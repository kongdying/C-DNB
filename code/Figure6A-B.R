#本段代码按运行环境分为四部分
#R部分
library(Seurat)
load("step1/NSC.rds")
data_use <- GetAssayData(NSC, assay = "RNA", layer = "counts")
rownames(data_use) <- rownames(NSC)
colnames(data_use) <- colnames(NSC)
write.csv(t(as.matrix(data_use)),file="step6/Pe_4_scenic_input.csv")

#python部分
#参考：python=3.8 pyscenic=0.12.1 numpy=1.21 pandas=1.5.3 loompy==3.0.8/3.0.9
import os, sys
import loompy as lp;
import numpy as np;
import scanpy as sc;
x=sc.read_csv("Pe_4_scenic_input.csv");#R中导出的表达矩阵
row_attrs = {"Gene": np.array(x.var_names),};
col_attrs = {"CellID": np.array(x.obs_names)};
lp.create("Pe_4.loom",x.X.transpose(),row_attrs,col_attrs)

#window命令行部分
#参考：python=3.8 pyscenic=0.12.1 numpy=1.21 pandas=1.5.3 loompy==3.0.8/3.0.9
conda activate pyscenic
pyscenic grn Pe_4.loom ./human/hs_hgnc_tfs.txt -o Pe_4_scenic_input.adj.tsv --num_workers 5
pyscenic ctx Pe_4_scenic_input.adj.tsv ./human/hg19-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather -o Pe_4_scenic_input.motif.tsv --annotations_fname ./human/motifs-v9-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname Pe_4.loom --mode dask_multiprocessing --num_workers 5
pyscenic aucell Pe_4.loom Pe_4_scenic_input.motif.tsv  --output Pe_4_SCENIC.loom --num_workers 10

#R部分
library(SCENIC)#1.3.1
library(Seurat)#5.4.0
library(SCopeLoomR)#0.13.0
library(AUCell)#1.32.0
library(dplyr)#1.2.0
rm(list=ls())
#Figure6A-B
run_SCENIC_TF_analysis<-function(regulons,group,geneall,n){
  net=c()
  dnb_tf=c()
  for(tf in names(regulons)){
    for(gene in regulons[[tf]]){
      if(is.element(gene,geneall[,1])){
        tf=gsub("[(+)]","",tf)
        dnb_tf=c(dnb_tf,tf)
        net=rbind(net,c(tf,gene))
      }
    }
  }
  dnb_tf<-unique(dnb_tf)
  regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
  sorted_cell_names <- names(group)[order(group)]
  regulonAUC <- regulonAUC[, sorted_cell_names, drop = FALSE]
  annotation_col = data.frame(DataType = group[order(group)])
  rownames(annotation_col)<-colnames(regulonAUC)
  tfname=names(regulonAUC)
  tfname_use=gsub("[(+)]","",tfname)
  index=is.element(tfname_use,dnb_tf)
  aucvalue=(regulonAUC@assays@data@listData$AUC)
  plot_value=aucvalue[index,]
  allpvaluegreater=c()#大于其他时间点的p值
  allpvalueless=c()#小于其他时间点的p值
  retain=c()#显著差异的tf
  up_tf=c()#显著高活性的tf，主要关注这个
  down_tf=c()#显著低活性的tf
  p0=c()#与T0时间点的p值
  p1=c()#以此类推
  p2=c()
  p3=c()
  p4=c()
  in_group<-annotation_col$DataType==paste0("T",n)
  for(i in 1:dim(plot_value)[1]){
    p_greater=wilcox.test(plot_value[i,in_group],plot_value[i,!in_group],alternative ="greater")
    p_less=wilcox.test(plot_value[i,in_group],plot_value[i,!in_group],alternative ="less")
    allpvaluegreater=c(allpvaluegreater,p_greater$p.value)
    allpvalueless=c(allpvalueless,p_less$p.value)
    p0=c(p0,wilcox.test(plot_value[i,in_group],plot_value[i,annotation_col$DataType=="T0"],alternative ="greater")$p.value)
    p1=c(p1,wilcox.test(plot_value[i,in_group],plot_value[i,annotation_col$DataType=="T1"],alternative ="greater")$p.value)
    p2=c(p2,wilcox.test(plot_value[i,in_group],plot_value[i,annotation_col$DataType=="T2"],alternative ="greater")$p.value)
    p3=c(p3,wilcox.test(plot_value[i,in_group],plot_value[i,annotation_col$DataType=="T3"],alternative ="greater")$p.value)
    p4=c(p4,wilcox.test(plot_value[i,in_group],plot_value[i,annotation_col$DataType=="T4"],alternative ="greater")$p.value)
    if(p_greater$p.value<0.01){
      retain=c(retain,i)
      up_tf=c(up_tf,rownames(plot_value)[i])
    }
    if(p_less$p.value<0.01){
      retain=c(retain,i)
      down_tf=c(down_tf,rownames(plot_value)[i])
    }
  }
  up_tf_use=gsub("[(+)]","",up_tf)
  write.table(up_tf_use,paste0("t",n,"_up_tf.txt"),col.names=F,row.names=F,quote=F)
  
  #做up_tf_gene表格
  net_up_tf=net[is.element(net[,1],up_tf_use),]
  colnames(net_up_tf)=c("tf","genes")
  write.table(net_up_tf,paste0("t",n,"_up_tf_gene_network.txt"),col.names=T,row.names=F,quote=F,sep=",")
  
  library("pheatmap")
  annotation_colors <- list(DataType= c("T0" = "#88DAFF", "T1" = "#8BFF88","T2" = "#FFF688", "T3" = "#FDC086","T4" = "#FF8888"))
  bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
  
  p <- pheatmap(plot_value[up_tf,],scale="row",
                show_colnames=F,cluster_cols=F,cluster_rows=T,
                color = c(colorRampPalette(colors = c("#366BAE","#619AC7","#FFFFFF"))(length(bk)/2),
                          colorRampPalette(colors = c("#FFFFFF","#F18C64","#B12424"))(length(bk)/2)),
                legend_breaks=seq(-2,2,1),breaks=bk,
                annotation_col=annotation_col,annotation_colors=annotation_colors)
  clustered_up_tf <- up_tf[p$tree_row$order]
  
  pheatmap(plot_value[clustered_up_tf,],scale="row",
           show_colnames=F,cluster_cols=F,cluster_rows=F,
           color = c(colorRampPalette(colors = c("#366BAE","#92BFDA","#FFFFFF"))(length(bk)/2),
                     colorRampPalette(colors = c("#FFFFFF","#D9A691","#B12424"))(length(bk)/2)),
           legend_breaks=seq(-2,2,1),breaks=bk,
           annotation_col=annotation_col,annotation_colors=annotation_colors,
           filename = paste0("t",n,"_scenic_heatmap_up_tf_cluster.pdf"))
  
  df<-data.frame(p0,p1,p2,p3,p4,allpvaluegreater,row.names = rownames(plot_value))
  df1<-apply(df,2,function(x){
    vn<-ifelse(x<0.0001,"****",x)
    vn<-ifelse(x>0.0001&x<0.001,"***",vn)
    vn<-ifelse(x>0.001&x<0.01,"**",vn)
    vn<-ifelse(x>0.01&x<0.05,"*",vn)
    vn<-ifelse(x>0.05,"ns",vn)
  })
  write.csv(df[clustered_up_tf,],paste0("t",n,"_wilcoxtest_pvalue.csv"),row.names=T,quote=F)
  write.csv(df1[clustered_up_tf,],paste0("t",n,"_wilcoxtest_p_state.csv"),row.names=T,quote=F)
}

setwd("step6")
loom <- open_loom("Pe_4_SCENIC.loom")
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
load('../step1/NSC.rds')
group<-NSC$stimgroup
t1_dnb <-  read.table("../step2/top50DNB_T1.txt",sep = '\t',header = T)
t3_dnb <-  read.table("../step2/top50DNB_T3.txt",sep = '\t',header = T)
outdir= "./result"
dir.create(outdir)
setwd(outdir)
run_SCENIC_TF_analysis(regulons,group,t1_dnb,1)
run_SCENIC_TF_analysis(regulons,group,t3_dnb,3)
setwd("..")
