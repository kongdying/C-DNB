library(dplyr)
library(reshape2)
library(psych)
library(Seurat)
library(igraph)
library(harmony)
library(clustree)
require(magrittr)
require(Matrix)
library(cowplot)
library("reshape")
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(ggsci)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(patchwork)
library(Mfuzz)


outdir= "./step1"
dir.create(outdir)
setwd(outdir)
sampleinfo=read.table("../data/sample.txt",sep="\t",header=T,comment.char = "*")
count <- read.table("../data/count.txt", sep = "\t", header = T)
All=CreateSeuratObject(counts = count,min.cells = 5)
identical(rownames(sampleinfo),colnames(All))
All$stim=factor(sampleinfo[,"sample.label"], levels = c("M0_1","M1_1","M1_2","M1_3","M1_4","M2_1","M2_2","End_1","End_2","End_SVZ_1","End_SVZ_2")) #样本标签
All$stimgroup=factor(sampleinfo[,"sample.time"], levels = c("M0","M1","M2","End","End_SVZ")) #分组标签
All[["percent.mt"]] <- PercentageFeatureSet(All, pattern = "^MT-")

All$stimgroup <- factor(All$stimgroup, 
                        levels = c("M0", "M1", "M2", "End", "End_SVZ"), 
                        labels = c("T0", "T1", "T2", "T3", "T4"))

all.genes <- rownames(All)
genes_to_remove  <- grep(pattern = "^MT-|^RP[SL]", x = all.genes, value = TRUE)
gene_names <- setdiff(all.genes, genes_to_remove )
All <- subset(All, features = gene_names)
save(All, file ="All.rds")

All=All%>%Seurat::NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000)
All=All %>%  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) 
All=All %>% ScaleData(verbose = FALSE)
All=All %>% RunPCA(features = VariableFeatures(object = All), verbose = FALSE)
Idents(All) <- sampleinfo$finalCellType
cell_type <- unique(sampleinfo$finalCellType)
color <- brewer.pal(n = length(cell_type), name = 'Paired') 
color <- colorRampPalette(brewer.pal(12, "Set3"))(length(cell_type))
Idents(All) = factor(Idents(All),levels=unique(sampleinfo$finalCellType))
#NSC <- subset(All, idents = "NSC1(G2/M)")
#######################run harmony
All <- All %>%  RunHarmony("stim", plot_convergence = FALSE)
harmony_embeddings <- Embeddings(All, 'harmony')
All <- All %>% RunTSNE(reduction = "harmony", dims = 1:30) #harmony
All <- All %>% RunUMAP(reduction = "harmony", dims = 1:30) #harmony
All <- All %>% FindNeighbors(reduction = "harmony", dims = 1:30)#harmony
NSC <- subset(All, idents = "NSC1(G2/M)")
##########################celltype annotation
pdf("after_harmony_tSNE_All.pdf")#Figure 2a
DimPlot(All, reduction = "tsne", label = TRUE,pt.size = 0.5,cols=color,raster=FALSE)
dev.off()

pdf("after_harmony_tsne_NSC.pdf")#Figure 2b
DimPlot(NSC,reduction = "tsne", pt.size = 0.5, group.by = "stimgroup",
        cols =color,label = FALSE, raster = FALSE) 
dev.off()

save(NSC, file ="NSC.rds")

#############################top20Marker Figure 2d
Idents(NSC) <- NSC$stimgroup
All_markers <- FindAllMarkers(NSC, min.pct=0.25, logfc.threshold=0.25, only.pos=T)
top20markers <- All_markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC)
top20markers <- top20markers[!duplicated(top20markers$gene), ]
NPCData <- GetAssayData(NSC, assay="RNA", layer="data")
linshi <- NPCData[rownames(NPCData) %in% top20markers$gene, ]
linshi <- apply(as.matrix(linshi), 1, scale)
linshi <- t(linshi)
colnames(linshi) <- colnames(NPCData)
linshi[linshi > 1] <- 1
linshi[linshi < -1] <- -1
group <- sampleinfo$sample.time[sampleinfo$X %in% colnames(NPCData)]
group <- factor(group, levels=c("M0","M1","M2","End","End_SVZ"),
                labels=c("T0","T1","T2","T3","T4"))
annotation_col <- data.frame(row.names=colnames(NPCData), condition=group)

annotation_color <- list(
  condition=c(T0="#FFFF96", T1="#FFC864", T2="#E7B800", T3="#FF9632", T4="#FF6400"),
  GeneGroup=c(T0="#FFFF96", T1="#FFC864", T2="#E7B800", T3="#FF9632", T4="#FF6400")
)

col_order <- order(annotation_col$condition)
annotation_col_sorted <- annotation_col[col_order, , drop=F]
gene_annotation <- data.frame(GeneGroup=top20markers$cluster, row.names=top20markers$gene)
gene_annotation$GeneGroup <- factor(gene_annotation$GeneGroup, levels=c("T0","T1","T2","T3","T4"))

time_order <- c("T0", "T1", "T2", "T3", "T4")
ordered_genes <- unlist(lapply(time_order, function(tp) {
  rownames(gene_annotation)[gene_annotation$GeneGroup == tp]
}))

gene_annotation_sorted <- gene_annotation[ordered_genes, , drop=F]
data_sorted <- linshi[ordered_genes, col_order]
annotation_row <- gene_annotation_sorted
pdf("pheatmap_NSC_markers.pdf", width=6, height=10)
pheatmap(data_sorted, border_color="NA",
         color=colorRampPalette(c("#587AA7", "white", "#F44336"))(100),
         scale="none", annotation_col=annotation_col_sorted,
         annotation_row=annotation_row, annotation_names_row=F,
         annotation_names_col=F, annotation_colors=annotation_color,
         cluster_rows=F, cluster_cols=F, legend=T,
         show_rownames=T, show_colnames=F, fontsize=10)
dev.off()

######################
T1 <- FindMarkers(NSC, ident.1="T1", ident.2="T0", group.by="stimgroup", 
                      min.pct=0.25, logfc.threshold=0.25, only.pos=F)
T2 <- FindMarkers(NSC, ident.1="T2", ident.2="T0", group.by="stimgroup", 
                      min.pct=0.25, logfc.threshold=0.25, only.pos=F)
T3 <- FindMarkers(NSC, ident.1="T3", ident.2="T0", group.by="stimgroup", 
                      min.pct=0.25, logfc.threshold=0.25, only.pos=F)
T4 <- FindMarkers(NSC, ident.1="T4", ident.2="T0", group.by="stimgroup", 
                      min.pct=0.25, logfc.threshold=0.25, only.pos=F)
add_direction <- function(df, fdr_cutoff = 0.05, logfc_cutoff = 1) {
  df$direction <- "NOT"
  df$direction[df$p_val_adj < fdr_cutoff & df$avg_log2FC > logfc_cutoff] <- "UP"
  df$direction[df$p_val_adj < fdr_cutoff & df$avg_log2FC < -logfc_cutoff] <- "DOWN"
  return(df)
}
T1 <- add_direction(T1, fdr_cutoff = 0.05)
T2 <- add_direction(T2, fdr_cutoff = 0.05)
T3 <- add_direction(T3, fdr_cutoff = 0.05)
T4 <- add_direction(T4, fdr_cutoff = 0.05)
write.table(T1,'deg_T1.txt',sep = '\t',quote = FALSE)
write.table(T2,'deg_T2.txt',sep = '\t',quote = FALSE)
write.table(T3,'deg_T3.txt',sep = '\t',quote = FALSE)
write.table(T4,'deg_T4.txt',sep = '\t',quote = FALSE)
prepare_plot_data <- function(df, label) {
  df.new <- df[df$direction != "NOT", ]
  df.new$label <- label
  df.new$ptname <- rownames(df.new)
  df.new$direction <- paste0(label, "-", df.new$direction)
  return(df.new)
}

T1.plot <- prepare_plot_data(T1, "T1")
T2.plot <- prepare_plot_data(T2, "T2")
T3.plot <- prepare_plot_data(T3, "T3")
T4.plot <- prepare_plot_data(T4, "T4")
merge.plot <- rbind(T1.plot, T2.plot, T3.plot, T4.plot)
merge.plot_filtered <- merge.plot %>% filter(p_val_adj < 0.01)
top5_points <- merge.plot_filtered %>%
  group_by(label) %>%
  top_n(5, avg_log2FC) %>%
  mutate(type="top") %>%
  ungroup()

bottom5_points <- merge.plot_filtered %>%
  group_by(label) %>%
  top_n(-5, avg_log2FC) %>%
  mutate(type="bottom") %>%
  ungroup()

extreme_points <- bind_rows(top5_points, bottom5_points)
dbar <- merge.plot_filtered %>%
  group_by(label) %>%
  summarise(avg_log2FC_min=min(avg_log2FC), avg_log2FC_max=max(avg_log2FC))
pdf("diff_gene_T0_vs_all.pdf", width=8, height=6)#Figure 2c
ggplot() +
  geom_col(data=dbar, aes(x=label, y=avg_log2FC_min),
           fill="#FFF9C4", alpha=0.6, width=0.7) +
  geom_col(data=dbar, aes(x=label, y=avg_log2FC_max),
           fill="#FFF9C4", alpha=0.6, width=0.7) +
  geom_jitter(data=merge.plot_filtered, aes(x=label, y=avg_log2FC, color=direction),
              width=0.35, size=0.1) +
  geom_tile(data=dbar, aes(x=label, y=0, fill=label),
            height=0.5, color="black", alpha=0.6, show.legend=F) +
  geom_text_repel(data=extreme_points, aes(x=label, y=avg_log2FC, label=ptname),
                  size=2.5, box.padding=0.5, point.padding=0.3,
                  max.overlaps=Inf, segment.color="grey50") +
  scale_fill_npg() +
  scale_color_manual(values=c('T1-UP'='#E63946', 'T2-UP'='#457B9D',
                              'T3-UP'='#2A9D8F', 'T4-UP'='#F4A261',
                              'T1-DOWN'='#800080', 'T2-DOWN'='#FF00FF',
                              'T3-DOWN'='#264653', 'T4-DOWN'='#06D6A0')) +
  ylab("log2(Fold-change)") + xlab("Time Point") +
  theme(axis.text.x=element_text(angle=45, hjust=1),
        panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        panel.background=element_blank(), axis.line=element_line(color="black"))
dev.off()

########################Figure 2F
time_points <- c("T1", "T2", "T3", "T4")
avg_log2FC_cutoff <- 1
pval_cutoff <- 0.05
for (tp in time_points) {
  df <- get(tp)
  all_degs <- rownames(df)[abs(df$avg_log2FC) > avg_log2FC_cutoff & df$p_val_adj < pval_cutoff]
  gene_ids <- bitr(all_degs, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)$ENTREZID
  options(timeout = 600)
  kegg_res <- enrichKEGG(gene = gene_ids, organism = 'hsa', pvalueCutoff = 0.05)
  kegg_df <- as.data.frame(kegg_res)
  top_kegg <- kegg_df %>% 
    arrange(desc(Count)) %>% 
    head(15)
 
  p2 <- ggplot(top_kegg, aes(x = Count, y = reorder(Description, Count))) +
    geom_col(fill = "#4A81AD", width = 0.8, color = "white", linewidth = 0.2) +
    scale_x_continuous(position = "bottom", expand = expansion(mult = c(0, 0.1))) + 
    labs(x = "Gene count", y = NULL) +
    theme_bw() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          axis.text.y = element_text(size = 10, color = "black"),
          axis.text.x = element_text(size = 10, color = "black"),
          axis.title.x = element_text(size = 12, vjust = 1)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 50))
  pdf(paste0('KEGG_DEGs_', tp, '.pdf'), width = 8, height = 6)
  print(p2)
  dev.off()
}
########################Figure 2E
mfuzz_plot_defined=function (eset, cl, mfrow = c(1, 1), colo, min.mem = 0, time.labels, 
                             time.points, ylim.set = c(0, 0), xlab = "Time", ylab = "Expression changes", 
                             x11 = TRUE, ax.col = "black", bg = "white", col.axis = "black", 
                             col.lab = "black", col.main = "black", col.sub = "black", 
                             col = "black", centre = FALSE, centre.col = "black", centre.lwd = 2, 
                             Xwidth = 5, Xheight = 5, single = FALSE, line.lwd=1,...) 
{
  clusterindex <- cl[[3]]
  memship <- cl[[4]]
  memship[memship < min.mem] <- -1
  colorindex <- integer(dim(exprs(eset))[[1]])
  if (missing(colo)) {
    colo <- c("#FF0000", "#FF1800", "#FF3000", "#FF4800", 
              "#FF6000", "#FF7800", "#FF8F00", "#FFA700", "#FFBF00", 
              "#FFD700", "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00", 
              "#AFFF00", "#97FF00", "#80FF00", "#68FF00", "#50FF00", 
              "#38FF00", "#20FF00", "#08FF00", "#00FF10", "#00FF28", 
              "#00FF40", "#00FF58", "#00FF70", "#00FF87", "#00FF9F", 
              "#00FFB7", "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", 
              "#00CFFF", "#00B7FF", "#009FFF", "#0087FF", "#0070FF", 
              "#0058FF", "#0040FF", "#0028FF", "#0010FF", "#0800FF", 
              "#2000FF", "#3800FF", "#5000FF", "#6800FF", "#8000FF", 
              "#9700FF", "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", 
              "#FF00EF", "#FF00D7", "#FF00BF", "#FF00A7", "#FF008F", 
              "#FF0078", "#FF0060", "#FF0048", "#FF0030", "#FF0018")
  }
  else {
    if (all(colo == "fancy")) {
      fancy.red <- c(c(0:255), rep(255, length(c(255:0))), 
                     c(255:150))
      colo <- rgb(b = fancy.blue/255, g = fancy.green/255, 
                  r = fancy.red/255)
    }
  }
  colorseq <- seq(0, 1, length = length(colo))
  for (j in 1:dim(cl[[1]])[[1]]) {
    if (single) 
      j <- single
    tmp <- exprs(eset)[clusterindex == j, , drop = FALSE]
    tmpmem <- memship[clusterindex == j, j]
    if (((j - 1)%%(mfrow[1] * mfrow[2])) == 0 | single) {
      if (x11) 
        X11(width = Xwidth, height = Xheight)
      if (sum(clusterindex == j) == 0) {
        ymin <- -1
        ymax <- +1
      }
      else {
        ymin <- min(tmp)
        ymax <- max(tmp)
      }
      if (sum(ylim.set == c(0, 0)) == 2) {
        ylim <- c(ymin, ymax)
      }
      else {
        ylim <- ylim.set
      }
      if (!is.na(sum(mfrow))) {
        par(mfrow = mfrow, bg = bg, col.axis = col.axis, 
            col.lab = col.lab, col.main = col.main, col.sub = col.sub, 
            col = col)
      }
      else {
        par(bg = bg, col.axis = col.axis, col.lab = col.lab, 
            col.main = col.main, col.sub = col.sub, col = col)
      }
      xlim.tmp <- c(1, dim(exprs(eset))[[2]])
      if (!(missing(time.points))) 
        xlim.tmp <- c(min(time.points), max(time.points))
      plot.default(x = NA, xlim = xlim.tmp, ylim = ylim, 
                   xlab = xlab, ylab = ylab, main = paste("Cluster", 
                                                          j), axes = FALSE, ...)
      if (missing(time.labels) && missing(time.points)) {
        axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]), 
             col = ax.col, ...)
        axis(2, col = ax.col, ...)
      }
      if (missing(time.labels) && !(missing(time.points))) {
        axis(1, time.points, 1:length(time.points), time.points, 
             col = ax.col, ...)
        axis(2, col = ax.col, ...)
      }
      if (missing(time.points) & !(missing(time.labels))) {
        axis(1, 1:dim(exprs(eset))[[2]], time.labels,
             col = ax.col, ...)
        axis(2, col = ax.col, ...)
      }
      if (!(missing(time.points)) & !(missing(time.labels))) {
        axis(1, time.points, time.labels, col = ax.col, 
             ...)
        axis(2, col = ax.col, ...)
      }
    }
    else {
      if (sum(clusterindex == j) == 0) {
        ymin <- -1
        ymax <- +1
      }
      else {
        ymin <- min(tmp)
        ymax <- max(tmp)
      }
      if (sum(ylim.set == c(0, 0)) == 2) {
        ylim <- c(ymin, ymax)
      }
      else {
        ylim <- ylim.set
      }
      xlim.tmp <- c(1, dim(exprs(eset))[[2]])
      if (!(missing(time.points)))
        xlim.tmp <- c(min(time.points), max(time.points))
      plot.default(x = NA, xlim = xlim.tmp, ylim = ylim, 
                   xlab = xlab, ylab = ylab, main = paste("Cluster", 
                                                          j), axes = FALSE, ...)
      if (missing(time.labels) && missing(time.points)) {
        axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]), 
             col = ax.col, ...)
        axis(2, col = ax.col, ...)
      }
      if (missing(time.labels) && !(missing(time.points))) {
        axis(1, time.points, 1:length(time.points), time.points, 
             col = ax.col, ...)
        axis(2, col = ax.col, ...)
      }
      if (missing(time.points) & !(missing(time.labels))) {
        axis(1, 1:dim(exprs(eset))[[2]], time.labels, 
             col = ax.col, ...)
        axis(2, col = ax.col, ...)
      }
      if (!(missing(time.points)) & !(missing(time.labels))) {
        axis(1, time.points, time.labels, col = ax.col, 
             ...)
        axis(2, col = ax.col, ...)
      }
    }
    if (length(tmpmem) > 0) {
      for(jj in 1:max(0, (length(colorseq) - 1))) {
        tmpcol <- (tmpmem >= colorseq[jj] & tmpmem <= 
                     colorseq[jj + 1])
        if (sum(tmpcol) > 0) {
          tmpind <- which(tmpcol)
          for (k in 1:length(tmpind)) {
            if (missing(time.points)) {
              lines(tmp[tmpind[k], ], col = colo[jj],lwd = line.lwd)
            }
            else lines(time.points, tmp[tmpind[k], ], 
                       col = colo[jj],lwd = line.lwd)
          }
        }
      }
    }
    if (centre) {
      lines(cl[[1]][j, ], col = centre.col, lwd = centre.lwd)
    }
    if (single) 
      return()
  }
}
marker_exp=AggregateExpression(NSC, group.by = "stimgroup", 
                               assays = "RNA", 
                               slot = "data")
diff=as.matrix(read.table("deg_T1.txt",header=T,sep="\t"))
diff1=as.matrix(read.table("deg_T2.txt",header=T,sep="\t"))
diff2=as.matrix(read.table("deg_T3.txt",header=T,sep="\t"))
diff3=as.matrix(read.table("deg_T4.txt",header=T,sep="\t"))
#################|logFC|>0.25 p<0.01
all=rbind(diff,diff1,diff2,diff3)
temp=all[which(as.numeric(all[,2])>=0.25),] 
up=temp[which(as.numeric(temp[,5])<= 0.01),]
temp=all[which(as.numeric(all[,2])<= -0.25),] 
down=temp[which(as.numeric(temp[,5])<= 0.01),]
gene=unique(c(rownames(up),rownames(down)))
data_4=as.matrix(marker_exp$RNA[gene,])
eset <- new("ExpressionSet",exprs = data_4)
tmp <- filter.std(eset,min.std=0) 
data.s <- standardise(tmp) 
m1 <- mestimate(data.s)
set.seed(2022101701)
cl_6 <- mfuzz(data.s,c=6,m=m1)
save(cl_6,file="Mfuzz_cl_6.Rdata")
gene_filter=c()
cluster_filter=c()
for(i in 1:dim(cl_6$membership)[1])
{
  if(max(cl_6$membership[i,])>0.6){
    gene_filter=c(gene_filter,rownames(cl_6$membership)[i])
    cluster_filter=c(cluster_filter,cl_6$cluster[i])
  }
}
peusodtime_gene=cbind(gene_filter,cluster_filter)#816
gene_cluster <- cbind(cl_6$cluster, data_4)
colnames(gene_cluster)[1] <- 'cluster'
plot_line=exprs(data.s)
plot_line=cbind(plot_line,cl_6[[3]])
colnames(plot_line)[5]="cluster"
color=c("#E3DCA3","#CADABD","#9AA193","#A1B2CC","#CB9C7A","#8A7B93")
for(i in 1:6)
{
  outfile=paste("deg_6_",i,".pdf",sep="")
  color1=rep(color[i],64)
  pdf(outfile)
  mfuzz_plot_defined(data.s,cl=cl_6,mfrow=c(1,1),single=i,
                     centre=T,line.lwd=2,colo=color1,
                     time.labels=c("T0","T1","T2","T3","T4"),
                     x11 = FALSE,centre.lwd=3,min.mem =0.75)
  dev.off()
}
#######################Figure S1
load('Mfuzz_cl_6.Rdata')
cluster_gene <- data.frame(gene=rownames(cl_6$membership),cluster=cl_6$cluster,mem=apply(cl_6$membership,1,max))
for (cl in sort(unique(cluster_gene$cluster))) {
  all_degs <- cluster_gene$gene[cluster_gene$cluster==cl]
  if(length(all_degs)<5) next
  go_res <- enrichGO(gene=all_degs,keyType="SYMBOL",OrgDb=org.Hs.eg.db,pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",qvalueCutoff=0.2,ont="BP",readable=T)
  GO_df <- as.data.frame(go_res)
  write.csv(GO_df,paste0("GO_cluster_",cl,".csv"),row.names=F)
  if(nrow(GO_df)>0){
    top_GO <- GO_df %>% arrange(desc(Count)) %>% head(10)
    p <- ggplot(top_GO,aes(x=Count,y=reorder(Description,Count))) +
      geom_col(fill="#4A81AD",width=0.8,color="white",linewidth=0.2) +
      scale_x_continuous(position="bottom",expand=expansion(mult=c(0,0.1))) +
      labs(x="Gene count",y=NULL) +
      theme_bw() +
      theme(panel.grid.major.y=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border=element_rect(color="black",fill=NA),
            axis.text.y=element_text(size=10,color="black"),
            axis.text.x=element_text(size=10,color="black"),
            axis.title.x=element_text(size=12,vjust=1)) +
      scale_y_discrete(labels=function(x) str_wrap(x,width=50))
    pdf(paste0("GO_BP_cluster_",cl,".pdf"),width=8,height=6)
    print(p)
    dev.off()
  }
}
