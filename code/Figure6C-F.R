setwd('./kdy/code')
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(dplyr)
outdir= "./step7"
dir.create(outdir)
setwd(outdir)
t1_up_tf<-read.table('../step6/result/t1_up_tf.txt',sep = '\t',header = F)
t3_up_tf<-read.table('../step6/result/t3_up_tf.txt',sep = '\t',header = F)
#Figure 6C
library(ggvenn)
T1_TF <- unique(t1_up_tf$V1)
T3_TF <- unique(t3_up_tf$V1)
gene_list <- list(
  T1 = T1_TF,
  T3 = T3_TF
)
p<-ggvenn(
  gene_list,
  fill_color = c("#FED976", "#74C476"),
  stroke_size = 0.5,
  set_name_size = 6,
  text_size = 1.3,
  show_elements = TRUE
)
pdf('TF_venn.pdf')
p
dev.off()

TF_T1 <- bitr(t1_up_tf$V1, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
write.table(TF_T1,'TF_T1.txt',sep = '\t')
TF_T3 <- bitr(t3_up_tf$V1, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
write.table(TF_T3,'TF_T3.txt',sep = '\t')

#run step7.2
#################Figure 6E
library(scTenifoldKnk)
library(Seurat)
library(ggplot2)
library(dplyr)
library(igraph)
library(ggrepel)
library(stringr)
set.seed(123)
load(NSC)
target_gene_name <- "NEK2"
count_matrix <- GetAssayData(NSC, layer = "data")
NSC <- FindVariableFeatures(object=NSC, selection.method="vst", nfeatures=8000)
high_variable_genes <- VariableFeatures(NSC)
input_data <- as.data.frame(count_matrix[unique(c(target_gene_name, high_variable_genes)),])
ko_analysis_result <- scTenifoldKnk(countMatrix = input_data,
                                    gKO = target_gene_name,      #需要敲除的基因
                                    qc = FALSE,                #是否进行QC
                                    qc_mtThreshold = 0.1,        #mt的阈值
                                    qc_minLSize = 1000,         #文库阈值(细胞测到的基因总数)
                                    nc_nNet = 10,              #子网络数量
                                    nc_nCells = 500            #每个网络中随机抽取的细胞数
)
saveRDS(ko_analysis_result,'kores.rds')

diff_regulation_df <- ko_analysis_result$diffRegulation
diff_regulation_df <- diff_regulation_df[diff_regulation_df$gene != target_gene_name, ]
significant_diff_table <- diff_regulation_df[diff_regulation_df$p.adj<0.05,]
write.table(diff_regulation_df, file="sigDiff_4000.txt", sep="\t", quote=F, row.names=F)

top20_diff_genes <- head(diff_regulation_df[order(-diff_regulation_df$FC), ], 20)
p1=ggplot(top20_diff_genes, aes(x=reorder(gene, FC), y=FC)) +
  geom_bar(stat='identity', fill='#5A9BD4') +
  coord_flip() + 
  labs(title="Top 20 Differentially Regulated Genes", x="Gene", y="FC") +
  theme_minimal() + theme(plot.title = element_text(hjust = 0.5))
pdf(file="barplot.pdf", width=6, height=5)
print(p1)
dev.off()
library(ggplot2)
library(dplyr)
library(ggrepel)
plot_data <- read.table("../../data/KO_Result_T1_SMC3.csv",sep = ",",header = TRUE,check.names = TRUE)
min_nonzero_p <- min(plot_data$p.value[plot_data$p.value > 0], na.rm = TRUE)

plot_data <- plot_data %>%
  mutate(
    pvalue_plot = ifelse(p.value == 0, min_nonzero_p / 10, p.value),
    log2FoldChange = log2(FC),
    negLog10P = -log10(pvalue_plot),
    group = case_when(
      p.adj < 0.05 & log2FoldChange > 1  ~ "Up",
      p.adj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "NotSig"
    ),
    label_name = ifelse(p.adj < 0.05, gene, "")
  )

p <- ggplot(plot_data, aes(x = log2FoldChange, y = negLog10P, color = group)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(
    values = c(
      "NotSig" = "grey70",
      "Up" = "red",
      "Down" = "blue"
    )
  ) +
  geom_text_repel(
    data = filter(plot_data, label_name != ""),
    aes(label = label_name),
    size = 2,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "black",
    segment.size = 0.4,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey40"
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "Virtual Knockout Volcano Plot: T1 SMC3",
    x = "log2(Fold Change)",
    y = "-log10(P value)",
    color = "Group"
  )+
   coord_cartesian(xlim = c(-50, 50))
pdf(file="KO_SMC3_vocanoplot.pdf")
p
dev.off()
##############################Figure 6
library(readxl)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
genelist<-plot_data$gene[plot_data$group=='Up']
ego <- enrichGO(gene          = genelist,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'SYMBOL',
                ont           = "BP",
                pAdjustMethod = "none",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 1,
                readable      = TRUE)
pdf(file="GO_SMC3.pdf")
print(barplot(ego, showCategory = 10))
dev.off()
write.csv(as.data.frame(ego), "All_Links_GOSMC3.csv", row.names = FALSE)

gene_convert <- bitr(genelist, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
options(timeout = 60000)
ekegg <- enrichKEGG(gene         = gene_convert$ENTREZID,
                    organism     = 'hsa', 
                    pvalueCutoff = 0.05)
ekegg@result <- ekegg@result[order(ekegg@result$Count, decreasing = TRUE), ]
p <- barplot(
  ekegg,
  color = "pvalue"
)
pdf('KO_function_result.pdf')
print(p)
dev.off()
write.csv(as.data.frame(ekegg), "All_Links_kegg_SMC3.csv", row.names = FALSE)

