library(GSVA)
library(limma)
library(ggplot2)
library(ggpubr) 
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(survival)

outdir= "./step3"
dir.create(outdir)
setwd(outdir)

DNB_PPIC<-read.csv('../step2/DNB_PPIC.csv')
############TCGA_702
# GSVA 需要 [基因(行) x 样本(列)]，所以需要转置 t()
outdir= "./TCGA_702_samples"
dir.create(outdir)
setwd(outdir)
data=read.table("../../../data/TCGA_702_samples/TCGA_mRNAseq_702.txt",sep="\t",header=T,check.names=F)
clinical=read.table("../../../data/TCGA_702_samples/TCGA_mRNAseq_702_clinical.txt",sep="\t",header=T)#没有primary
rownames(data)=data[,1]
data=data[,-1]
gene_sets_list<-list(
  t1 = DNB_PPIC$t1,
  t3 = DNB_PPIC$t3
)
param <- ssgseaParam(
  expr = as.matrix(data),
  geneSets = gene_sets_list)
ssgsea_score <- gsva(param)
clin_var<-c("Histology","Grade","Gender","Age","IDH_mutation_status","X1p19q_codeletion_status")
# gsva(param) 返回的是矩阵 (Pathway x Sample)，需要转置为 (Sample x Pathway)
ssgsea_df <- as.data.frame(t(ssgsea_score))
ssgsea_df$SampleID <- rownames(ssgsea_df)
#将 clinical 的第一列名改为 "SampleID" 以便合并
colnames(clinical)[1] <- "SampleID"
plot_data <- merge(clinical, ssgsea_df, by = "SampleID")

# 处理年龄分组 (40岁前后)
plot_data$Age <- as.numeric(plot_data$Age)
plot_data$Age_Group <- ifelse(plot_data$Age >= 40, ">=40", "<40")

vars <- c("Histology","Grade", "Gender", "IDH_mutation_status", "X1p19q_codeletion_status", "Age_Group")
valid_vars <- intersect(vars, colnames(plot_data))

batch_save_boxplots <- function(data,target_pathways,clinical_vars) {
  draw_boxplot <- function(data, var_group, var_score) {
    sub_data <- data[!is.na(data[[var_group]]) & !is.na(data[[var_score]]), ]
    # 自动判断统计方法 (2组 Wilcox, >2组 Kruskal)
    n_groups <- length(unique(sub_data[[var_group]]))
    if(n_groups < 2) return(NULL) # 如果分组少于2，无法画图，跳过
    stat_method <- ifelse(n_groups > 2, "kruskal.test", "wilcox.test")
    # 绘图
    colors <- c(brewer.pal(8,"Set2"), brewer.pal(12,"Set3"))
    p <- ggplot(sub_data, aes(x = .data[[var_group]], y = .data[[var_score]], fill = .data[[var_group]])) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 1, alpha = 0.4) +
      theme_bw() +
      labs(x = "", y = "ssGSEA Score", title = paste(var_score, "by", var_group)) +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none",
            plot.title = element_text(hjust = 0.5, size = 10)) +
      scale_fill_manual(values = colors) +
      stat_compare_means(method = stat_method, label = "p.signif") # 显示显著性(*号)
    return(p)
  }
  target_pathways <- c("t1", "t3") 
  all_plots <- list()
  counter <- 1
  for (pathway in target_pathways) {
    for (clin in valid_vars) {
      p <- draw_boxplot(plot_data, clin, pathway)
      file_name <- paste0( "Boxplot_", pathway, "_", clin, ".pdf")
      n_groups <- length(ggplot_build(p)$data[[1]]$x)
      ggsave(filename = file_name, plot = p, width = 0.5+1*n_groups, height = 6)
      all_plots[[counter]] <- p
      counter <- counter + 1
    }
  }
}

my_plots <- batch_save_boxplots(plot_data,target_pathways,valid_vars)

###########独立预后
bioForest <- function(coxFile = NULL, forestFile = NULL, forestCol = NULL) {
  rt <- read.table(coxFile, header = T, sep = "\t", check.names = F, row.names = 1)
  gene <- rownames(rt)
  hr <- sprintf("%.3f", rt$"HR")
  hrLow  <- sprintf("%.3f", rt$"HR.95L")
  hrHigh <- sprintf("%.3f", rt$"HR.95H")
  Hazard.ratio <- paste0(hr, "(", hrLow, "-", hrHigh, ")")
  pVal <- ifelse(rt$pvalue < 0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  pdf(file = forestFile, width = 6, height = 4.3)
  n <- nrow(rt)
  nRow <- n + 1
  ylim <- c(1, nRow)
  layout(matrix(c(1, 2), nc = 2), width = c(3, 2.5))
  
  xlim = c(0, 3)
  par(mar = c(4, 2.5, 2, 1))
  plot(1, xlim = xlim, ylim = ylim, type = "n", axes = F, xlab = "", ylab = "")
  text.cex = 0.8
  text(0, n:1, gene, adj = 0, cex = text.cex)
  text(1.5 - 0.5 * 0.2, n:1, pVal, adj = 1, cex = text.cex)
  text(1.5 - 0.5 * 0.2, n + 1, 'pvalue', cex = text.cex, font = 2, adj = 1)
  text(3.1, n:1, Hazard.ratio, adj = 1, cex = text.cex)
  text(3.1, n + 1, 'Hazard ratio', cex = text.cex, font = 2, adj = 1)
  
  # 绘制森林图
  par(mar = c(4, 1, 2, 1), mgp = c(2, 0.5, 0))
  xlim = c(0, max(as.numeric(hrLow), as.numeric(hrHigh)))
  plot(1, xlim = xlim, ylim = ylim, type = "n", axes = F, ylab = "", xaxs = "i", xlab = "Hazard ratio")
  arrows(as.numeric(hrLow), n:1, as.numeric(hrHigh), n:1, angle = 90, code = 3, length = 0.05, col = "darkblue", lwd = 2.5)
  abline(v = 1, col = "black", lty = 2, lwd = 2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex = 1.5)
  axis(1)
  dev.off()
}
dulicox_custom <- function(input_data, dataset_name, target_score) {
  input_data$os_time <- as.numeric(input_data$os_time)
  input_data$os_status <- as.numeric(input_data$os_status)
  base_vars <- c("Gender","Age_Group", "Grade", "IDH_mutation_status", "X1p19q_codeletion_status","TCGA_subtypes","PRS_type","Histology","MGMTp_methylation_status")
  exist_vars <- intersect(base_vars, colnames(input_data))
  use_cols <- c("os_status", "os_time", exist_vars, target_score)
  rt <- input_data[, use_cols]
  if ("Gender" %in% colnames(rt)) {
    rt <- rt %>% mutate(Gender = case_when(
      Gender %in% c("Male", "male", "M") ~ 0,
      Gender %in% c("Female", "female", "F") ~ 1,
      TRUE ~ NA_real_
    ))
  }
  
  if ("IDH_mutation_status" %in% colnames(rt)) {
    rt <- rt %>% mutate(IDH_mutation_status = case_when(
      IDH_mutation_status %in% c("Mutant", "mutant") ~ 0,
      IDH_mutation_status %in% c("Wildtype", "WT", "wildtype") ~ 1,
      TRUE ~ NA_real_
    ))
  }
  
  if ("X1p19q_codeletion_status" %in% colnames(rt)) {
    rt <- rt %>% mutate(X1p19q_codeletion_status = case_when(
      X1p19q_codeletion_status %in% c("Non-codel", "non-codel") ~ 0,
      X1p19q_codeletion_status %in% c("Codel", "codel") ~ 1,
      TRUE ~ NA_real_
    ))
  }
  
  if ("Grade" %in% colnames(rt)) {
    # 转换为字符再处理，防止因子类型报错
    rt$Grade <- as.character(rt$Grade)
    rt <- rt %>% mutate(Grade = case_when(
      Grade %in% c("WHO I", "Grade I") ~ 1,
      Grade %in% c("WHO II", "Grade II") ~ 2,
      Grade %in% c("WHO III", "Grade III") ~ 3,
      Grade %in% c("WHO IV", "Grade IV", "GBM") ~ 4,
      TRUE ~ NA_real_ # 其他情况设为NA
    ))
  }
  
  if ("Age_Group" %in% colnames(rt)) {
    rt$Age_Group <- ifelse(rt$Age_Group %in% '<40',0,1) # 确保是数值
  }
  # TCGA Subtypes
  if ("TCGA_subtypes" %in% colnames(rt)) {
    rt$TCGA_subtypes <- as.character(rt$TCGA_subtypes)
    rt <- rt %>% mutate(TCGA_subtypes = case_when(
      TCGA_subtypes == "Proneural" ~ 1,
      TCGA_subtypes == "Neural" ~ 2,
      TCGA_subtypes == "Classical" ~ 3,
      TCGA_subtypes == "Mesenchymal" ~ 4,
      TRUE ~ NA_real_
    ))
  }
  
  # PRS Type
  if ("PRS_type" %in% colnames(rt)) {
    rt$PRS_type <- as.character(rt$PRS_type)
    rt <- rt %>% mutate(PRS_type = case_when(
      PRS_type == "Primary" ~ 1,
      PRS_type == "Recurrent" ~ 2,
      PRS_type == "Secondary" ~ 3,
      TRUE ~ NA_real_
    ))
  }
  
  # Histology
  if ("Histology" %in% colnames(rt)) {
    rt$Histology <- as.character(rt$Histology)
    # 这里需要根据你数据中实际的 Histology 名称来写
    # 这是一个比较通用的映射，你可以根据需要修改
    rt <- rt %>% mutate(Histology = case_when(
      grepl("Oligodendroglioma|O", Histology, ignore.case = T) & !grepl("Anaplastic", Histology, ignore.case = T) ~ 1,
      grepl("Astrocytoma|A", Histology, ignore.case = T) & !grepl("Anaplastic", Histology, ignore.case = T) ~ 2,
      grepl("Anaplastic", Histology, ignore.case = T) ~ 3,
      grepl("Glioblastoma|GBM", Histology, ignore.case = T) ~ 4,
      TRUE ~ NA_real_
    ))
  }
  
  # MGMT
  if ("MGMTp_methylation_status" %in% colnames(rt)) {
    rt$MGMTp_methylation_status <- as.character(rt$MGMTp_methylation_status)
    rt <- rt %>% mutate(MGMTp_methylation_status = case_when(
      grepl("Un-methylated|Unmethylated", MGMTp_methylation_status, ignore.case = T) ~ 1,
      grepl("Methylated", MGMTp_methylation_status, ignore.case = T) ~ 0,
      TRUE ~ NA_real_
    ))
  }
  rt <- na.omit(rt)
  
  # 单因素 Cox 分析
  uniTab <- data.frame()
  #  (临床变量 + target_score)
  analysis_vars <- colnames(rt)[3:ncol(rt)]
  for (i in analysis_vars) {
    # 构建公式
    form <- as.formula(paste0("Surv(os_time, os_status) ~ ", i))
    cox <- coxph(form, data = rt)
    coxSummary <- summary(cox)
    
    uniTab <- rbind(uniTab, cbind(
      id = i,
      HR = coxSummary$conf.int[, "exp(coef)"],
      HR.95L = coxSummary$conf.int[, "lower .95"],
      HR.95H = coxSummary$conf.int[, "upper .95"],
      pvalue = coxSummary$coefficients[, "Pr(>|z|)"]
    ))
  }
  
  outname_uni <- paste0(dataset_name, "_", target_score, "_uniCox.txt")
  write.table(uniTab, file = outname_uni, sep = "\t", row.names = F, quote = F)
  outpdf_uni <- paste0(dataset_name, "_", target_score, "_uniForest.pdf")
  if(nrow(uniTab) > 0) {
    bioForest(coxFile = outname_uni, forestFile = outpdf_uni, forestCol = "green")
  }
  
  
  # 多因素 Cox 分析
  multiCox <- coxph(Surv(os_time, os_status) ~ ., data = rt)
  multiCoxSum <- summary(multiCox)
  
  multiTab <- data.frame()
  multiTab <- cbind(
    HR = multiCoxSum$conf.int[, "exp(coef)"],
    HR.95L = multiCoxSum$conf.int[, "lower .95"],
    HR.95H = multiCoxSum$conf.int[, "upper .95"],
    pvalue = multiCoxSum$coefficients[, "Pr(>|z|)"]
  )
  multiTab <- cbind(id = row.names(multiTab), multiTab)
  
  outname_multi <- paste0(dataset_name, "_", target_score, "_multiCox.txt")
  write.table(multiTab, file = outname_multi, sep = "\t", row.names = F, quote = F)
  outpdf_multi <- paste0(dataset_name, "_", target_score, "_multiForest.pdf")
  
  if(nrow(multiTab) > 0){
    bioForest(coxFile = outname_multi, forestFile = outpdf_multi, forestCol = "red")
  }
}

colnames(plot_data)[which(colnames(plot_data) == "Censor")] <- "os_status"
colnames(plot_data)[which(colnames(plot_data) == "OS")] <- "os_time"

# 运行独立预后分析 
dulicox_custom(plot_data, "TCGA_702", "t1")
dulicox_custom(plot_data, "TCGA_702","t3")





###########################################################mRNA_seq_693_samples
outdir= "../mRNA_seq_693_samples"
dir.create(outdir)
setwd(outdir)
data=read.table("../../../data/mRNA_seq_693_samples/CGGA.mRNAseq_693.RSEM-genes.20200506.txt",sep="\t",header=T)
clinical=read.table("../../../data/mRNA_seq_693_samples/CGGA.mRNAseq_693_clinical.20200506.txt",sep="\t",header=T)#Primary
rownames(data)=data[,1]
data=data[,-1]
data=log2(data+1)
colnames(clinical)
param <- ssgseaParam(
  expr = as.matrix(data),
  geneSets = gene_sets_list)
ssgsea_score <- gsva(param)
# gsva(param) 返回的是矩阵 (Pathway x Sample)，需要转置为 (Sample x Pathway)
ssgsea_df <- as.data.frame(t(ssgsea_score))
ssgsea_df$SampleID <- rownames(ssgsea_df)
#将 clinical 的第一列名改为 "SampleID" 以便合并
colnames(clinical)[1] <- "SampleID"
plot_data <- merge(clinical, ssgsea_df, by = "SampleID")
# 处理年龄分组 (40岁前后)
plot_data$Age <- as.numeric(plot_data$Age)
plot_data$Age_Group <- ifelse(plot_data$Age >= 40, ">=40", "<40")
plot_data$Gender <- trimws(as.character(plot_data$Gender))
clin_var<-c("PRS_type","Histology","MGMTp_methylation_status","Grade","Gender","Age_Group","IDH_mutation_status","X1p19q_codeletion_status")
valid_vars <- intersect(clin_var, colnames(plot_data))
plots_693 <- batch_save_boxplots(plot_data,target_pathways,valid_vars)

colnames(plot_data)[which(colnames(plot_data) == "Censor..alive.0..dead.1.")] <- "os_status"
colnames(plot_data)[which(colnames(plot_data) == "OS")] <- "os_time"
dulicox_custom(plot_data, "mRNA_seq_693_samples", "t1")
dulicox_custom(plot_data, "mRNA_seq_693_samples", "t3")

###########################################################mRNA_seq_325_samples
outdir= "../mRNA_seq_325_samples"
dir.create(outdir)
setwd(outdir)
data=read.table("../../../data/mRNA_seq_325_samples/CGGA.mRNAseq_325.RSEM-genes.20200506.txt",sep="\t",header=T)
clinical=read.table("../../../data/mRNA_seq_325_samples/CGGA.mRNAseq_325_clinical.20200506.txt",sep="\t",header=T)#Primary
rownames(data)=data[,1]
data=data[,-1]
data=log2(data+1)
param <- ssgseaParam(
  expr = as.matrix(data),
  geneSets = gene_sets_list)
ssgsea_score <- gsva(param)
# gsva(param) 返回的是矩阵 (Pathway x Sample)，需要转置为 (Sample x Pathway)
ssgsea_df <- as.data.frame(t(ssgsea_score))
ssgsea_df$SampleID <- rownames(ssgsea_df)
colnames(clinical)[1] <- "SampleID"
plot_data <- merge(clinical, ssgsea_df, by = "SampleID")

# 处理年龄分组 (40岁前后)
plot_data$Age <- as.numeric(plot_data$Age)
plot_data$Age_Group <- ifelse(plot_data$Age >= 40, ">=40", "<40")
plot_data$Gender <- trimws(as.character(plot_data$Gender))
clin_var<-c("PRS_type","Histology","MGMTp_methylation_status","Grade","Gender","Age_Group","IDH_mutation_status","X1p19q_codeletion_status")
valid_vars <- intersect(clin_var, colnames(plot_data))
plots_325 <- batch_save_boxplots(plot_data,target_pathways,valid_vars)

colnames(plot_data)[which(colnames(plot_data) == "Censor..alive.0..dead.1.")] <- "os_status"
colnames(plot_data)[which(colnames(plot_data) == "OS")] <- "os_time"
dulicox_custom(plot_data, "mRNA_seq_325_samples", "t1")
dulicox_custom(plot_data, "mRNA_seq_325_samples", "t3")

###########################################################array_301_samples
outdir= "../array_301_samples"
dir.create(outdir)
setwd(outdir)
data=read.table("../../../data/array_301_samples/CGGA.mRNA_array_301_gene_level.20200506.txt",sep="\t",header=T)
clinical=read.table("../../../data/array_301_samples/CGGA.mRNA_array_301_clinical.20200506.txt",sep="\t",header=T)#Primary
colnames(clinical)[13]<-"X1p19q_codeletion_status"
rownames(data)=data[,1]
data=data[,-1]
param <- ssgseaParam(
  expr = as.matrix(data),
  geneSets = gene_sets_list)
ssgsea_score <- gsva(param)
# gsva(param) 返回的是矩阵 (Pathway x Sample)，需要转置为 (Sample x Pathway)
ssgsea_df <- as.data.frame(t(ssgsea_score))
ssgsea_df$SampleID <- rownames(ssgsea_df)
colnames(clinical)[1] <- "SampleID"
plot_data <- merge(clinical, ssgsea_df, by = "SampleID")

# 处理年龄分组 (40岁前后)
plot_data$Age <- as.numeric(plot_data$Age)
plot_data$Age_Group <- ifelse(plot_data$Age >= 40, ">=40", "<40")
plot_data$Gender <- trimws(as.character(plot_data$Gender))
clin_var<-c("TCGA_subtypes","PRS_type","Histology","MGMTp_methylation_status","Grade","Gender","Age_Group","IDH_mutation_status","X1p19q_codeletion_status")
valid_vars <- intersect(clin_var, colnames(plot_data))
plots_301 <- batch_save_boxplots(plot_data,target_pathways,valid_vars)

colnames(plot_data)[which(colnames(plot_data) == "Censor..alive.0..dead.1.")] <- "os_status"
colnames(plot_data)[which(colnames(plot_data) == "OS")] <- "os_time"
dulicox_custom(plot_data, "array_301_samples", "t1")
dulicox_custom(plot_data, "array_301_samples", "t3")


###########################################################Rembrandt_475
outdir= "../Rembrandt_microarray_475"
dir.create(outdir)
setwd(outdir)
data=read.table("../../../data/Rembrandt_microarray_475/Rembrandt_mRNA_array_475.txt",sep="\t",header=T)
clinical=read.table("../../../data/Rembrandt_microarray_475/Rembrandt_mRNA_array_475_clinical.txt",sep="\t",header=T)#Primary
rownames(data)=data[,1]
data=data[,-1]
param <- ssgseaParam(
  expr = as.matrix(data),
  geneSets = gene_sets_list)
ssgsea_score <- gsva(param)
# gsva(param) 返回的是矩阵 (Pathway x Sample)，需要转置为 (Sample x Pathway)
ssgsea_df <- as.data.frame(t(ssgsea_score))
ssgsea_df$SampleID <- rownames(ssgsea_df)
colnames(clinical)[1] <- "SampleID"
plot_data <- merge(clinical, ssgsea_df, by = "SampleID")
plot_data$Age_numeric <- as.numeric(sub("-.*", "", plot_data$Age))
plot_data$Age_Group <- ifelse(plot_data$Age_numeric >= 40, ">=40", "<40")
plot_data$Gender <- trimws(as.character(plot_data$Gender))
clin_var<-c("Histology","Grade","Gender","Age_Group","X1p19q_Codeletion_status")
valid_vars <- intersect(clin_var, colnames(plot_data))
plots_475 <- batch_save_boxplots(plot_data,target_pathways,valid_vars)
colnames(plot_data)[which(colnames(plot_data) == "Censor")] <- "os_status"
colnames(plot_data)[which(colnames(plot_data) == "OS")] <- "os_time"
dulicox_custom(plot_data, "Rembrandt_475", "t1")
dulicox_custom(plot_data, "Rembrandt_475", "t3")
