####生存分析验证####
library(GSVA)
library(survival)
library(survminer)
library(ggplot2)
# --- 1. 导入数据 ---
getdata<-function(dataset){
  if(dataset=="TCGA_702_samples"){#FPKM
    data=read.table("TCGA_702_samples/TCGA_mRNAseq_702.txt",sep="\t",header=T,check.names=F)
    clinical=read.table("TCGA_702_samples/TCGA_mRNAseq_702_clinical.txt",sep="\t",header=T)#没有primary
    info=rbind(clinical[,7],clinical[,6])
    rownames(data)=data[,1]
    data=data[,-1]
    ###stage
    stage=clinical[,3]
    unique(stage)
    #data=log2(data+1)
  }
  #609 90
  #663  89
  if(dataset=="mRNA_seq_693_samples"){ #FPKM
    data=read.table("mRNA_seq_693_samples/CGGA.mRNAseq_693.RSEM-genes.20200506.txt",sep="\t",header=T)
    clinical=read.table("mRNA_seq_693_samples/CGGA.mRNAseq_693_clinical.20200506.txt",sep="\t",header=T)#Primary
    info=rbind(clinical[,8],clinical[,7])
    ###stage
    stage=clinical[,4]
    unique(stage)
    rownames(data)=data[,1]
    data=data[,-1]
    data=log2(data+1)
  }
  #316 88
  if(dataset=="mRNA_seq_325_samples")#FPKM
  {
    data=read.table("mRNA_seq_325_samples/CGGA.mRNAseq_325.RSEM-genes.20200506.txt",sep="\t",header=T)
    clinical=read.table("mRNA_seq_325_samples/CGGA.mRNAseq_325_clinical.20200506.txt",sep="\t",header=T)#Primary
    info=rbind(clinical[,8],clinical[,7])
    ###stage
    stage=clinical[,4]
    unique(stage)
    rownames(data)=data[,1]
    data=data[,-1]
    data=log2(data+1)
  }
  #287 89
  if(dataset=="array_301_samples") #microarray
  {
    data=read.table("array_301_samples/CGGA.mRNA_array_301_gene_level.20200506.txt",sep="\t",header=T)
    clinical=read.table("array_301_samples/CGGA.mRNA_array_301_clinical.20200506.txt",sep="\t",header=T)#Primary
    info=rbind(clinical[,9],clinical[,8])
    ###stage
    stage=clinical[,5]
    unique(stage)
    rownames(data)=data[,1]
    data=data[,-1]
  }
  if(dataset=="Rembrandt_475")#microarray
  {
    data=read.table("Rembrandt_microarray_475/Rembrandt_mRNA_array_475.txt",sep="\t",header=T)
    clinical=read.table("Rembrandt_microarray_475/Rembrandt_mRNA_array_475_clinical.txt",sep="\t",header=T)#Primary
    info=rbind(clinical[,7],clinical[,6])
    ###stage
    stage=clinical[,3]
    rownames(data)=data[,1]
    data=data[,-1]
    unique(stage)
  }
  train_x=data
  
  identical(clinical[,1],colnames(data)) 
  rownames(info)[1:2]=c("os_status","os_time")
  colnames(info)=colnames(data)
  filter=union(which(is.na(info[1,])),which(is.na(info[2,])))
  if(length(filter)>0){
    train_x=t(train_x[,-filter])
    info=t(info[,-filter])
    stage=stage[-filter]
  }
  names(stage)=rownames(train_x)
  return(list(train=train_x,info_x=info,stage=stage))
}
setwd("../data")
train_702=getdata("TCGA_702_samples")
train_693=getdata("mRNA_seq_693_samples")
train_325=getdata("mRNA_seq_325_samples")
train_301=getdata("array_301_samples")
train_475=getdata("Rembrandt_475")
setwd("../code")

t1 <- read.table('step2/top50DNB_T1.txt',sep = '\t',header = T)
t3 <-  read.table('step2/top50DNB_T3.txt',sep = '\t',header = T)
gene_sets_list1 <- list(t1 = t1$x)
gene_sets_list3 <- list(t3 = t3$x)

# --- 2. ssGSEA 分析与生存绘图 ---
run_ssgsea_survival <- function(train_obj, outname,dir_path="./step5/cox_survial_result/") {
  # getdata 返回的 train 是 [样本(行) x 基因(列)]
  # GSVA 需要 [基因(行) x 样本(列)]，所以需要转置 t()
  expr_mat <- t(train_obj$train)
  
  print(paste("Running ssGSEA for:", outname))
  param1 <- ssgseaParam(
    expr = as.matrix(expr_mat),
    geneSets = gene_sets_list1)
  param3 <- ssgseaParam(
    expr = as.matrix(expr_mat),
    geneSets = gene_sets_list3)
  ssgsea_res1 <- gsva(param1)
  ssgsea_res3 <- gsva(param3)
  t1_score <- ssgsea_res1["t1", ]
  t3_score <- ssgsea_res3["t3", ]
  
  info <- train_obj$info_x # 包含 os_status, os_time
  survdata <- as.data.frame(info)
  survdata$os_time <- as.numeric(as.character(survdata$os_time))
  survdata$os_status <- as.numeric(as.character(survdata$os_status))
  
  # 合并分数
  merge_data1 <- data.frame(time = survdata$os_time, 
                           event = survdata$os_status, 
                           t1_score = t1_score)
  merge_data3 <- data.frame(time = survdata$os_time, 
                            event = survdata$os_status, 
                            t3_score = t3_score)
  
  # 使用 surv_cutpoint 寻找最佳截断值
  print("Calculating optimal cutpoint...")
  res_cut1 <- surv_cutpoint(merge_data1, time = "time", event = "event", variables = c("t1_score"))
  res_cat1 <- surv_categorize(res_cut1)
  res_cut3 <- surv_cutpoint(merge_data3, time = "time", event = "event", variables = c("t3_score"))
  res_cat3 <- surv_categorize(res_cut3)
  # 拟合生存曲线
  fit1 <- surv_fit(Surv(time, event) ~ t1_score, data = res_cat1)
  fit3 <- surv_fit(Surv(time, event) ~ t3_score, data = res_cat3)
  
  # 绘图
  gg1 <- ggsurvplot(fit1, 
                    pval = TRUE, 
                    conf.int = TRUE,
                    pval.coord = c(max(merge_data1$time, na.rm=T)*0.1, 0.2), # 自动调整p值位置
                    palette = c("#A21017", "#002D63"), # 你的配色: High=红, Low=蓝
                    xlab = "Time(Days)",
                    risk.table = TRUE,
                    linetype = "strata",
                    ggtheme = theme_bw() + theme(text = element_text(size = 12)),
                    title = paste0(outname, " (ssGSEA t1 Score)")
  )
  gg3 <- ggsurvplot(fit3, 
                    pval = TRUE, 
                    conf.int = TRUE,
                    pval.coord = c(max(merge_data3$time, na.rm=T)*0.1, 0.2), # 自动调整p值位置
                    palette = c("#A21017", "#002D63"), # 你的配色: High=红, Low=蓝
                    xlab = "Time(Days)",
                    risk.table = TRUE,
                    linetype = "strata",
                    ggtheme = theme_bw() + theme(text = element_text(size = 12)),
                    title = paste0(outname, " (ssGSEA t3 Score)")
  )
  # 保存 PDF
  if(!dir.exists(dir_path)) dir.create(dir_path, recursive = T)
  outfile <- paste0(dir_path, outname, "_ssgsea_survival_t1.pdf")
  pdf(file = outfile, height = 6, width = 5) # 稍微宽一点以显示 risk table
  print(gg1, newpage = FALSE)
  dev.off()
  outfile <- paste0(dir_path, outname, "_ssgsea_survival_t3.pdf")
  pdf(file = outfile, height = 6, width = 5) # 稍微宽一点以显示 risk table
  print(gg3, newpage = FALSE)
  dev.off()
  
  # 返回打分结果以便后续使用
  return(list(merge_data1,merge_data3))
}

score_702 <- run_ssgsea_survival(train_702, "TCGA_702")
score_693 <- run_ssgsea_survival(train_693, "CGGA_693")
score_325 <- run_ssgsea_survival(train_325, "CGGA_325")
score_301 <- run_ssgsea_survival(train_301, "Array_301")
score_475 <- run_ssgsea_survival(train_475, "Rembrandt_475")

####生存分析验证-分级####
#除了Rembrandt_475均只有WHO II、WHO III、WHO IV分期，Rembrandt_475包含的WHO I分期仅2例，不参与
run_ssgsea_survival_by_stage <- function(train_obj,outname,dir_path="./step5/cox_survial_result_by_stage/"){
  table(train_obj$stage)
  for(i in unique(train_obj$stage)){
    if(is.na(i)||i=="WHO I"){next}
    keep<-train_obj$stage==i
    train_n<-list()
    train_n$info_x<-train_obj$info_x[keep,]
    train_n$train<-train_obj$train[keep,]
    score_obj<-run_ssgsea_survival(train_n, paste(outname,i),dir_path)
  }
}

run_ssgsea_survival_by_stage(train_702, "TCGA_702")
run_ssgsea_survival_by_stage(train_693, "CGGA_693")
run_ssgsea_survival_by_stage(train_325, "CGGA_325")
run_ssgsea_survival_by_stage(train_301, "Array_301")
run_ssgsea_survival_by_stage(train_475, "Rembrandt_475")

####ROC####
library(survival)
library(timeROC)
library(ggplot2)

draw_ROC<- function(train_obj,score_obj,outname,dir_path="./step5/ROC_result/") {
  #t1
  roc_obj <- timeROC(
    T = train_obj$info_x[,2],#time
    delta = train_obj$info_x[,1],#state
    marker = score_obj[[1]]$t1_score,
    cause = 1,
    times = c(365,3*365,5*365)
  )
  plot_df <- rbind(
    data.frame(x = roc_obj$FP[,1], y = roc_obj$TP[,1], time = paste0("1-year (AUC=",roc_obj$AUC[1],")")),
    data.frame(x = roc_obj$FP[,2], y = roc_obj$TP[,2], time = paste0("3-year (AUC=",roc_obj$AUC[2],")")),
    data.frame(x = roc_obj$FP[,3], y = roc_obj$TP[,3], time = paste0("5-year (AUC=",roc_obj$AUC[3],")"))
  )
  gg1<-ggplot(plot_df, aes(x = x, y = y, color = time)) +
    geom_step(linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") + #参考线
    labs(x = "1 - Specificity", y = "Sensitivity", color = "") +
    scale_color_manual(values = c("#E63946","#457B9D","#F39C12")) +
    theme_bw() +
    coord_equal()
  #t3
  roc_obj <- timeROC(
    T = train_obj$info_x[,2],#time
    delta = train_obj$info_x[,1],#state
    marker = score_obj[[2]]$t3_score,
    cause = 1,
    times = c(365,3*365,5*365)
  )
  plot_df <- rbind(
    data.frame(x = roc_obj$FP[,1], y = roc_obj$TP[,1], time = paste0("1-year (AUC=",roc_obj$AUC[1],")")),
    data.frame(x = roc_obj$FP[,2], y = roc_obj$TP[,2], time = paste0("3-year (AUC=",roc_obj$AUC[2],")")),
    data.frame(x = roc_obj$FP[,3], y = roc_obj$TP[,3], time = paste0("5-year (AUC=",roc_obj$AUC[3],")"))
  )
  gg3<-ggplot(plot_df, aes(x = x, y = y, color = time)) +
    geom_step(linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") + #参考线
    labs(x = "1 - Specificity", y = "Sensitivity", color = "") +
    scale_color_manual(values = c("#E63946","#457B9D","#F39C12")) +
    theme_bw() +
    coord_equal()
  #保存
  if(!dir.exists(dir_path)) dir.create(dir_path, recursive = T)
  outfile <- paste0(dir_path, outname, "_ROC_t1.pdf")
  pdf(file = outfile, height = 4, width = 6)
  print(gg1)
  dev.off()
  outfile <- paste0(dir_path, outname, "_ROC_t3.pdf")
  pdf(file = outfile, height = 4, width = 6) 
  print(gg3)
  dev.off()
}

draw_ROC(train_702,score_702,"TCGA_702")
draw_ROC(train_693,score_693,"CGGA_693")
draw_ROC(train_325,score_325,"CGGA_325")
draw_ROC(train_301,score_301,"Array_301")
draw_ROC(train_475,score_475,"Rembrandt_475")
