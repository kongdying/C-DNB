library(dplyr)
library(reshape2)
library(psych)
library(Seurat)
library(igraph)
library(parallel)
library(ggplot2)
library(tidyr)
library(purrr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(WGCNA)


outdir= "./step2"
dir.create(outdir)
setwd(outdir)
load('../step1/NSC.rds')
##########figure 3A-B
expression_matrix <- as(as.matrix(NPC@assays$RNA@layers$counts), 'sparseMatrix')
cell_metadata <- new('AnnotatedDataFrame',data=NPC@meta.data)# 132510
gene_annotation <- new('AnnotatedDataFrame',data=data.frame(gene_short_name = row.names(NPC), row.names = row.names(NPC)))
monocle_cds <- newCellDataSet(expression_matrix,
                              phenoData =cell_metadata,
                              featureData = gene_annotation)
cds <- monocle_cds
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)  
diff_test_res <- differentialGeneTest(cds,fullModelFormulaStr = "~ time")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))# 8234
cds <- setOrderingFilter(cds, ordering_genes)  
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2,method = 'DDRTree')
cds <- orderCells(cds)
pdf("pseudutime.state.pre.order.pdf")
plot_cell_trajectory(cds) 
dev.off()
pdf("pseudutime.time.pre.order.pdf")
plot_cell_trajectory(cds, color_by = "time",show_branch_points=F) + scale_color_manual(breaks = c("M0", "M1", "M2","End","End_SVZ"), values=c("#8491B4","#B5D6FD","#FDC086","#B0CBA4","#E88BC5")) + theme(legend.position = "right")
dev.off()
cds <- orderCells(cds,root_state=1)
pdf("pseudutime.Pseudotime.pdf")
plot_cell_trajectory(cds, color_by = "Pseudotime")+scale_color_gradientn(colours=c("#440154","#3C508B","#228E8D","#53C569","#EEE51C")) 
dev.off()
##########Figure 3J-L
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
diff_test <- differentialGeneTest(mycds[disp.genes,], cores = 4, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test, qval < 1e-04))
write.csv(sig_gene_names, "pseudotime/sig_gene_names.csv")
T1_dnb<-read.csv('time_1.csv')
T3_dnb<-read.csv('time_3.csv')
T1_pseudnb<-intersect(T1_dnb[[1]],sig_gene_names)
T3_pseudnb<-intersect(T3_dnb[[1]],sig_gene_names)
pseudnb<-c(T1_pseudnb,T3_pseudnb)
p1 = plot_pseudotime_heatmap(mycds[T1_pseudnb,], num_clusters=5,
                             show_rownames=T,return_heatmap=T)
ggsave("p1.pdf", plot = p1, width = 5, height = 3)
p3 = plot_pseudotime_heatmap(mycds[T3_pseudnb,], num_clusters=5,
                             show_rownames=T,return_heatmap=T)
ggsave("p3.pdf", plot = p3, width = 5, height = 3)
####################HVG
nsc_list <- SplitObject(NSC, split.by = "stimgroup")
nsc_list <- lapply(nsc_list, function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  return(x)
})
for (tp in names(nsc_list)) {
  hvg_genes <- VariableFeatures(nsc_list[[tp]])
  expr_matrix<-GetAssayData(nsc_list[[tp]],assay = "RNA",layer = "data")
  write.table(as.matrix(expr_matrix), file = paste0("Data_", tp, ".txt"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
  write.table(hvg_genes, file = paste0("HVG_", tp, ".txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
gene_names<-rownames(NSC@assays$RNA@features@.Data)
####################
precompute_cor_matrix <- function(expr_matrix, cor_threshold) {
  expr_t <- t(expr_matrix)
  gene_vars <- apply(expr_t, 2, sd, na.rm = TRUE)
  zero_var_genes <- names(gene_vars)[gene_vars == 0 | is.na(gene_vars)]
  if (length(zero_var_genes) > 0) {
    expr_t <- expr_t[, !colnames(expr_t) %in% zero_var_genes, drop = FALSE]
  }
  if (ncol(expr_t) == 0) {
    return(matrix(NA, nrow = 0, ncol = 0))
  }
  cor_matrix <- WGCNA::cor(expr_t, method = "pearson", use = "pairwise.complete.obs")
  cor_matrix[abs(cor_matrix) < cor_threshold] <- NA
  edges_indices <- which(!is.na(cor_matrix) & upper.tri(cor_matrix), arr.ind = TRUE)
  edges <- data.frame(
    Gene1 = rownames(cor_matrix)[edges_indices[, 1]],
    Gene2 = colnames(cor_matrix)[edges_indices[, 2]],
    Correlation = cor_matrix[edges_indices]
  )
  filename <- paste0("GBM_cor_edges_", sprintf("%.1f", cor_threshold), ".txt")
  write.table(edges, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)
}

expression_matrix<-GetAssayData(NSC,assay = "RNA",layer = "data")
cor_threshold=0.2
precompute_cor_matrix(expression_matrix,cor_threshold)

network_genes <- read.table("../../data/PPI.tsv",sep = '\t',header = F,stringsAsFactors = F)
edges<-read.table(paste0("GBM_cor_edges_",cor_threshold,'.txt'),sep = "\t",header = T)
ppi_ref <- network_genes %>%
  mutate(key = ifelse(V1 < V2, paste0(V1, "_", V2), paste0(V2, "_", V1))) %>%
  distinct(key)
network_edges <- edges %>%
  mutate(key = ifelse(Gene1 < Gene2, paste0(Gene1, "_", Gene2), paste0(Gene2, "_", Gene1))) %>%
  filter(key %in% ppi_ref$key) %>%
  dplyr::select(-key)
filename <- paste0("net_cor_threshold_",cor_threshold, ".txt")
write.table(network_edges, file = filename, sep = "\t", quote = FALSE, row.names = FALSE)

network_gene_set <- unique(c(network_edges$Gene1, network_edges$Gene2))
length(network_gene_set)
########################################################
Subigraph <- function(x) {
  HVG <-  read.table(paste0('HVG_T', x,'.txt'), sep = "\t")
  net <- read.table(filename, header = T, sep = "\t")
  v1 <-  net$Gene1%in%HVG$V1 | net$Gene2%in%HVG$V1
  filtered_net <- net[v1, ]
  all_nodes <- unique(c(filtered_net$Gene1, filtered_net$Gene2))
  v2 <- (net$Gene1 %in% all_nodes) | (net$Gene2 %in% all_nodes)
  final_edges <- net[v2, ]
  g_final <- graph_from_data_frame(final_edges, directed = FALSE)
  components <- igraph::components(g_final)
  largest <- which.max(components$csize)
  largest_nodes <- V(g_final)$name[components$membership == largest]
  largest_edges <- final_edges[final_edges$Gene1%in% largest_nodes & final_edges$Gene2 %in% largest_nodes, ]
  write.table(largest_edges, file = paste0("largest_", x, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
}
list<-0:4
for (str in list) {
  Subigraph(str)
}

library(WGCNA)
CIc_score <- function(y) {
  module_metrics_df <- data.frame(
    Module = character(),
    length = integer(),
    secondary = integer(),
    SD_in = numeric(),
    PCC_in = numeric(),
    PCC_out = numeric(),
    PCC_ratio = numeric(),
    DNB_score = numeric(),
    time = numeric(),
    stringsAsFactors = FALSE
  )
  network_file <- paste0("largest_", y, ".txt")
  GPCM <- read.table(network_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  target_genes <- read.table(paste0('HVG_T', y, '.txt'), sep = "\t")
  target_genes <- target_genes$V1  
  expr <- read.table(paste0('Data_T', y,'.txt'), sep = "\t",header=T, row.names=1)
  for (x in target_genes) {
    s <- Sys.time()
    neighbors_v1 <- GPCM[GPCM$Gene1 == x, "Gene2"]
    neighbors_v2 <- GPCM[GPCM$Gene2 == x, "Gene1"]
    x1 <- unique(c(neighbors_v1, neighbors_v2))
    x1 <- x1[x1 %in% gene_names]
    if (length(x1) == 0) next  
    secondary_v1 <- GPCM[GPCM$Gene1 %in% x1, "Gene2"]
    secondary_v2 <- GPCM[GPCM$Gene2 %in% x1, "Gene1"]
    x2 <- unique(c(secondary_v1, secondary_v2))
    x2 <- x2[x2 %in% gene_names]
    x2 <- setdiff(x2, c(x, x1))
   
    if (length(x2) == 0) {
      cat("基因", x, "没有二级邻居，跳过计算\n")
      next
    }
    
    module_genes <- c(x, x1)
    valid_module_genes <- module_genes[module_genes %in% rownames(expr)]
    if (length(valid_module_genes) == 0) {
      cat("基因", x, "的模块基因不在表达数据中，跳过计算\n")
      next
    }
    module_expr <- expr[valid_module_genes, , drop = FALSE]
    sd_in <- apply(module_expr, 1, sd, na.rm = TRUE)
    SD_in <- mean(sd_in, na.rm = TRUE)
    
    if (nrow(module_expr) > 1) {
      module_expr_t <- t(module_expr)
      cor_matrix_in <- WGCNA::cor(module_expr_t, use = "pairwise.complete.obs")
      diag(cor_matrix_in) <- NA
      PCC_in <- mean(abs(cor_matrix_in), na.rm = TRUE)
    } else {
      PCC_in <- NA
      cat("基因", x, "的模块内只有一个基因，无法计算相关性\n")
      next
    }
    
    valid_x2_genes <- x2[x2 %in% rownames(expr)]
    if (length(valid_x2_genes) == 0) {
      cat("基因", x, "的二级邻居不在表达数据中，跳过计算\n")
      next
    }
    x2_expr <- expr[valid_x2_genes, , drop = FALSE]
    combined_expr <- rbind(module_expr, x2_expr)
    combined_expr_t <- t(combined_expr)
    cor_matrix_all <- WGCNA::cor(combined_expr_t, use = "pairwise.complete.obs")
    module_indices <- 1:nrow(module_expr)
    x2_indices <- (nrow(module_expr) + 1):(nrow(module_expr) + nrow(x2_expr))
    cross_cor_matrix <- cor_matrix_all[module_indices, x2_indices, drop = FALSE]
    PCC_out <- mean(abs(cross_cor_matrix), na.rm = TRUE)
    if (is.na(PCC_out) || PCC_out == 0) {
      cat("基因", x, "的模块间相关性无效，跳过计算\n")
      next
    }
    #PCCc
    cor_matrix_cells <- WGCNA::cor(module_expr, use = "pairwise.complete.obs")
    diag(cor_matrix_cells) <- NA
    PCC_c<-mean(abs(cor_matrix_cells),na.rm = T)
    PCC_ratio <- PCC_in / PCC_out
    DNB_score <-(SD_in/PCC_c) * PCC_ratio
    elapsed <- as.numeric(difftime(Sys.time(), s, units = "secs"))
    new_row <- data.frame(
      Module = x,
      length = length(x1) + 1,
      secondary = length(x2),
      SD_in = round(SD_in, 4),
      PCC_c=round(PCC_c, 4),
      PCC_in = round(PCC_in, 4),
      PCC_out = round(PCC_out, 4),
      PCC_ratio = round(PCC_ratio, 4),
      DNB_score = round(DNB_score, 4),
      time = round(elapsed, 2)
    )
    
    module_metrics_df <- rbind(module_metrics_df, new_row)
  }
  output_file <- paste0("CItime_", y,".csv")
  write.csv(module_metrics_df, output_file, row.names = FALSE)
  
  return(module_metrics_df)
  
}



cores <-6
cl <- makeCluster(cores)
clusterEvalQ(cl, {
  library(WGCNA)
})
clusterExport(cl, varlist = "gene_names", envir = .GlobalEnv)
results <- parLapply(cl, 0:4, CIc_score)
stopCluster(cl)

##################
ci_results <- data.frame(Time = factor(paste0("T", 0:4), levels = paste0("T", 0:4)), CI = numeric(5))
for (y in 0:4){
  dnb<-read.table(paste0( "CItime_", y, ".csv"),sep=',',header = T)
  dnb<-dnb[order(dnb$DNB_score,decreasing=T),]
  write.table(dnb$Module[1:50],paste0( "top50DNB_T", y, ".txt"),sep='\t')
  top50CI<-mean(dnb$DNB_score[1:50])
  #print(top50CI)
  ci_results$CI[y+1] <- mean(dnb$DNB_score[1:50])
}
p <- ggplot(ci_results, aes(x = Time, y = CI, group = 1)) +
  geom_line(color = "#C0392B", linewidth = 1.2) +
  #geom_point(color = "#C0392B", linewidth = 4, stroke = 1, shape = 21, fill = "white") +
  labs(x = "Time Point", y = "Average DNB Score (Top 50)", title = "C-DNB Score") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(linewidth = 12, color = "black"),
        axis.title = element_text(linewidth = 14, face = "bold"),
        plot.title = element_text(hjust = 0.5, linewidth = 16, face = "bold"),
        axis.line = element_line(color = "black")) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))
pdf('DNB Score.pdf', width = 8, height = 6)#Figure 3C
print(p)
dev.off()


DNBscore <- function(y2) {
  module_metrics_df <- data.frame(
    Module = character(),
    length = integer(),
    secondary = integer(),
    SD_in = numeric(),
    PCC_in = numeric(),
    PCC_out = numeric(),
    PCC_ratio = numeric(),
    DNB_score = numeric(),
    time = numeric(),
    stringsAsFactors = FALSE
  )
  network_file <- paste0("largest_", y, ".txt")
  GPCM <- read.table(network_file, header = TRUE, sep = "\t")
  target_genes <- read.table(paste0('HVG_T', y, '.txt'),  sep = "\t")
  target_genes <- target_genes$V1  
  expr <- read.table(paste0('Data_T', y2,'.txt'), sep = "\t",header=T, row.names=1)
  for (x in target_genes) {
    s <- Sys.time()
    neighbors_v1 <- GPCM[GPCM$Gene1 == x, "Gene2"]
    neighbors_v2 <- GPCM[GPCM$Gene2 == x, "Gene1"]
    x1 <- unique(c(neighbors_v1, neighbors_v2))
    x1 <- x1[x1 %in% gene_names]
    
    if (length(x1) == 0) next  
    secondary_v1 <- GPCM[GPCM$Gene1 %in% x1, "Gene2"]
    secondary_v2 <- GPCM[GPCM$Gene2 %in% x1, "Gene1"]
    x2 <- unique(c(secondary_v1, secondary_v2))
    x2 <- x2[x2 %in% gene_names]
    x2 <- setdiff(x2, c(x, x1))
    
    if (length(x2) == 0) {
      cat("基因", x, "没有二级邻居，跳过计算\n")
      next
    }
    module_genes <- c(x, x1)
    
    valid_module_genes <- module_genes[module_genes %in% rownames(expr)]
    if (length(valid_module_genes) == 0) {
      cat("基因", x, "的模块基因不在表达数据中，跳过计算")
      next
    }
    module_expr <- expr[valid_module_genes, , drop = FALSE]
    sd_in <- apply(module_expr, 1, sd, na.rm = TRUE)
    SD_in <- mean(sd_in, na.rm = TRUE)
    
    if (nrow(module_expr) > 1) {
      module_expr_t <- t(module_expr)
      cor_matrix_in <- WGCNA::cor(module_expr_t, use = "pairwise.complete.obs")
      diag(cor_matrix_in) <- NA
      PCC_in <- mean(abs(cor_matrix_in), na.rm = TRUE)
    } else {
      PCC_in <- NA
      cat("基因", x, "的模块内只有一个基因，无法计算相关性")
      next
    }
    
    valid_x2_genes <- x2[x2 %in% rownames(expr)]
    if (length(valid_x2_genes) == 0) {
      cat("基因", x, "的二级邻居不在表达数据中，跳过计算\n")
      next
    }
    
    x2_expr <- expr[valid_x2_genes, , drop = FALSE]
    combined_expr <- rbind(module_expr, x2_expr)
    combined_expr_t <- t(combined_expr)
    cor_matrix_all <- WGCNA::cor(combined_expr_t, use = "pairwise.complete.obs")
    module_indices <- 1:nrow(module_expr)
    x2_indices <- (nrow(module_expr) + 1):(nrow(module_expr) + nrow(x2_expr))
    cross_cor_matrix <- cor_matrix_all[module_indices, x2_indices, drop = FALSE]
    PCC_out <- mean(abs(cross_cor_matrix), na.rm = TRUE)
    if (is.na(PCC_out) || PCC_out == 0) {
      cat("基因", x, "的模块间相关性无效，跳过计算\n")
      next
    }
    #PCCc
    cor_matrix_cells <- WGCNA::cor(module_expr, use = "pairwise.complete.obs")
    diag(cor_matrix_cells) <- NA
    PCC_c<-mean(abs(cor_matrix_cells),na.rm = T)
    
    PCC_ratio <- PCC_in / PCC_out
    DNB_score <-(SD_in/PCC_c) * PCC_ratio
    elapsed <- as.numeric(difftime(Sys.time(), s, units = "secs"))
    new_row <- data.frame(
      Module = x,
      length = length(x1) + 1,
      secondary = length(x2),
      SD_in = round(SD_in, 4),
      PCC_c=round(PCC_c, 4),
      PCC_in = round(PCC_in, 4),
      PCC_out = round(PCC_out, 4),
      PCC_ratio = round(PCC_ratio, 4),
      DNB_score = round(DNB_score, 4),
      time = round(elapsed, 2)
    )
    
    module_metrics_df <- rbind(module_metrics_df, new_row)
    
  }
  output_file <- paste0(y,"dnbtime_", y2, ".csv")
  write.csv(module_metrics_df, output_file, row.names = FALSE)
  
  return(module_metrics_df)
}

y=1
cores=8
cl <- makeCluster(cores)
clusterEvalQ(cl, {
  library(WGCNA)
})
clusterExport(cl, varlist = c("gene_names", 'y'), envir = .GlobalEnv)
results <- parLapply(cl, 0:4, DNBscore)
stopCluster(cl)

y=3
cl <- makeCluster(cores)
clusterEvalQ(cl, {
  library(WGCNA)
})
clusterExport(cl, varlist = c("gene_names", 'y'), envir = .GlobalEnv)
results <- parLapply(cl, 0:4, DNBscore)
stopCluster(cl)

###########Figure 3G-H
library(plotly)
library(dplyr)
library(htmlwidgets)
t_list <- 0:4
time_labels <- c("T0", "T1", "T2", "T3", "T4")

get_data_mtx <- function(ref_tp) {
  v4 <- read.table(paste0("top50DNB_T", ref_tp, ".txt"), sep = "\t")
  colnames(v4) <- "Module"
  t_data <- data.frame(Module = v4$Module)
  for (y in t_list) {
    ref_file <- paste0("../step2/", ref_tp, "dnbtime_", y, ".csv")
    t1 <- read.table(ref_file, sep = ",", header = TRUE)
    t1[is.na(t1)] <- 0
    tmp_sub <- data.frame(
      Module = t1$Module,
      DNB_score = t1$DNB_score
    )
    colnames(tmp_sub)[2] <- paste0("t", y, "_score")
    t_data <- left_join(t_data, tmp_sub, by = "Module")
  }
  data_mtx <- t_data[, -1, drop = FALSE]
  rownames(data_mtx) <- t_data$Module
  data_mtx[is.na(data_mtx)] <- 0
  return(data_mtx)
}

get_cum_mtx <- function(mtx) {
  cum_mtx <- as.data.frame(apply(mtx, 2, function(x) cumsum(x) / seq_along(x)))
  rownames(cum_mtx) <- 1:nrow(cum_mtx)
  return(cum_mtx)
}
draw_dnb_3d <- function(data_mtx,t) {
  n_rows <- nrow(data_mtx)
  Z <- matrix(-0.0001, nrow = n_rows, ncol = length(t_list), dimnames = list(NULL, time_labels))
  v_idx <- (1:n_rows) 
  Z[v_idx, ] <- apply(data_mtx[1:n_rows,], 2, as.numeric)
  T=t+1
  p <- plot_ly(x = 1:n_rows, y = time_labels, z = t(Z), type = "surface",
               lighting = list(specular = 0, ambient = 1),
               colorscale = list(c(0, "#300070"), c(0.4, "#30a0f0"), c(0.6, "#b0f060"), c(0.8, "#f0f070"), c(1, "#f04000"))) %>%
    add_trace(
      x = 1:nrow(Z),
      y = rep(paste0("T",t), nrow(Z)),
      z = Z[,T ],
      type = "scatter3d",
      mode = "lines",
      line = list(width = 3, color = "red")
    ) %>%
    layout(showlegend = FALSE,scene = list(camera = list(projection = list(type = "orthographic")),aspectmode = "cube"
    ))
  p
}

data_mtx_t1 <- get_cum_mtx(get_data_mtx(1))
data_mtx_t3 <- get_cum_mtx(get_data_mtx(3))

# 绘制并保存
p_t1 <- draw_dnb_3d(data_mtx_t1,1)
p_t3 <- draw_dnb_3d(data_mtx_t3,3)

saveWidget(p_t1, "DNB_3D_T1nolegend.html")
saveWidget(p_t3, "DNB_3D_T3nolegend.html")

p_t1
p_t3
###################Figure 3E
DNB_T1<- read.table("top50DNB_T1.txt", sep = "\t" )
colnames(DNB_T1)='t1'
DNB_T3<- read.table("top50DNB_T3.txt", sep = "\t" )
colnames(DNB_T3)='t3'
DNB_PPIC<-cbind(DNB_T1,DNB_T3)
write.csv(DNB_PPIC, 'DNB_PPIC.csv')

sce_outdir <- "SCE_compare"
dir.create(sce_outdir, showWarnings = FALSE)

time_points <- 0:4
top_n <- 50

sce_results <- data.frame(
  Time = factor(paste0("T", time_points), levels = paste0("T", time_points)),
  SCE_top50 = NA
)

for (y in time_points) {
  tt <- paste0("T", y)
  expr <- read.table(paste0("step2/Data_T", y, ".txt"),sep = "\t",header = TRUE,row.names = 1,check.names = FALSE)
  expr <- as.matrix(expr)
  mode(expr) <- "numeric"
  expr <- expr[rowSums(expr, na.rm = TRUE) > 0, , drop = FALSE]
  
  GPCM <- read.table(paste0("step2/largest_", y, ".txt"),header = TRUE,sep = "\t",stringsAsFactors = FALSE,check.names = FALSE)
  colnames(GPCM)[1:2] <- c("Gene1", "Gene2")
  GPCM$Gene1 <- as.character(GPCM$Gene1)
  GPCM$Gene2 <- as.character(GPCM$Gene2)
  
  GPCM <- GPCM[GPCM$Gene1 %in% rownames(expr) &GPCM$Gene2 %in% rownames(expr),]
  
  ppi1 <- data.frame(center = GPCM$Gene1,neighbor = GPCM$Gene2,stringsAsFactors = FALSE)
  ppi2 <- data.frame(center = GPCM$Gene2,neighbor = GPCM$Gene1,stringsAsFactors = FALSE)
  ppi_edges <- unique(rbind(ppi1, ppi2))
  neighbor_list <- split(ppi_edges$neighbor, ppi_edges$center)
  
  target_genes <- read.table(paste0("step2/HVG_T", y, ".txt"),sep = "\t",stringsAsFactors = FALSE)$V1
  neighbor_list <- neighbor_list[names(neighbor_list) %in% target_genes]
  neighbor_count <- sapply(neighbor_list, length)
  neighbor_list <- neighbor_list[neighbor_count >= 2]
  score_t <- data.frame()
  
  for (center_gene in names(neighbor_list)) {
    neighbors <- neighbor_list[[center_gene]]
    neighbors <- intersect(neighbors, rownames(expr))
    if (length(neighbors) < 2) next
    center_expr <- as.numeric(expr[center_gene, ])
    neighbor_expr <- expr[neighbors, , drop = FALSE]
    sd_center <- sd(center_expr, na.rm = TRUE)
    if (is.na(sd_center) || sd_center == 0) next
    
    pcc_vec <- as.numeric(WGCNA::cor( x = matrix(center_expr, ncol = 1),y = t(neighbor_expr),use = "pairwise.complete.obs",method = "pearson"))
    names(pcc_vec) <- neighbors
    neighbor_mean <- rowMeans(neighbor_expr, na.rm = TRUE)
    weight <- abs(pcc_vec) * abs(neighbor_mean)
    weight[is.na(weight)] <- 0
    
    if (sum(weight) == 0) next
    prob <- weight / sum(weight)
    prob <- prob[prob > 0]
    S <- length(prob)
    if (S < 2) next
    entropy <- -sum(prob * log(prob)) / log(S)
    SCE <- entropy * sd_center
    
    tmp <- data.frame(
      Time = tt,
      Gene = center_gene,
      n_neighbors = length(neighbors),
      SD_center = sd_center,
      mean_abs_PCC = mean(abs(pcc_vec), na.rm = TRUE),
      mean_neighbor_expr = mean(neighbor_mean, na.rm = TRUE),
      entropy = entropy,
      SCE = SCE,
      stringsAsFactors = FALSE
    )
    
    score_t <- rbind(score_t, tmp)
  }
  
  score_t <- score_t[order(score_t$SCE, decreasing = TRUE), ]
  score_t$Rank <- 1:nrow(score_t)
  write.csv(
    score_t,
    file = file.path(sce_outdir, paste0("SCE_score_all_genes_", tt, ".csv")),
    row.names = FALSE
  )
}
time_points <- 0:4
top_percent <- 0.05

Ht_result <- data.frame(
  Time = factor(paste0("T", time_points), levels = paste0("T", time_points)),
  n_valid_genes = NA,
  n_top5_genes = NA,
  Ht_top5_sum = NA
)

for (y in time_points) {
  
  tt <- paste0("T", y)
  
  score_t <- read.csv(
    file.path(sce_outdir, paste0("SCE_score_all_genes_", tt, ".csv")),
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  score_t <- score_t[order(score_t$SCE, decreasing = TRUE), ]
  
  n_top <- ceiling(nrow(score_t) * top_percent)
  top_t <- score_t[1:n_top, ]
  
  Ht <- sum(top_t$SCE, na.rm = TRUE)
  
  Ht_result$n_valid_genes[Ht_result$Time == tt] <- nrow(score_t)
  Ht_result$n_top5_genes[Ht_result$Time == tt] <- n_top
  Ht_result$Ht_top5_sum[Ht_result$Time == tt] <- Ht
}

write.csv(
  Ht_result,
  file = file.path(sce_outdir, "SCE_Ht_top5percent_sum.csv"),
  row.names = FALSE
)
p_sce_ht <- ggplot(Ht_result, aes(x = Time, y = Ht_top5_sum, group = 1)) +
  geom_line(color = "#4D9A51", linewidth = 1.2) +
  #geom_point(color = "#C0392B", size = 3) +
  labs(
    x = "Time Point",
    y = "Global SCE Score Ht (Top 5% Sum)",
    title = "SCE-DNB Score"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.line = element_line(color = "black")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

pdf(
  file.path(sce_outdir, "SCE.pdf"),
  width = 8,
  height = 6
)
print(p_sce_ht)
dev.off()
###################Figure 3F
score_df <- data.frame(
  Method = c("C-DNB", "Tranditional DNB", "SCE-DNB"),
  T0 = c(2.761082, 1.862000, 24.47601),
  T1 = c(4.253434, 2.540134, 31.34131),
  T2 = c(3.496520, 2.166690, 26.66610),
  T3 = c(4.146434, 2.248834, 29.99829),
  T4 = c(3.221412, 1.863066, 28.61240)
)

rate_df <- data.frame(
  Method = score_df$Method,
  T1_growth_rate = (score_df$T1 - score_df$T0) / score_df$T0 * 100,
  T1_decline_rate = (score_df$T1 - score_df$T2) / score_df$T1 * 100,
  T3_growth_rate = (score_df$T3 - score_df$T2) / score_df$T2 * 100,
  T3_decline_rate = (score_df$T3 - score_df$T4) / score_df$T3 * 100
)

rate_df[, -1] <- round(rate_df[, -1], 2)

print(rate_df)

write.csv(
  rate_df,
  file = "SCE_compare/T1_T3_growth_decline_rate_three_methods.csv",
  row.names = FALSE
)
library(ggplot2)
library(dplyr)
library(tidyr)

rate_df <- read.csv(
  "SCE_compare/T1_T3_growth_decline_rate_three_methods.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)

plot_df <- rate_df %>%
  pivot_longer(
    cols = c(
      T1_growth_rate,
      T1_decline_rate,
      T3_growth_rate,
      T3_decline_rate
    ),
    names_to = "Metric",
    values_to = "Rate"
  )

plot_df$Method <- factor(
  plot_df$Method,
  levels = c("C-DNB","Tranditional DNB","SCE-DNB"),
  labels =  c("C-DNB","Tranditional DNB","SCE-DNB"))

plot_df$Metric <- factor(
  plot_df$Metric,
  levels = c(
    "T1_growth_rate",
    "T1_decline_rate",
    "T3_growth_rate",
    "T3_decline_rate"
  ),
  labels = c(
    "T1 Growth",
    "T1 Decline",
    "T3 Growth",
    "T3 Decline"
  )
)

plot_df$Label <- sprintf("%.2f%%", plot_df$Rate)

p <- ggplot(
  plot_df,
  aes(x = Metric, y = Rate, fill = Method)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = Label),
    position = position_dodge(width = 0.8),
    vjust = -0.4,
    size = 3.8,
    color = "black"
  ) +
  scale_fill_manual(
    values = c(
      "C-DNB" = "#C0392B",
      "Tranditional DNB" = "#2E5A88",
      "SCE-DNB" = "#4D9A51"
    )
  ) +
  labs(
    x = "Tipping-point rate metric",
    y = "Rate (%)",
    title = "Relative Changes Rates around T1 and T3"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.line = element_line(color = "black"),
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15))
  )

pdf("SCE_compare/T1_T3_growth_decline_rate_barplot.pdf")
print(p)
dev.off()

########################################
DNB_PPIC<-read.csv('DNB_PPIC.csv',row.names = 1)
#各时间点内部打分
load_dnb <- function(p) {
  df_list <- lapply(0:4, function(i) {
    d <- read.csv(paste0(p, "dnbtime_", i, ".csv"))[, c("Module", "DNB_score")]
    setNames(d, c("Gene", paste0("T", i)))
  })
  as.data.frame(reduce(df_list, inner_join, by = "Gene"))
}
group1 <- load_dnb("1")
group3 <- load_dnb("3")
plot_data_long <- group1 %>%
  pivot_longer(cols = starts_with("T"), names_to = "Time", values_to = "Score") %>%
  mutate(Time = factor(Time, levels = c("T0", "T1", "T2", "T3", "T4"))) %>%
  drop_na() 
subset_data <- plot_data_long %>% filter(Gene %in% DNB_PPIC$t1)
p <- ggplot(subset_data, aes(x = Time, y = Score)) +
  geom_line(aes(group = Gene, color = Gene), alpha = 0.3, linewidth = 0.6) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "black", linewidth = 1.5) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", color = "black", size = 3) +
  labs(x = "Time Point", y = "DNB Score", title = "Individual Gene Trends & Mean") +
  theme_bw() +
  theme(legend.position = "none")
pdf('T1DNB_Score.pdf', width = 8, height = 6)
print(p)
dev.off()

plot_data_long <- group3 %>%
  pivot_longer(cols = starts_with("T"), names_to = "Time", values_to = "Score") %>%
  mutate(Time = factor(Time, levels = c("T0", "T1", "T2", "T3", "T4"))) %>%
  drop_na() 
subset_data <- plot_data_long %>% filter(Gene %in% DNB_PPIC$t3)
p <- ggplot(subset_data, aes(x = Time, y = Score)) +
  geom_line(aes(group = Gene, color = Gene), alpha = 0.3, linewidth = 0.6) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "black", linewidth = 1.5) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", color = "black", size = 3) +
  labs(x = "Time Point", y = "DNB Score", title = "Individual Gene Trends & Mean") +
  theme_bw() +
  theme(legend.position = "none")
print(p)
pdf('T3DNB_Score.pdf', width = 8, height = 6)
print(p)
dev.off()
##########################################Figure 3K-M
T1<-DNB_PPIC$t1
T3<-DNB_PPIC$t3
time_points <- c("T1", "T3")
for (tp in time_points) {
  df <- get(tp)
  gene_ids <- bitr(df, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)$ENTREZID
  options(timeout = 1600)
  kegg_res <- enrichKEGG(gene = gene_ids, organism = 'hsa', pvalueCutoff = 0.5)
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
  pdf(paste0('KEGG_DNBs_', tp, '.pdf'))
  print(p2)
  dev.off()
}
############################################Figure 3I+Cytoscape
extract_dnb_neighbors <- function(tp) {
  dnb_core <- read.table(paste0("step2/top50DNB_T", tp, ".txt"), sep = "\t", header = TRUE, stringsAsFactors = T)
  net <- read.table(paste0("step2/largest_", tp, ".txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  neighbor_edges <- net[(net$Gene1 %in% dnb_core$x) | (net$Gene2 %in% dnb_core$x), ]
  nodes <- unique(c(neighbor_edges$Gene1, neighbor_edges$Gene2))
  #group <- load_dnb("1")
  #cyto_res<-group[group$Gene%in%nodes,]#DNB值
  rescale_01 <- function(x) (x - min(x, na.rm=T)) / (max(x, na.rm=T) - min(x, na.rm=T))
  exp_list <- lapply(0:4, function(i) {
    df <- read.table(paste0("step2/Data_T", i, ".txt"), header=T, row.names=1, sep="\t")
    rowMeans(df[intersect(nodes, rownames(df)), , drop=FALSE], na.rm=T)
  })
  exp_combined <- do.call(cbind, exp_list)
  colnames(exp_combined) <- paste0("Exp_T", 0:4)
  exp_final <- as.data.frame(t(apply(exp_combined, 1, rescale_01)))
  exp_final$Gene <- rownames(exp_final)
  exp_final <- exp_final[, c("Gene", paste0("Exp_T", 0:4))]
  write.table(exp_final, paste0("step2/DNB_cyto_res_T", tp, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
}
extract_dnb_neighbors(tp = 1)
extract_dnb_neighbors(tp = 3)

