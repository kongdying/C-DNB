setwd("./cmap/filter_t1/")
data = read.table("filter_pvalue_2.txt", sep="\t")
result_df <- data[which(data$V7 < 0.01 & data$V5 < 0), ]#
result_df <- result_df[order(result_df$V5, decreasing = FALSE), ]

clean_ids <- c()
for(full_id in result_df$V1){
  txt <- gsub("\\[|\\]|'", "", full_id)
  parts <- strsplit(txt, ":")[[1]]
  
  if(length(parts) >= 2){
    brd_id <- parts[2]
    
    if(grepl("BRD", brd_id) && grepl("-.*-", brd_id)){
      brd_parts <- strsplit(brd_id, "-")[[1]]
      brd_id <- paste(brd_parts[1], brd_parts[2], sep="-")
    }
    clean_ids <- c(clean_ids, brd_id)
  } else {
    clean_ids <- c(clean_ids, NA) 
  }
}

result_df$pert_id <- clean_ids
small_mole <- read.table("../../../../data/compoundinfo_beta.txt.used", sep="\t", quote="", header=T)
final_out <- merge(result_df, small_mole, by="pert_id", all.x=TRUE)
final_out <- final_out[order(final_out$V5, decreasing = FALSE), ]

##########################################Figure 7C
setwd("./cmap/filter_t3/")
data = read.table("filter_pvalue_2.txt", sep="\t")
result_df <- data[which(data$V7 < 0.01 & data$V5 < 0), ]#
result_df <- result_df[order(result_df$V5, decreasing = FALSE), ]

clean_ids <- c()
for(full_id in result_df$V1){
  txt <- gsub("\\[|\\]|'", "", full_id)
  parts <- strsplit(txt, ":")[[1]]
  
  if(length(parts) >= 2){
    brd_id <- parts[2]
    
    if(grepl("BRD", brd_id) && grepl("-.*-", brd_id)){
      brd_parts <- strsplit(brd_id, "-")[[1]]
      brd_id <- paste(brd_parts[1], brd_parts[2], sep="-")
    }
    clean_ids <- c(clean_ids, brd_id)
  } else {
    clean_ids <- c(clean_ids, NA)
  }
}

result_df$pert_id <- clean_ids
small_mole <- read.table("../../../../data/compoundinfo_beta.txt.used", sep="\t", quote="", header=T)
final_out <- merge(result_df, small_mole, by="pert_id", all.x=TRUE)
final_out <- final_out[order(final_out$V5, decreasing = FALSE), ]
output_cols <- c("pert_id", "cmap_name", "V5", "V6", "V1")

if("cmap_name" %in% colnames(final_out)){
  write_data <- final_out[, output_cols]
  colnames(write_data) <- c("Drug_ID", "Drug_Name", "NES", "P_Value", "Original_ID")
} else {
  write_data <- final_out
}
write_data<-write_data[!is.na(write_data$Drug_Name),]
write_data<-write_data[!duplicated(write_data$Drug_ID),]
write_data <- write_data %>% filter(!grepl("^BRD", Drug_Name))

write.table(write_data, "Final_drug_3.txt", 
            sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

drug_list <- data.frame(
  drug = write_data[1:10,2],
  score = write_data[1:10,3], 
  moa_t3 <- c(
    "CDK1/2/9 inhibitor",                # 细胞周期蛋白依赖性激酶1、2、9抑制剂
    "5-HT2C receptor agonist",           # 5-羟色胺2C受体激动剂
    "SARI",  # 5-羟色胺再摄取抑制剂及拮抗剂（SARI）
    "kinase inhibitor",                  # 激酶抑制剂
    "kinase inhibitor",                  # 激酶抑制剂
    "Topo II inhibitor",                # DNA拓扑异构酶II抑制剂
    "Gyrase/Topo IV inhibitor",  # 细菌DNA促旋酶及拓扑异构酶IV抑制剂
    "ChE inhibitor",                    # 胆碱酯酶抑制剂
    "kinase inhibitor",                 # 激酶抑制剂
    "PI3K/Akt inhibitors"# PI3K/Akt小分子抑制剂
  )
)

drug_list <- drug_list %>% arrange(score)
drug_list$drug <- factor(drug_list$drug, levels = drug_list$drug)
pill_width <- 0.6
pill_height <- 0.4 
pdf('Drug_function.pdf')
ggplot(drug_list, aes(x = drug, y = score)) +
  geom_point(aes(fill = moa_t3), size = 5, shape = 21, color = "black", stroke = 1) +
  geom_text(aes(label = sprintf("%.2f", score)), 
            vjust = -1.2, size = 3.5, fontface = "bold") +
  scale_fill_brewer(palette = "Set3") +
  labs(y = "NES", 
       x = "Drug Candidate", 
       fill = "Mechanism of Action") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  coord_cartesian(ylim =c(-2.32,-2.05) )

dev.off()


############################Figure 7B
setwd("D:/A055/大创/26年1月/药物/filter_1/")
write_data = read.table("Final_drug_1.txt", sep="\t",header=T,quote="", fill=T)
drug_list_t1 <- data.frame(
  drug = write_data[1:10,2],
  score = write_data[1:10,3], 
  moa_t1 <- c(
    "Non-anticancer drug",      # phentermine (拟交感胺/食欲抑制)
    "DNA adduct inducer",         # aflatoxin-b1 (DNA损伤/加合物形成)
    "MDM2/p53 inhibitor",         # carnosol (p53激活/干性抑制)
    "RAR agonist",                # tretinoin (维甲酸受体激动剂/Notch抑制)
    "PI3K/Akt inhibitor",         # tianeptine (PI3K通路抑制)
    "Non-anticancer drug",  # desoxycorticosterone (盐皮质激素受体激动剂)
    "Non-anticancer drug",      # nateglinide (钾通道阻滞剂)
    "IKK2 inhibitor",             # TPCA-1 (NF-κB信号抑制)
    "mGluR1 antagonist",          # JNJ-16259685 (代谢型谷氨酸受体拮抗)
    "Ferroptosis inducer"         # triptolide (铁死亡诱导/NF-κB抑制)
  )
)

colnames(drug_list_t1)[3]<-'action'
colnames(drug_list)[3]<-'action'

wholedrug<-rbind(drug_list_t1,drug_list)
whole_drug<-wholedrug[!duplicated(wholedrug$drug),]

whole_drug <- whole_drug %>% arrange(score)
whole_drug$drug <- factor(whole_drug$drug, levels = whole_drug$drug)


pdf('whole_Drug_function.pdf')
ggplot(whole_drug, aes(x = drug, y = score)) +
  geom_point(aes(fill = action), size = 5, shape = 21, color = "black", stroke = 1) +
  geom_text(aes(label = sprintf("%.2f", score)), 
            vjust = -1.2, size = 3.5, fontface = "bold") +
  scale_fill_brewer(palette = "Set3") +
  labs(y = "NES", 
       x = "Drug Candidate", 
       fill = "Mechanism of Action") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "right",
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  coord_cartesian(ylim =c(-2.42,-2.17) )

dev.off()
