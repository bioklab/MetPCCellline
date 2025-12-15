# Fig.5a PCA ATAC-seq -----------------------------------------------------


load(file = '/home/liuxueying/met_pc_cell_line/pcdata/paper/data/section55/atac.count.RData')


pca.rs                       <- prcomp(atac.count %>% t)    
pca                          <- pca.rs$x[,1:2]
pca                          <- as.data.frame(pca)
pca$cell.line                <- gsub(x = rownames(pca), pattern = '_PROSTATE',replacement = '')

ggplot(pca,aes(x= PC1, y = PC2, color = cell.line))+
  geom_point(size = 8)+
  geom_text_repel(data = subset(pca), 
                  aes(label = gsub("VCAP", "VCaP", 
                                   gsub("NCIH660", "NCI-H660", 
                                        gsub("MDAPCA2B", "MDA-PCa-2b", 
                                             gsub("LNCAPCLONEFGC", "LNCaP", cell.line))))),
                  nudge_y = 0.02,  
                  color = "black",
                  size = 6,  # Increase the font size
                  # Set the angle of the text (e.g., 45 degrees)
                  arrow = arrow(length = unit(0.02, "npc")))+
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )+
  theme(legend.position = 'none')
ggsave('/home/liuxueying/met_pc_cell_line/pcdata/paper/plot/section4/pca.pdf', width = 8, height = 8)


#Fig.5b Top 1000 peaks enrichment MSigDB Hallmark -------------------------------
### FUNCTION: getGRanges() ---
getGRanges <- function(dat){
  dat <- subset(dat, select=c("seqnames","start","end"))
  gr <- GenomicRanges::makeGRangesFromDataFrame(dat, keep.extra.columns=FALSE)
  return(gr)
}


### FUNCTION: processEnrichment() ---
processEnrichment <- function(dat, pval){
  dat <- subset(dat, select=c("id","fold_enrichment","p_value","p_adjust"))
  dat <- dat[which(dat$p_adjust <= pval),]
  dat <- dat[order(dat$p_adjust, decreasing=FALSE),]
  dat$nlogp <- -log10(dat$p_adjust)
  
  return(dat)
}

### FUNCTION: plotBar() ----
plotBar <- function(dat, analysis_id){
  # PREPARE DATA ---
  dat <- dat[order(dat$nlogp, decreasing=FALSE),]
  dat$id <- factor(dat$id, levels=dat$id)
  
  # PLOT ---
  p <- ggplot(dat, aes(x=nlogp, y=id)) +
    geom_bar(stat="identity", fill="#006094", width=0.8) +
    coord_cartesian(xlim=c(0,12)) +
    scale_x_continuous(breaks=seq(0,12,by=4)) +
    theme(
      axis.text.x = element_text(size = 11, color="#000000"),
      axis.text.y = element_text(size = 11, color="#000000", ),
      axis.title = element_text(size = 11, color="#000000"),
      plot.title = element_text(size = 11, color="#000000", hjust=0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_line(size=0.2, color="#000000"), 
      strip.text = element_text(size=11, color="#000000"),
      strip.background = element_rect(fill="#FFFFFF", color="#FFFFFF"),
      panel.background = element_rect(fill="#FFFFFF", color="#000000"),
      legend.text = element_text(size = 11, color="#000000"),
      legend.title = element_blank(),
      legend.key.size = unit(0.3, "cm"),
      legend.position = "none") +
    ylab("") +
    xlab("-log10(FDR)") + 
    ggtitle(analysis_id) 
  
  return(p)
}


load(file = '/home/liuxueying/met_pc_cell_line/pcdata/paper/data/section55/feature_annotation.RData')

gr_all             <- getGRanges(dat = feature_annotation)

### SET SEED ---
set.seed(12345)

### GREAT ENRICHMENT ---

obj_all        <- rGREAT::great(gr=gr_all, gene_sets="MSigDB:H", tss_source="txdb:hg38", 
                                biomart_dataset = NULL, min_gene_set_size = 5, 
                                mode = "basalPlusExt", basal_upstream = 5000, basal_downstream = 1000, extension = 1000000, 
                                extended_tss = NULL, background = NULL, exclude = "gap", cores = 50, verbose = TRUE)
### PARSE ENRICHMENT DATA ---

dat_all        <- obj_all@table

df_all         <- processEnrichment(dat=dat_all, pval=0.05)

p4             <- plotBar(dat=df_all, analysis_id="MSigDB Hallmark")

pdf('/home/liuxueying/met_pc_cell_line/pcdata/paper/plot/section4/enrich.pdf',width = 8 ,height = 6)
p4
dev.off()




# ATAC organoid correlation with MSPC samples -----------------------------


targets <- c(
  "MSK-PCa15_1",
  "MSK-PCa12_1",
  "MSK-PCa8_1",
  "MSK-PCa9_1",
  "MSK-PCa11_1",
  "MSK-PCa3_1",
  "MSK-PCa17_1",
  "MSK-PCa20_1",
  "MSK-PCa13_1",
  "PC3_1",
  "MSK-PCa18_1"
  
)
#save(list = c('metadata_tang2022_wcdt','atac.three.subtype.mcrpc', 'atac.organoid'), file = '/home/liuxueying/met_pc_cell_line/pcdata/paper/data/section55/atac.organoid.RData')

organoid.metadata             <- filter(metadata_tang2022_wcdt, Dataset == 'TANG2022')
organoid.metadata$Sample_Name <- make.unique(organoid.metadata$Sample_Name)
atac.organoid                 <- atacseq_read_counts_tang2022_wcdt_combined
colnames(atac.organoid)       <- make.unique(colnames(atac.organoid))
organoid.metadata             <- filter(organoid.metadata,  Sample_Name_2 %in% targets)
organoid.metadata$Sample_Name_3 <- sub("_.*", "", organoid.metadata$Sample_Name_2)
# atac organoid and PC3 sample ---------------------------------------------------

organoid.sample.name     <- organoid.metadata$Sample_Name
atac.organoid            <- atac.organoid[,organoid.sample.name]###这样提取不对，永远都是第一个，重复值
colnames(atac.organoid)

colnames(atac.organoid)  <- organoid.metadata$Sample_Name_3


# 颜色对应 --------------------------------------------------------------------

# organoid_levels <- levels(cor.stem.res$)
# color_map <- setNames(
#   ifelse(organoid_levels == "PC3_PROSTATE", "firebrick", 'steelblue'),
#   organoid_levels
# )

my.cols <- c(
  "MSK-PCa15"   = 'steelblue',
  "MSK-PCa12"   = 'steelblue',
  "MSK-PCa8"    = 'steelblue',
  "MSK-PCa9"    = 'steelblue',
  "PC3"         = "firebrick",
  "MSK-PCa11"   = 'steelblue',
  'MSK-PCa3'    = 'steelblue',
  "MSK-PCa17"   = 'steelblue',
  "MSK-PCa20"   = 'steelblue',
  "MSK-PCa13"   = 'steelblue',
  "MSK-PCa18"   = 'steelblue'
  
)


cor.stem.res             <- cor(atac.three.subtype.mcrpc[rownames(mat_features),stem.id], atac.organoid[rownames(mat_features),], method = 'spearman')
cor.stem.res %>% reshape2::melt() %>% 
  ggplot(aes(x = reorder(Var2,value, median), y = value , fill = Var2, color =Var2))+
  geom_boxplot(alpha = 0.5)+
  scale_fill_manual(values = my.cols)+
  scale_color_manual(values = my.cols)+
  geom_jitter(
    aes(color = Var2),                # 取消 Var2 默认颜色映射
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
    size = 1.2,
    # 点的颜色
    # 点透明度
  )+
  
  
  stat_compare_means(comparisons = list(c('PC3','MSK-PCa12')
  ))+xlab('')+ylab('Correlation')+theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )+RotatedAxis()

ggsave('./pcdata/paper/plot/section5/organoid.atac.pdf', width = 12, height = 6)


# ATAC subtype samples correlation with CCLE  Fig.5c-e-----------------------------


## good vcap  high
mat_features <- readRDS("~/met_pc_cell_line/pcdata/atac/data/tang2022_atacseq_dqnorm_combat_top1percent.rds")
load(file = '/home/liuxueying/met_pc_cell_line/pcdata/atac/data/three.type.meta.data.RData')##AR.id, neuroendocrine.id, stem.id
load(file = '/home/liuxueying/met_pc_cell_line/pcdata/paper/data/section55/atac.count.RData')##cell line atac 
load(file = '/home/liuxueying/met_pc_cell_line/pcdata/paper/data/section1/atac.three.subtype.mcrpc.RData')#ATAC 37 samples ATAC count


cor.AR.res                <- cor(atac.three.subtype.mcrpc[rownames(mat_features),AR.id], atac.count[rownames(mat_features),], method = 'spearman')
boxplot(cor.AR.res)
cor.AR.res %>% reshape2::melt() %>% 
  ggplot(aes(x = reorder(Var2,value, median), y = value, fill = Var2, colour = Var2))+
  
  geom_boxplot(alpha = 0.5)+
  scale_fill_manual(values = my.cols)+
  scale_colour_manual(values = my.cols)+
  geom_jitter(
    aes(color = Var2),                # 取消 Var2 默认颜色映射
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
    size = 1.2,
    # 点的颜色
    # 点透明度
  )+
  stat_compare_means(comparisons = list(
    c('VCaP','LNCaP')), method = 'wilcox.test', paired = T)+
  xlab('')+
  ylab('')+
  # scale_x_discrete(labels = c('NCI-H660' = paste('NCI-H660', sep = '\n', '(n = 26)'),
  #                             '22RV1' = paste('22RV1', sep = '\n', '(n = 26)'),
  #                             'DU145' = paste('DU145', sep = '\n', '(n = 26)'),
  #                             'LNCaP'=paste('LNCaP', sep = '\n','(n= 26)'),
  #                           
  #                             'PC3'=paste('PC3',sep = '\n','(n = 26)'),
  #                             'VCaP'=paste('VCaP',sep = '\n','(n = 26)')
  # ))+
  ylab('Correlation')+
  theme_bw() +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 35),
    axis.title = element_text(size = 35, hjust = 0.5, color = 'black'),
    axis.text = element_text(size = 35, color = 'black'),
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA),
    axis.ticks.length = unit(0.6, "cm")
    #axis.ticks = element_line(size = 1.8)
  )+
  theme(legend.position = "none")+
  ggtitle('ARPC')
ggsave('/home/liuxueying/met_pc_cell_line/pcdata/paper/plot/section4/ARPC.atac.pdf', width = 15, height = 8)
# ATAC NEPC samples correlation with cell line ----------------------------
##good NCIH660 high
cor.NE.res                <- cor(atac.three.subtype.mcrpc[rownames(mat_features),neuroendocrine.id], atac.count[rownames(mat_features),], method = 'spearman')
#boxplot(cor.NE.res)

cor.NE.res %>% reshape2::melt() %>% 
  ggplot(aes(x = reorder(Var2,value, median), y = value, fill = Var2, colour = Var2))+
  geom_boxplot(alpha = 0.5)+
  scale_fill_manual(values = my.cols)+
  scale_colour_manual(values = my.cols)+
  geom_jitter(
    aes(color = Var2),                # 取消 Var2 默认颜色映射
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
    size = 1.2,
    # 点的颜色
    # 点透明度
  )+
  stat_compare_means(comparisons = list(c('NCI-H660','DU145'),c('NCI-H660','PC3')), method = 'wilcox.test',paired = T)+
  scale_x_discrete(labels = c('NCI-H660' = 'NCI-H660',
                              '22RV1' = '22RV1',
                              'DU145' = 'DU145',
                              'LNCaP'='LNCaP',
                              
                              'PC3'='PC3',
                              'VCaP'='VCaP')
  )+
  # scale_x_discrete(labels = c('NCI-H660' = paste('NCI-H660', sep = '\n', '(n = 4)'),
  #                             '22RV1' = paste('22RV1', sep = '\n', '(n = 4)'),
  #                             'DU145' = paste('DU145', sep = '\n', '(n = 4)'),
  #                             'LNCaP'=paste('LNCaP', sep = '\n','(n= 4)'),
  #                             
  #                             'PC3'=paste('PC3',sep = '\n','(n = 4)'),
  #                             'VCaP'=paste('VCaP',sep = '\n','(n = 4)')
  # ))+
  ylab('Correlation')+
  theme_bw() +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 35),
    axis.title = element_text(size = 35, hjust = 0.5, color = 'black'),
    axis.text = element_text(size = 35, color = 'black'),
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA),
    axis.ticks.length = unit(0.6, "cm")
    #axis.ticks = element_line(size = 1.8)
  )+
  theme(legend.position = "none")+
  xlab('')+
  
  ggtitle('NEPC')
ggsave('/home/liuxueying/met_pc_cell_line/pcdata/paper/plot/section4/NEPC.atac.pdf', width = 15, height = 8)


# ATAC MSPC samples correlation with cell line ----------------------------
## stem DU145 PC3 NCIH660 high

cor.stem.res              <- cor(atac.three.subtype.mcrpc[rownames(mat_features),stem.id], atac.count[rownames(mat_features),], method = 'spearman')
#boxplot(cor.stem.res)
cor.stem.res %>% reshape2::melt() %>% 
  ggplot(aes(x = reorder(Var2,value, median), y = value, fill = Var2, colour = Var2))+
  geom_boxplot(alpha = 0.5)+
  scale_fill_manual(values = my.cols)+
  scale_colour_manual(values = my.cols)+
  geom_jitter(
    aes(color = Var2),                # 取消 Var2 默认颜色映射
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
    size = 1.2,
    # 点的颜色
    # 点透明度
  )+
  stat_compare_means(comparisons = list(c('NCI-H660','DU145'),c('NCI-H660','PC3')), method = 'wilcox.test',paired = T)+
  # scale_x_discrete(labels = c('NCI-H660' = 'NCI-H660',
  #                             '22RV1' = '22RV1',
  #                             'DU145' = 'DU145',
  #                             'LNCaP'='LNCaP',
  #                             
  #                             'PC3'='PC3',
  #                             'VCaP'='VCaP')
  # )+
  ylab('Correlation')+
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )+
  theme(legend.position = "none")+
  
  xlab('')+
  ggtitle('MSPC')
ggsave('/home/liuxueying/met_pc_cell_line/pcdata/paper/plot/section4/MSPC.atac.pdf', width = 15, height = 8)

# ATAC ARPC samples correlation with cell line -------------------------------------------------------------


## good vcap  high
cor.AR.res                <- cor(atac.three.subtype.mcrpc[rownames(mat_features),AR.id], atac.count[rownames(mat_features),], method = 'spearman')
boxplot(cor.AR.res)
cor.AR.res %>% reshape2::melt() %>% 
  ggplot(aes(x = reorder(Var2,value, median), y = value, fill = Var2, colour = Var2))+
  geom_boxplot(alpha = 0.5)+
  scale_fill_manual(values = my.cols)+
  scale_colour_manual(values = my.cols)+
  geom_jitter(
    aes(color = Var2),                # 取消 Var2 默认颜色映射
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
    size = 1.2,
    # 点的颜色
    # 点透明度
  )+
  stat_compare_means(comparisons = list(c('NCI-H660','DU145'),c('NCI-H660','PC3')), method = 'wilcox.test',paired = T)+
  scale_x_discrete(labels = c('NCI-H660' = 'NCI-H660',
                              '22RV1' = '22RV1',
                              'DU145' = 'DU145',
                              'LNCaP'='LNCaP',
                              
                              'PC3'='PC3',
                              'VCaP'='VCaP')
  )+
  ylab('Correlation')+
  ggplot.style+
  xlab('')+
  ylab('Correlation')+
  ggtitle('ARPC')

ggsave('/home/liuxueying/met_pc_cell_line/pcdata/paper/plot/section55/MSPC.atac.pdf', width = 8, height = 4)




