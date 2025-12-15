# Supp Fig6a-c LNCAPenza correlation with CCLE  -----------------------------------

load("~/met_pc_cell_line/pcdata/GSE162225/RData/gene.expression.RData")
head(log2.tpm.matrix)
lncap.expr <- log2.tpm.matrix
colnames(lncap.expr) <- c(rep('enza_d0', 3), rep('enza_d14',3))


lncap.expr <- lncap.expr[,c(4,5,6)]
colnames(lncap.expr) <- c('LNCaP.enza.1','LNCaP.enza.2','LNCaP.enza.3')
gene.id <- intersect(rownames(CCLE.log2.tpm.matrix), rownames(lncap.expr))
gene.id <- intersect(gene.id, CCLE.rna.seq.marker.gene.1000)

lncap.ccle.cor <- cor(CCLE.log2.tpm.matrix[gene.id,], lncap.expr[gene.id,],method = 'spearman')
lncap.ccle.cor <- as.data.frame(lncap.ccle.cor)
lncap.enza.1 <- lncap.ccle.cor[,1,drop=F]
df <- lncap.enza.1

df$name      <- rownames(df)

df$prostate  <- ifelse(rownames(df) %in% pc.cell.line.name, 
                       "prostate", 
                       'others')

colnames(df)[1] <- "correlation"


df_ranked <- df %>% arrange(correlation)
df_ranked$rank <- 1:nrow(df_ranked)

ggplot(df_ranked, aes(x = rank, y = correlation, color = prostate)) +
  geom_point(size = 2.5) +
  
  # 高亮 prostate
  geom_point(
    data = subset(df_ranked, name %in% pc.cell.line.name),
    color = "red",
    size = 4
  ) +
  
  # 高亮 prostate标签
  geom_text_repel(
    data = subset(df_ranked, name %in% pc.cell.line.name),
    aes(label = gsub(x = name, pattern = '_PROSTATE','')),
    nudge_y = 0.03,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.size = 0.3,
    size = 4,
    color = "black"
  ) +
  
  scale_color_manual(values = c("prostate" = "red", "others" = "black")) +
  
  labs(
    x = "Rank",
    y = "Transcriptome similarity",
    title = "LNCaP-enza1"
  )+theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )

ggsave('./pcdata/paper/plot/section5/lncap.enza.1.pdf', width = 9, height = 8)


# lncap-enza2
lncap.enza.2 <- lncap.ccle.cor[,2,drop=F]
df <- lncap.enza.2

df$name      <- rownames(df)

df$prostate  <- ifelse(rownames(df) %in% pc.cell.line.name, 
                       "prostate", 
                       'others')

colnames(df)[1] <- "correlation"


df_ranked <- df %>% arrange(correlation)
df_ranked$rank <- 1:nrow(df_ranked)

ggplot(df_ranked, aes(x = rank, y = correlation, color = prostate)) +
  geom_point(size = 2.5) +
  
  # 高亮 prostate
  geom_point(
    data = subset(df_ranked, name %in% pc.cell.line.name),
    color = "red",
    size = 4
  ) +
  
  # 高亮 prostate标签
  geom_text_repel(
    data = subset(df_ranked, name %in% pc.cell.line.name),
    aes(label = gsub(x = name, pattern = '_PROSTATE','')),
    nudge_y = 0.03,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.size = 0.3,
    size = 4,
    color = "black"
  ) +
  
  scale_color_manual(values = c("prostate" = "red", "others" = "black")) +
  
  labs(
    x = "Rank",
    y = "Transcriptome similarity",
    title = "LNCaP-enza2"
  )+theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )

ggsave('./pcdata/paper/plot/section5/lncap.enza.2.pdf', width = 9, height = 8)


# lncap-enza3
lncap.enza.3 <- lncap.ccle.cor[,3,drop=F]
df <- lncap.enza.3

df$name      <- rownames(df)

df$prostate  <- ifelse(rownames(df) %in% pc.cell.line.name, 
                       "prostate", 
                       'others')

colnames(df)[1] <- "correlation"


df_ranked <- df %>% arrange(correlation)
df_ranked$rank <- 1:nrow(df_ranked)

ggplot(df_ranked, aes(x = rank, y = correlation, color = prostate)) +
  geom_point(size = 2.5) +
  
  # 高亮 prostate
  geom_point(
    data = subset(df_ranked, name %in% pc.cell.line.name),
    color = "red",
    size = 4
  ) +
  
  # 高亮 prostate标签
  geom_text_repel(
    data = subset(df_ranked, name %in% pc.cell.line.name),
    aes(label = gsub(x = name, pattern = '_PROSTATE','')),
    nudge_y = 0.03,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.size = 0.3,
    size = 4,
    color = "black"
  ) +
  
  scale_color_manual(values = c("prostate" = "red", "others" = "black")) +
  
  labs(
    x = "Rank",
    y = "Transcriptome similarity",
    title = "LNCaP-enza3"
  )+theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )

ggsave('./pcdata/paper/plot/section5/lncap.enza.3.pdf', width = 9, height = 8)





#Supp6 d-e LNCAP TP53/RB1 knock out ------------------------------------------------------
load(file = './pcdata/GSE175975_LNCAP_scrna/lncap.six.sample.ave.cpm.RData')
lncap.six.sample.ave.cpm 

lncap.ko  <- lncap.six.sample.ave.cpm[, 5:6]
lncap.ko.expr <- lncap.ko
gene.id <- intersect(rownames(lncap.ko.expr), rownames(CCLE.log2.tpm.matrix))
gene.id <- intersect(gene.id, CCLE.rna.seq.marker.gene.1000)

lncap.ko.cor <- cor(CCLE.log2.tpm.matrix[gene.id,], lncap.ko.expr[gene.id,], method = 'spearman')
lncap.ko.cor <- as.data.frame(lncap.ko.cor)
lncap.ko1    <- lncap.ko.cor[,2,drop=F]

##lncap.ko.1
df <- lncap.ko1
df$name      <- rownames(df)
df$prostate  <- ifelse(rownames(df) %in% pc.cell.line.name, 
                       "prostate", 
                       'others')

colnames(df)[1] <- "correlation"


df_ranked <- df %>% arrange(correlation)
df_ranked$rank <- 1:nrow(df_ranked)

ggplot(df_ranked, aes(x = rank, y = correlation, color = prostate)) +
  geom_point(size = 2.5) +
  
  # 高亮 prostate
  geom_point(
    data = subset(df_ranked, name %in% pc.cell.line.name),
    color = "red",
    size = 4
  ) +
  
  # 高亮 prostate标签
  geom_text_repel(
    data = subset(df_ranked, name %in% pc.cell.line.name),
    aes(label = gsub(x = name, pattern = '_PROSTATE','')),
    nudge_y = 0.03,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.size = 0.3,
    size = 4,
    color = "black"
  ) +
  
  scale_color_manual(values = c("prostate" = "red", "others" = "black")) +
  
  labs(
    x = "Rank",
    y = "Transcriptome similarity",
    title = "LNCaP-TP53/RB1-KO1"
  )+theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )

ggsave('./pcdata/paper/plot/section5/lncap.ko.1.pdf', width = 9, height = 8)

###lncap-ko2

lncap.ko2    <- lncap.ko.cor[,1,drop=F]
df <- lncap.ko2


df$name      <- rownames(df)
df$prostate  <- ifelse(rownames(df) %in% pc.cell.line.name, 
                       "prostate", 
                       'others')

colnames(df)[1] <- "correlation"


df_ranked <- df %>% arrange(correlation)
df_ranked$rank <- 1:nrow(df_ranked)

ggplot(df_ranked, aes(x = rank, y = correlation, color = prostate)) +
  geom_point(size = 2.5) +
  
  # 高亮 prostate
  geom_point(
    data = subset(df_ranked, name %in% pc.cell.line.name),
    color = "red",
    size = 4
  ) +
  
  # 高亮 prostate标签
  geom_text_repel(
    data = subset(df_ranked, name %in% pc.cell.line.name),
    aes(label = gsub(x = name, pattern = '_PROSTATE','')),
    nudge_y = 0.03,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.size = 0.3,
    size = 4,
    color = "black"
  ) +
  
  scale_color_manual(values = c("prostate" = "red", "others" = "black")) +
  
  labs(
    x = "Rank",
    y = "Transcriptome similarity",
    title = "LNCaP-TP53/RB1-KO2"
  )+theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )

ggsave('./pcdata/paper/plot/section5/lncap.ko.2.pdf', width = 9, height = 8)



#Supp6 f-g MET-oe ------------------------------------------------------------------

load("~/met_pc_cell_line/pcdata/vcap_cMET_high_bulk_data/gene.expression.RData")
cMET_high_expr <- log2.tpm.matrix
cMET_high_expr <- cMET_high_expr[,-6]
cMET_high_expr <- cMET_high_expr[,-6]
cMET_high_expr <- cMET_high_expr[, 5:6]
colnames(cMET_high_expr) <- c('MET_OE1', 'MET_OE2')

gene.id <- intersect(rownames(cMET_high_expr), rownames(CCLE.log2.tpm.matrix))
gene.id <- intersect(CCLE.rna.seq.marker.gene.1000, gene.id)
vcap.ccle.cor <- cor(CCLE.log2.tpm.matrix[gene.id,], cMET_high_expr[gene.id,], method = 'spearman')
vcap.ccle.cor <- as.data.frame(vcap.ccle.cor)
VCAP.OE1 <- vcap.ccle.cor[,2,drop=F]

df <- VCAP.OE1


df$name      <- rownames(df)
df$prostate  <- ifelse(rownames(df) %in% pc.cell.line.name, 
                       "prostate", 
                       'others')

colnames(df)[1] <- "correlation"


df_ranked <- df %>% arrange(correlation)
df_ranked$rank <- 1:nrow(df_ranked)

ggplot(df_ranked, aes(x = rank, y = correlation, color = prostate)) +
  geom_point(size = 2.5) +
  
  # 高亮 prostate
  geom_point(
    data = subset(df_ranked, name %in% pc.cell.line.name),
    color = "red",
    size = 4
  ) +
  
  # 高亮 prostate标签
  geom_text_repel(
    data = subset(df_ranked, name %in% pc.cell.line.name),
    aes(label = gsub(x = name, pattern = '_PROSTATE','')),
    nudge_y = 0.03,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.size = 0.3,
    size = 4,
    color = "black"
  ) +
  
  scale_color_manual(values = c("prostate" = "red", "others" = "black")) +
  
  labs(
    x = "Rank",
    y = "Transcriptome similarity",
    title = "VCaP-MET-OE1"
  )+theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )

ggsave('./pcdata/paper/plot/section5/vcap.oe.1.pdf', width = 9, height = 8)


####vcap oe 2

VCAP.OE2 <- vcap.ccle.cor[,1,drop=F]

df <- VCAP.OE2


df$name      <- rownames(df)
df$prostate  <- ifelse(rownames(df) %in% pc.cell.line.name, 
                       "prostate", 
                       'others')

colnames(df)[1] <- "correlation"


df_ranked <- df %>% arrange(correlation)
df_ranked$rank <- 1:nrow(df_ranked)

ggplot(df_ranked, aes(x = rank, y = correlation, color = prostate)) +
  geom_point(size = 2.5) +
  
  # 高亮 prostate
  geom_point(
    data = subset(df_ranked, name %in% pc.cell.line.name),
    color = "red",
    size = 4
  ) +
  
  # 高亮 prostate标签
  geom_text_repel(
    data = subset(df_ranked, name %in% pc.cell.line.name),
    aes(label = gsub(x = name, pattern = '_PROSTATE','')),
    nudge_y = 0.03,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.size = 0.3,
    size = 4,
    color = "black"
  ) +
  
  scale_color_manual(values = c("prostate" = "red", "others" = "black")) +
  
  labs(
    x = "Rank",
    y = "Transcriptome similarity",
    title = "VCaP-MET-OE2"
  )+theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )
ggsave('./pcdata/paper/plot/section5/vcap.oe.2.pdf', width = 9, height = 8)


#Fig.6h atac organoid -----------------------------------------------------------

load(file = '/home/liuxueying/met_pc_cell_line/pcdata/atac/data/three.type.meta.data.RData')##stem.id
load(file = '/home/liuxueying/met_pc_cell_line/pcdata/paper/data/section55/atac.organoid.RData')


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
# atac organoid and PC3 sample 

organoid.sample.name     <- organoid.metadata$Sample_Name
atac.organoid            <- atac.organoid[,organoid.sample.name]###这样提取不对，永远都是第一个，重复值
colnames(atac.organoid)

colnames(atac.organoid)  <- organoid.metadata$Sample_Name_3


# 颜色对应

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







