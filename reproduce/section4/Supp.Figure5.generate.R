
# Supp Fig5a --------------------------------------------------------------

load("~/met_pc_cell_line/pcdata/paper/data/section4/mcrpc.three.sample.seob.RData")

FeaturePlot(mcrpc.three.sample.seob, features = c('AR'))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('AR')+
  NoLegend()

ggsave('./pcdata/paper/plot/section4/ar.mcrpc.pdf', width = 5, height = 5)

FeaturePlot(mcrpc.three.sample.seob, features = c('SYP'))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('SYP')+
  NoLegend()

ggsave('./pcdata/paper/plot/section4/syp.crpc.pdf', width = 5, height = 5)

FeaturePlot(mcrpc.three.sample.seob, features = c('CD44'))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('CD44')+
  NoLegend()

ggsave('./pcdata/paper/plot/section4/cd44.crpc.pdf', width = 5, height = 5)

DimPlot(mcrpc.three.sample.seob, group.by = 'three.type',label = T, 
        cols = c('MSPC'='#D65813', 'ARPC' = '#2CAD3F', 'NEPC'='#6792CD'))+
  ggtitle('')+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('Subtype')+
  NoLegend()

ggsave('./pcdata/paper/plot/section4/mcrpc.subtype.pdf', width = 5, height = 5) 

pdf('./pcdata/paper/plot/section4/mcrpc.pdf',width = 20, height = 5)
p4|p1|p2|p3
dev.off()



#Supp Fig5b metastatic MSPC correlation with ccle  -----------------------------------


stem.cell.pesudobulk <- AggregateExpression(
  mcrpc.three.sample.seob,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "three.type",
  add.ident = NULL,
  slot = "counts",
  verbose = TRUE
)

stem.cell.pesudobulk                     <- stem.cell.pesudobulk$RNA
stem.cell.pesudobulk.cpm                 <- apply(stem.cell.pesudobulk, 2, function(x) { x/sum(x)*1000000 })
stem.cell.pesudobulk.log2cpm             <- log2(stem.cell.pesudobulk.cpm +1)

stem.cell.pesudobulk.log2cpm             <- change_rownames_symbol_to_ensemble(stem.cell.pesudobulk.log2cpm)




marker.gene           <- intersect(rownames(stem.cell.pesudobulk.log2cpm),CCLE.rna.seq.marker.gene.1000)  
marker.gene           <- intersect(rownames(CCLE.log2.tpm.matrix),marker.gene) 
ccle.stem.cor         <- cor(stem.cell.pesudobulk.log2cpm[marker.gene,],CCLE.log2.tpm.matrix[marker.gene,],method='spearman')
ccle.stem.cor         <- ccle.stem.cor %>% t %>% as.data.frame()
ccle.stem.cor[pc.cell.line.name,]###MDA-PCA-2b PC3

df           <- ccle.stem.cor[, 'MSPC', drop = FALSE]
df$name      <- rownames(df)
df$prostate  <- ifelse(rownames(df) %in% pc.cell.line.name, 
                       "prostate", 
                       "others")
colnames(df)[1] <- "correlation"

df$rank     <- rank(df$correlation)
df_ranked   <- df[order(df$rank), ]

#write.xlsx(df_ranked, file = "mspc.tc.xlsx")


#--- B. plot ---
highlight_points <- df_ranked[df_ranked$rank %in% c(1018, 1019,1017), ]
ggplot(df_ranked, aes(x = rank, y = correlation, color = prostate)) +
  geom_point(size = 3) +  
  labs(
    # 标题使用该列名，或你想替换的名字
    x = "Rank",
    y = "Transcriptome similarity"
  ) +
  scale_color_manual(values = c("prostate" = "red", "others" = "black")) +
  # 在图中高亮 pc.cell.line.name 这几行
  geom_point(
    data = subset(df, name %in% pc.cell.line.name),
    color = "red",
    size = 4
  )+
  geom_text_repel(data = subset(df, name %in% pc.cell.line.name), 
                  aes(label = gsub(x= name, pattern = '_PROSTATE','')), 
                  nudge_y = 0.02,  # 标签在纵坐标上的偏移量
                  color = "black",
                  size = 4,
                  arrow = arrow(length = unit(0.02, "npc")))+
  geom_text_repel(data = highlight_points ,
                  aes(label = name),
                  size = 4,
                  
                  color = "blue",
                  nudge_y = 0.02)+
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )+ggtitle('MSPC')###

ggsave('./pcdata/paper/plot/section4/ARPC.supp.MCRPC.pdf', width = 8, height = 8)
#Supp Fig5cmetastatic ARPC correlation with CCLE--------------------------------------------------------------------

df           <- ccle.stem.cor[, 'ARPC', drop = FALSE]
df$name      <- rownames(df)
df$prostate  <- ifelse(rownames(df) %in% pc.cell.line.name, 
                       "prostate", 
                       "others")
colnames(df)[1] <- "correlation"

df$rank     <- rank(df$correlation)
df_ranked   <- df[order(df$rank), ]

#write.xlsx(df_ranked, file = "mspc.tc.xlsx")


#--- B. plot ---
highlight_points <- df_ranked[df_ranked$rank %in% c(1018, 1019,1017), ]
ggplot(df_ranked, aes(x = rank, y = correlation, color = prostate)) +
  geom_point(size = 3) +  
  labs(
    # 标题使用该列名，或你想替换的名字
    x = "Rank",
    y = "Transcriptome similarity"
  ) +
  scale_color_manual(values = c("prostate" = "red", "others" = "black")) +
  # 在图中高亮 pc.cell.line.name 这几行
  geom_point(
    data = subset(df, name %in% pc.cell.line.name),
    color = "red",
    size = 4
  )+
  geom_text_repel(data = subset(df, name %in% pc.cell.line.name), 
                  aes(label = gsub(x= name, pattern = '_PROSTATE','')), 
                  nudge_y = 0.02,  # 标签在纵坐标上的偏移量
                  color = "black",
                  size = 4,
                  arrow = arrow(length = unit(0.02, "npc")))+
  
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )+ggtitle('ARPC')###
ggsave('./pcdata/paper/plot/section4/ARPC.supp.mcrpc.pdf', width = 8, height = 8)

#Supp Fig.5dmetastatic NEPC correlation with CCLE--------------------------------------------------------------------


df           <- ccle.stem.cor[, 'NEPC', drop = FALSE]
df$name      <- rownames(df)
df$prostate  <- ifelse(rownames(df) %in% pc.cell.line.name, 
                       "prostate", 
                       "others")
colnames(df)[1] <- "correlation"

df$rank     <- rank(df$correlation)
df_ranked   <- df[order(df$rank), ]

#write.xlsx(df_ranked, file = "mspc.tc.xlsx")


#--- B. plot ---
highlight_points <- df_ranked[df_ranked$rank %in% c(1018, 1019,1017), ]
ggplot(df_ranked, aes(x = rank, y = correlation, color = prostate)) +
  geom_point(size = 3) +  
  labs(
    # 标题使用该列名，或你想替换的名字
    x = "Rank",
    y = "Transcriptome similarity"
  ) +
  scale_color_manual(values = c("prostate" = "red", "others" = "black")) +
  # 在图中高亮 pc.cell.line.name 这几行
  geom_point(
    data = subset(df, name %in% pc.cell.line.name),
    color = "red",
    size = 4
  )+
  geom_text_repel(data = subset(df, name %in% pc.cell.line.name), 
                  aes(label = gsub(x= name, pattern = '_PROSTATE','')), 
                  nudge_y = 0.02,  # 标签在纵坐标上的偏移量
                  color = "black",
                  size = 4,
                  arrow = arrow(length = unit(0.02, "npc")))+
  
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )+ggtitle('NEPC')###
ggsave('./pcdata/paper/plot/section4/NEPC.supp.mcrpc.pdf', width = 8, height = 8)




#Supp Fig.4g check stem samples -------------------------------------------------------
load( file = '/home/liuxueying/met_pc_cell_line/pcdata/paper/data/section4/ATAC.sample.metadata.RData')
load(file = '/home/liuxueying/met_pc_cell_line/pcdata/atac/data/three.type.meta.data.RData')
load(file = '/home/liuxueying/met_pc_cell_line/pcdata/paper/data/section55/ssGSEA.score.df.RData')
ssGSEA.score.df$Sample_Name <- rownames(ssGSEA.score.df)
wcdt.metadata <- left_join(wcdt.metadata, ssGSEA.score.df, by = 'Sample_Name')

my.cols <- c(
  "ARPC" = '#2CAD3F',
  "MSPC"    = '#D65813',
  "NEPC"    = '#6792CD'
  
)


wcdt.metadata %>% filter(three.type != 'others') %>% 
  mutate(
    three.type = factor(
      three.type,
      levels = c("ARPC", "NEPC", "MSPC")   # ⭐ 指定顺序
    )
  ) %>% 
  ggplot(aes(x = three.type, y = EMT.marker.gene, fill = three.type))+
  geom_boxplot()+
  scale_fill_manual(values = my.cols)+
  geom_jitter(
    aes(color = NULL),                # 取消 Var2 默认颜色映射
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
    size = 1.2,
    color = "black",                  # 点的颜色
    # 点透明度
  )+
  stat_compare_means(comparisons = list( c('MSPC','ARPC')), paired = F)+
  xlab('')+
  ylab('EMT score')+
  scale_x_discrete(labels = c('ARPC' = paste('ARPC', sep = '\n', '(n = 26)'),
                              'MSPC' = paste('MSPC', sep = '\n', '(n = 7)'),
                              'NEPC' = paste('NEPC', sep = '\n', '(n = 4)')
  ))+
  theme_bw() +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 35),
    axis.title = element_text(size = 35, hjust = 0.5, color = 'black'),
    axis.text = element_text(size = 35, color = 'black'),
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA),
    axis.ticks.length = unit(0.4, "cm")
    #axis.ticks = element_line(size = 1.8)
  )+
  theme(legend.position = "none")
ggsave('/home/liuxueying/met_pc_cell_line/pcdata/paper/plot/section55/check.MSPC.patients.pdf', width = 8, height = 8)

#Supp Fig.4f check neuroendocrine samples --------------------------------------------
wcdt.metadata %>% filter(three.type != 'others') %>% 
  mutate(
    three.type = factor(
      three.type,
      levels = c("ARPC", "NEPC", "MSPC")   # ⭐ 指定顺序
    )
  ) %>% 
  ggplot(aes(x = three.type, y = neuroendocrine.marker.gene,fill = three.type))+
  geom_boxplot()+
  scale_fill_manual(values = my.cols)+
  geom_jitter(
    aes(color = NULL),                # 取消 Var2 默认颜色映射
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
    size = 1.2,
    color = "black",                  # 点的颜色
    # 点透明度
  )+
  stat_compare_means(comparisons = list(c("MSPC", "NEPC")), paired = F)+
  scale_x_discrete(labels = c('ARPC' = paste('ARPC', sep = '\n', '(n = 26)'),
                              'MSPC' = paste('MSPC', sep = '\n', '(n = 7)'),
                              'NEPC' = paste('NEPC', sep = '\n', '(n = 4)')
  ))+
  xlab('')+
  ylab('NE score')+
  theme_bw() +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 35),
    axis.title = element_text(size = 35, hjust = 0.5, color = 'black'),
    axis.text = element_text(size = 35, color = 'black'),
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA),
    axis.ticks.length = unit(0.4, "cm")
    #axis.ticks = element_line(size = 1.8)
  )+
  theme(legend.position = "none")
ggsave('/home/liuxueying/met_pc_cell_line/pcdata/paper/plot/section55/check.NEPC.patients.pdf', width = 8, height = 8)


#Supp Fig4e check AR samples --------------------------------------------------------
wcdt.metadata %>% filter(three.type != 'others') %>% 
  mutate(
    three.type = factor(
      three.type,
      levels = c("ARPC", "NEPC", "MSPC")   # ⭐ 指定顺序
    )
  ) %>% 
  ggplot(aes(x = three.type, y = AR.marker.gene, fill = three.type))+
  geom_boxplot()+
  scale_fill_manual(values = my.cols)+
  geom_jitter(
    aes(color = NULL),                # 取消 Var2 默认颜色映射
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
    size = 1.2,
    color = "black",                  # 点的颜色
    # 点透明度
  )+
  stat_compare_means(comparisons = list( c('MSPC','ARPC')), paired = F)+
  xlab('')+
  ylab('AR score')+
  scale_x_discrete(labels = c('ARPC' = paste('ARPC', sep = '\n', '(n = 26)'),
                              'MSPC' = paste('MSPC', sep = '\n', '(n = 7)'),
                              'NEPC' = paste('NEPC', sep = '\n', '(n = 4)')))+
  theme_bw() +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 35),
    axis.title = element_text(size = 35, hjust = 0.5, color = 'black'),
    axis.text = element_text(size = 35, color = 'black'),
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA),
    axis.ticks.length = unit(0.4, "cm")
    #axis.ticks = element_line(size = 1.8)
  )+
  theme(legend.position = "none")
ggsave('/home/liuxueying/met_pc_cell_line/pcdata/paper/plot/section55/check.ARPC.patients.pdf', width = 8, height = 8)





