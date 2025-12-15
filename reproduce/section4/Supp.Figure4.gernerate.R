# Supp Fig4a-c AR CD44 SYP expr external bulk pc3 data------------------------------------------------------------
###good result
### bulk
load("~/met_pc_cell_line/pcdata/GSE116668/RData/gene.expression.RData")
log2.tpm.matrix %>% colnames()
head(log2.tpm.matrix)


colnames(log2.tpm.matrix) <- c('pc3.1','pc3.2','pc3.3', rep('treat', 3))

log2.pc3.matrix <- log2.tpm.matrix[, 1:3]

gene.id <- intersect(rownames(CCLE.log2.tpm.matrix), rownames(log2.pc3.matrix))
ccle.tpm <- CCLE.log2.tpm.matrix[gene.id,]
pc3 <- log2.pc3.matrix[gene.id,]
identical(rownames(ccle.tpm), rownames(pc3))
ccle.pc3         <- cbind(pc3, ccle.tpm)


id                                <- c(pc.cell.line.name, 'pc3.1','pc3.2','pc3.3')
ccle.pc3.data                     <- 2^ccle.pc3-1
ccle.pc3.data.cpm                 <- apply(ccle.pc3.data, 2, function(x) { x/sum(x)*1000000 })
ccle.pc3.data.log2cpm             <- log2(ccle.pc3.data.cpm +1)
ccle.pc3.data.log2cpm.symbol <- change_rownames_ensemble_to_symbol(ccle.pc3.data.log2cpm)
ccle.pc3.data.log2cpm[symbol.to.ensemble('CD44'), c(pc.cell.line.name,colnames(log2.pc3.matrix))] %>% sort()
ccle.pc3.data.log2cpm[symbol.to.ensemble('SYP'), c(pc.cell.line.name,colnames(log2.pc3.matrix))] %>% sort()
ccle.pc3.data.log2cpm[symbol.to.ensemble('AR'), c(pc.cell.line.name,colnames(log2.pc3.matrix))] %>% sort()

ccle.pc3.data.log2cpm.symbol <- change_rownames_ensemble_to_symbol(ccle.pc3.data.log2cpm)
expr <- ccle.pc3.data.log2cpm.symbol[c('CD44','AR','SYP'), id] %>% t %>% as.data.frame()
expr$cell.line <- gsub(x = rownames(expr), pattern = '_PROSTATE', replacement = '')



ggplot(expr, aes(x = reorder(cell.line, CD44, decreasing = T), y = CD44))+
  geom_col(width = 0.6, fill = "#D8CDEB")+
  scale_x_discrete(labels = c('NCIH660'='NCI-H660',
                              'VCAP'='VCaP',
                              'PC3' = 'PC3 (CCLE)',
                              'MDAPCA2B'='MDA-PCa-2b',
                              'LNCAPCLONEFGC'='LNCaP',
                              'pc3.1'='PC3',
                              'pc3.2' = 'PC3',
                              'pc3.3' = 'PC3'))+
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
  RotatedAxis()+
  ylab('CD44 expression')+xlab('')

ggsave('./pcdata/paper/plot/section4/CD44.expr.pdf', width = 12, height = 12)

ggplot(expr, aes(x = reorder(cell.line, AR, decreasing = T), y = AR))+
  geom_col(width = 0.6, fill = "#D8CDEB")+
  scale_x_discrete(labels = c('NCIH660'='NCI-H660',
                              'VCAP'='VCaP',
                              'PC3' = 'PC3 (CCLE)',
                              'MDAPCA2B'='MDA-PCa-2b',
                              'LNCAPCLONEFGC'='LNCaP',
                              'pc3.1'='PC3',
                              'pc3.2' = 'PC3',
                              'pc3.3' = 'PC3'))+
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
  RotatedAxis()+
  ylab('AR expression')+xlab('')

ggsave('./pcdata/paper/plot/section4/AR.expr.pdf', width = 12, height = 12)

ggplot(expr, aes(x = reorder(cell.line, SYP, decreasing = T), y = SYP))+
  geom_col(width = 0.6, fill = "#D8CDEB")+
  scale_x_discrete(labels = c('NCIH660'='NCI-H660',
                              'VCAP'='VCaP',
                              'PC3' = 'PC3 (CCLE)',
                              'MDAPCA2B'='MDA-PCa-2b',
                              'LNCAPCLONEFGC'='LNCaP',
                              'pc3.1'='PC3',
                              'pc3.2' = 'PC3',
                              'pc3.3' = 'PC3'))+
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
  RotatedAxis()+
  ylab('SYP expression')+xlab('')



# supp 4d-e expression in scrna-pc3 -----------------------------------------

pc3.seob                    <- readRDS("~/met_pc_cell_line/pcdata/met_pc/pc3.seob.rds")
p4 <- FeaturePlot(pc3.seob, features = c('CD44'), pt.size = 2.5)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        plot.title = element_text(hjust = 0.5, size = 16))+
  ggtitle('CD44')+
  NoLegend()
p5 <-  FeaturePlot(pc3.seob, features = c('AR'), pt.size = 2.5)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        plot.title = element_text(hjust = 0.5, size = 16))+
  ggtitle('AR')+
  NoLegend()
p6 <-  FeaturePlot(pc3.seob, features = c('SYP'), pt.size = 2.5)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        plot.title = element_text(hjust = 0.5, size = 16))+
  ggtitle('SYP')+
  NoLegend()
pdf(file = './pcdata/paper/plot/section4/CD44.AR.SYP.supp.scRNApc3.pdf', width = 15, height = 6)
p5|p6|p4
dev.off()






load("~/met_pc_cell_line/pcdata/GSE140440_smart-seq_pc3/GSE140440.seob.RData")
FeaturePlot(GSE140440.seob, features = c('CD44','SYP','AR','TP63'))
p4 <- FeaturePlot(GSE140440.seob, features = c('CD44'), pt.size = 2.5)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        plot.title = element_text(hjust = 0.5, size = 16))+
  ggtitle('CD44')+
  NoLegend()
p5 <-  FeaturePlot(GSE140440.seob, features = c('AR'), pt.size = 2.5)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        plot.title = element_text(hjust = 0.5, size = 16))+
  ggtitle('AR')+
  NoLegend()
p6 <-  FeaturePlot(GSE140440.seob, features = c('SYP'), pt.size = 2.5)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        plot.title = element_text(hjust = 0.5, size = 16))+
  ggtitle('SYP')+
  NoLegend()
pdf(file = './pcdata/paper/plot/section4/CD44.AR.SYP.supp.scRNA.pdf', width = 15, height = 6)
p5|p6|p4
dev.off()


# Supp Fig4f-g ------------------------------------------------------------

# Supp 4g NEPC-----------------------------------------------------------------


load(file = './pcdata/paper/data/section4/CRPC.confident.seob.RData')
CRPC.tumor.seob <- subset(CRPC.confident.seob, tumor_subtype != 'normal')
stem.cell.pesudobulk <- AggregateExpression(
  CRPC.tumor.seob,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "tumor_subtype",
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

df           <- ccle.stem.cor[, 'NEPC', drop = FALSE]
df$name      <- rownames(df)
df$prostate  <- ifelse(rownames(df) %in% pc.cell.line.name, 
                       "prostate", 
                       "others")
colnames(df)[1] <- "correlation"

df$rank     <- rank(df$correlation)
df_ranked   <- df[order(df$rank), ]

#--- B. plot ---
highlight_points <- df_ranked[df_ranked$rank %in% c(1018, 1019), ]
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
    
  )+ggtitle('NEPC')
ggsave('./pcdata/paper/plot/section4/NEPC.pdf', width = 8, height = 8)


#4f ARPC --------------------------------------------------------------------

stem.cell.pesudobulk <- AggregateExpression(
  CRPC.confident.seob,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "tumor_subtype",
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

df           <- ccle.stem.cor[, 'ARPC', drop = FALSE]
df$name      <- rownames(df)
df$prostate  <- ifelse(rownames(df) %in% pc.cell.line.name, 
                       "prostate", 
                       "others")
colnames(df)[1] <- "correlation"

df$rank     <- rank(df$correlation)
df_ranked   <- df[order(df$rank), ]

#--- B. plot ---
highlight_points <- df_ranked[df_ranked$rank %in% c(1018, 1019), ]
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
    
  )+ggtitle('ARPC')
ggsave('./pcdata/paper/plot/section4/ARPC.pdf', width = 8, height = 8)


# Supp 4h -----------------------------------------------------------------

# ARPC-vcap ---------------------------------------------------------------
load(file = './pcdata/paper/data/section1/CCLE.log2.tpm.matrix.symbol.RData')
load("~/met_pc_cell_line/pcdata/CCLE.rna.seq.marker.gene.1000.symbol.RData")
ARPC_seob <- subset(CRPC.confident.seob, tumor_subtype =='ARPC')

gene.id <- intersect(rownames(ARPC_seob), rownames(CCLE.log2.tpm.matrix.symbol))
gene.id <- intersect(gene.id, CCLE.rna.seq.marker.gene.1000.symbol)


ARPC.vcap.cor <- cor(ARPC_seob@assays$RNA@data[gene.id, ] %>% as.matrix() ,CCLE.log2.tpm.matrix.symbol[gene.id, 'VCAP_PROSTATE', drop=F], method = 'spearman')

ARPC.vcap.cor <- ARPC.vcap.cor %>% reshape2::melt() %>% as.data.frame()
colnames(ARPC.vcap.cor) <- c('barcode', 'type','value')

# MSPC-PC3 ----------------------------------------------------------------


MSPC_seob     <- subset(CRPC.confident.seob, tumor_subtype =='MSPC')

gene.id       <- intersect(rownames(MSPC_seob), rownames(CCLE.log2.tpm.matrix.symbol))
gene.id       <- intersect(gene.id, CCLE.rna.seq.marker.gene.1000.symbol)

MSPC.pc3.cor  <- cor(MSPC_seob@assays$RNA@data[gene.id, ] %>% as.matrix() ,CCLE.log2.tpm.matrix.symbol[gene.id, 'PC3_PROSTATE', drop=F], method = 'spearman')


MSPC.pc3.cor            <- MSPC.pc3.cor %>% reshape2::melt() %>% as.data.frame()
colnames(MSPC.pc3.cor)  <- c('barcode', 'type','value')
MSPC.ARPC.cor           <- rbind(MSPC.pc3.cor, ARPC.vcap.cor)

colnames(MSPC.ARPC.cor)[2] <- 'variable'
MSPC.ARPC.cor$variable <- factor(MSPC.ARPC.cor$variable, levels = c("VCAP_PROSTATE", "PC3_PROSTATE"))
ggplot(MSPC.ARPC.cor, aes(x = variable, y = value)) +
  geom_violin(aes(fill = variable, colour = variable), alpha = 0.5) +
  # alpha控制不透明度
  geom_boxplot(aes(colour = variable), width = 0.2)+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )+
  scale_color_manual(values = c("PC3_PROSTATE" = 'firebrick', "VCAP_PROSTATE" = 'blue'))+
  scale_fill_manual(values = c("PC3_PROSTATE" = 'firebrick', "VCAP_PROSTATE" = 'blue'))+
  stat_compare_means()+
  scale_x_discrete(labels = c('PC3_PROSTATE' = paste('PC3', sep = '\n', '(n = 2785)'),
                              'VCAP_PROSTATE' = paste('VCAP', sep = '\n', '(n = 686)')))+
  
  xlab('')+
  ylab('Transcriptone similarity')+
  theme_bw() +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),
    axis.text = element_text(size = 25, color = 'black'),
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA),
    axis.ticks.length = unit(0.4, "cm")
  )
ggsave('./pcdata/paper/plot/section4/PC3.vcap.cor.pdf', width = 8, height = 8)


# Supp 4j -----------------------------------------------------------------

#MSPC pheatmap annotation and boxplot -----------------------------------------------------
rownames(row_anno_df) <- rownames(mat)  ###矩阵行名和注释信息行名需要对应


load(file = './pcdata/paper/data/section3/CRPC.confident.seob.RData')
MSPC.seob                   <- subset(CRPC.confident.seob, tumor_subtype == 'MSPC')
norm.seob                   <- readRDS("~/met_pc_cell_line/pcdata/norm.seob.rds")

three.cell.type.seob <- subset(norm.seob, cell_type_new =='BE'| cell_type_new=='LE'|cell_type_new=='NE')
three.cell.type.aver <- AverageExpression(
  three.cell.type.seob,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "cell_type_new",
  add.ident = NULL,
  slot = 'data',
  verbose = TRUE)

three.cell.type.aver.matrix                                                     <- as.matrix(three.cell.type.aver$RNA)

###筛选基因使用三种亚型5000高变基因
three.cell.type.seob <- FindVariableFeatures(three.cell.type.seob, selection.method = "vst", nfeatures = 5000)

norm.5000.marker.gene <- head(VariableFeatures(three.cell.type.seob), 5000)

nor.data                                                        <- as.matrix(GetAssayData(MSPC.seob,layer='data'))
gg                 <- intersect(rownames(nor.data),norm.5000.marker.gene)
it.data            <- nor.data[gg,]  ##scrna expression matrix
norm.data          <- three.cell.type.aver.matrix[gg,]
cor.matrix <- foreach(it = iterators::iter(obj = it.data,by = 'column',chunksize= 100) ,.combine = 'rbind') %dopar% {
  rs         <- cor(it,norm.data,method = 'spearman')
  rs  
}
cor.matrix  ## result get  one cell - one cell line  correlation
cor.max                             <- colnames(cor.matrix)[apply(cor.matrix, 1, function(x) which(x == max(x))[1])]
table(cor.max)
cor.anno.matrix                     <- cbind(cor.matrix, cor.max)
p1 <- pheatmap(cor.matrix, show_rownames = F, cluster_rows = F, cluster_cols = F)


###### boxplot
cor.res <- cor.matrix
cor.res <- cor.res %>% reshape2::melt()
colnames(cor.res)[2] <- 'variable'
## theme
theme_niwot <- function(){
  theme_bw() +
    theme(text = element_text(family = "Helvetica"),
          axis.text = element_text(size = 16), 
          axis.title = element_text(size = 18),
          axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 18, vjust = 1, hjust = 0),
          legend.text = element_text(size = 12),          
          legend.title = element_blank(),                              
          legend.position = c(0.95, 0.15), 
          legend.key = element_blank(),
          legend.background = element_rect(color = "black", 
                                           fill = "transparent", 
                                           size = 2, linetype = "blank"))
}

p1 <- ggplot(cor.res, aes(x = variable, y = value)) +
  geom_violin(aes(fill = variable, colour = variable), alpha = 0.5) +
  # alpha控制不透明度
  geom_boxplot(aes(colour = variable), width = 0.2)+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )+
  scale_color_manual(values = c("BE" = '#D65813', "LE" = '#2CAD3F', 'NE'='#6792CD'))+
  scale_fill_manual(values = c("BE" = '#D65813', "LE" = '#2CAD3F', 'NE'='#6792CD'))+
  theme_niwot()+
  stat_compare_means(paired = T, comparisons = list(c('BE','LE')))+
  
  scale_x_discrete(labels = c('BE' = paste('Basal', sep = '\n', '(n = 2785)'),
                              'LE' = paste('Luminal', sep = '\n', '(n = 2785)'),
                              'NE' = paste('Neuroendocrine', sep = '\n', '(n = 2785)')))+
  
  xlab('')+
  ylab('Correlation')
pdf('./pcdata/paper/plot/section4/three.type.anno.cor.pdf', width = 8, height = 8)
p1
dev.off()


####heatmap
mat <- cor.matrix

row_df <- cor.anno.matrix[,'cor.max', drop=F]
row_df <- as.data.frame(row_df)
#type_colors    <- setNames(RColorBrewer::brewer.pal(2, "Set1"), unique(row_df$cor.max))
type_colors    <- list(cor.max = c('LE'='#050507' , 'BE'='#B22222'))
row_anno <- rowAnnotation(
  df = data.frame(cor.max = row_df$cor.max),
  col = type_colors,
  annotation_name_side = "top",
  show_annotation_name = TRUE
)
boxplot(mat)
boxplot( t(scale(t(mat))))
colnames(mat) <- c('Basal', 'Luminal', 'Neuroendocrine')
col.fun    <- circlize::colorRamp2(seq(-1,1, length.out = 101),colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(101))
p1 <- Heatmap(
  t(scale(t(mat))),
  name = "z_score",
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = FALSE,
  right_annotation = row_anno,
  col = col.fun,
  column_names_rot = 315
  
)
pdf('./pcdata/paper/plot/section4/heatmap.annotation.pdf', width = 7, height = 7)
p1
dev.off()




# Supp Fig 4 i ------------------------------------------------------------

p1 <- FeaturePlot(CRPC.tumor.seob, features = c('CAV2'))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        plot.title = element_text(hjust = 0.5))+
  NoLegend()

# supp4 k-l PC3 basal marker --------------------------------------------------------
load("./pcdata/GSE140440_smart-seq_pc3/GSE140440.seob.RData")
pc3.seob                    <- readRDS("./pcdata/met_pc/pc3.seob.rds")

p1 <- FeaturePlot(pc3.seob,
                  features = "TP63", pt.size = 2,
                  cols = c("#D3D3D3", "#D3D3D3"),
                  min.cutoff = 0,
                  max.cutoff = 0)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('TP63')+
  NoLegend()

pdf('./pcdata/paper/plot/section4/supp.tp63.pc3.expression.pdf', width = 5, height = 5)
p1
dev.off()


p1 <- FeaturePlot(GSE140440.seob,
                  features = "TP63", pt.size = 2,
                  cols = c("#D3D3D3", "#D3D3D3"),
                  min.cutoff = 0,
                  max.cutoff = 0)+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('TP63')+
  NoLegend()
#4m MSPC correlation with integrate ccle.pc3  -----------------------------------

load(file = './pcdata/paper/data/section4/CRPC.confident.seob.RData')
CRPC.tumor.seob <- subset(CRPC.confident.seob, tumor_subtype != 'normal')

stem.cell.pesudobulk <- AggregateExpression(
  CRPC.tumor.seob,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "tumor_subtype",
  add.ident = NULL,
  slot = "counts",
  verbose = TRUE
)

stem.cell.pesudobulk                     <- stem.cell.pesudobulk$RNA
stem.cell.pesudobulk.cpm                 <- apply(stem.cell.pesudobulk, 2, function(x) { x/sum(x)*1000000 })
stem.cell.pesudobulk.log2cpm             <- log2(stem.cell.pesudobulk.cpm +1)

stem.cell.pesudobulk.log2cpm             <- change_rownames_symbol_to_ensemble(stem.cell.pesudobulk.log2cpm)




marker.gene           <- intersect(rownames(stem.cell.pesudobulk.log2cpm),CCLE.rna.seq.marker.gene.1000)  
marker.gene           <- intersect(rownames(ccle.pc3.data.log2cpm),marker.gene) 
ccle.stem.cor         <- cor(stem.cell.pesudobulk.log2cpm[marker.gene,],ccle.pc3.data.log2cpm[marker.gene,],method='spearman')
ccle.stem.cor         <- ccle.stem.cor %>% t %>% as.data.frame()
ccle.stem.cor[c(pc.cell.line.name,'pc3.1','pc3.2','pc3.3'),]###MDA-PCA-2b PC3

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
highlight_points <- df_ranked[df_ranked$name %in% c('pc3.1','pc3.2','pc3.3'), ]
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
    data = subset(df, name %in% c(pc.cell.line.name, 'pc3.1','pc3.2','pc3.3')),
    color = "red",
    size = 4
  )+
  geom_text_repel(data = subset(df, name %in% c(pc.cell.line.name, 'pc3.1','pc3.2','pc3.3')), 
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
    
  )+ggtitle('MSPC')
ggsave('./pcdata/paper/plot/section4/MSPC.supp.pdf', width = 8, height = 8)


