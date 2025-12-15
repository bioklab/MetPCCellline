#Fig6a cMET high volcano -----------------------------------------------------------------

cMET_vs_wild_two_res <- read.csv('./pcdata/paper/data/section5/cMET_vs_wild_two_res.csv')
cMET_vs_wild_two_res$gene <- cMET_vs_wild_two_res$X
cMET_vs_wild_two_res$sig  <- ifelse(cMET_vs_wild_two_res$padj<0.05 & cMET_vs_wild_two_res$log2FoldChange > 1, 'up',
                                    ifelse(cMET_vs_wild_two_res$padj <0.05 & cMET_vs_wild_two_res$log2FoldChange < -1, 'dn', 'no'))
de_genes                  <- c('EP300', 'CREBBP', 'BMPR2','TGFBR2', 'JAK2','MET','KLK2','KLK4','KLK3')

ggplot(cMET_vs_wild_two_res, aes(x = log2FoldChange, y = -log10(padj),color=sig))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "black")+
  geom_point()+
  scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
  scale_size_continuous(range = c(0,1))+
  theme(legend.title = element_blank() )+
  geom_label_repel(data = subset(cMET_vs_wild_two_res, gene %in% de_genes),
                   max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                   aes(label = gene),
                   size = 5, 
                   color = 'black',fontface = "bold.italic",
                   segment.color = "black",   
                   segment.linetype = "solid", 
                   segment.size = 0.5 ,
                   box.padding = 3,        
                   point.padding = 0.2 ) +
  xlab("log2FC")+
  ylab("-log10 p-adj")+
  theme_bw()+
  ggtitle('')+
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),     # 坐标轴标题加粗
    axis.text = element_text(size = 14, face = "bold"),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.2, "cm"),         # 增加刻度线长度
    axis.ticks = element_line(size = 1.2)        # 增加刻度线粗细
  )

ggsave('./pcdata/paper/plot/section5/cMET.volcano.pdf', width = 8, height = 8)

#Fig6b MET JAK/STAT score ------------------------------------------------------
load(file = './pcdata/paper/data/section3/MET.df.RData')
tmp            <- dplyr::select(MET.df, c('median.pc3.cor', 'JAK_STAT' ,'basal.sig','CD44','KLK2','KLK4','type'))
tmp <- tmp[c(-6,-7),] %>% dplyr::select(c('JAK_STAT','type')) %>% reshape2::melt() 
tmp$type <- factor(tmp$type, levels = c('wild_type','cMET_high'))


tmp$rank <- rank(tmp$value)

tmp %>% ggplot(aes(x= type, y = value, color = type))+
  geom_point(size = 7)+
  scale_color_manual(values =c(wild_type="Blue",cMET_high="Red"))+
  ylab('JAK/STAT score')+
  xlab('')+
  theme_bw()+
  scale_x_discrete(labels=(c('wild_type'='Control',
                             c('cMET_high'='VCAP-MET-OE'))))+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )

####engineered correlation with MSPC

#Fig6c 单个MSPC 细胞与 工程细胞算相似性 -----------------------------------------------------

load(file = './pcdata/paper/data/section5/engineered.cell.line.RData')
load(file = './pcdata/paper/data/section3/CRPC.confident.seob.RData')
engineered.cell.line.tpm   <- all.sample.tpm[,c(lncap.ko.id, lncap.enza.id, vcap.oe.id, 'PC3_PROSTATE')]
engineered.cell.line.tpm.symbol <- change_rownames_ensemble_to_symbol(engineered.cell.line.tpm)

MSPC_seob     <- subset(CRPC.confident.seob, tumor_subtype =='MSPC')


engineered.cell.line.tpm.symbol %>% rownames()
gene.id       <- intersect(rownames(MSPC_seob), rownames(engineered.cell.line.tpm.symbol))
gene.id       <- intersect(gene.id, CCLE.rna.seq.marker.gene.1000.symbol)



MSPC.stem.cor    <- cor(MSPC_seob@assays$RNA@data[gene.id, ] %>% as.matrix() ,engineered.cell.line.tpm.symbol[gene.id, , drop=F], method = 'spearman')

#### 按照median进行排序
engineered.order <- apply(MSPC.stem.cor, 2, median) %>% sort(decreasing = T) %>% names()
MSPC.stem.cor    <- MSPC.stem.cor[,engineered.order]
head(MSPC.stem.cor)

MSPC.stem.cor <- MSPC.stem.cor %>% reshape2::melt()
head(MSPC.stem.cor)
colnames(MSPC.stem.cor)[2] <- 'engineered.cell.line'
MSPC.stem.cor %>% ggplot(aes(x = engineered.cell.line, y = value))+
  
  geom_violin(aes(fill = engineered.cell.line, colour = engineered.cell.line), alpha = 0.5) +
  # alpha控制不透明度
  geom_boxplot(aes(colour = engineered.cell.line), width = 0.2)+
  scale_color_manual(values = c("PC3_PROSTATE" = "firebrick", 'MET_OE2' = '#FBCE6A', 'MET_OE1'='#FBCE6A', 'LNCAPKO2'='#974F9F',
                                'LNCAPKO1'='#974F9F','LNCaP.enza.1'='#354898','LNCaP.enza.2'='#354898','LNCaP.enza.3'='#354898'))+
  scale_fill_manual(values = c("PC3_PROSTATE" = "firebrick", 'MET_OE2' = '#FBCE6A', 'MET_OE1'='#FBCE6A', 'LNCAPKO2'='#974F9F',
                               'LNCAPKO1'='#974F9F','LNCaP.enza.1'='#354898','LNCaP.enza.2'='#354898','LNCaP.enza.3'='#354898'))+
  scale_x_discrete(labels = c('PC3_PROSTATE' = paste('PC3', sep = '\n', '(n = 2785)'),
                              'MET_OE2' = paste('VCaP-MET-OE1', sep = '\n', '(n = 2785)'),
                              'MET_OE1' = paste('VCaP-MET-OE2', sep = '\n', '(n = 2785)'),
                              'LNCAPKO2'=paste('LNCaP-TP53/RB1-KO1', sep = '\n','(n= 2785)'),
                              'LNCAPKO1'=paste('LNCaP-TP53/RB1-KO2', sep = '\n','(n= 2785)'),
                              'LNCaP.enza.1'=paste('LNCaP-enza1',sep = '\n','(n = 2785)'),
                              'LNCaP.enza.2'=paste('LNCaP-enza2',sep = '\n','(n = 2785)'),
                              'LNCaP.enza.3'=paste('LNCaP-enza3',sep = '\n','(n = 2785)')))+
  
  theme_niwot()+
  xlab('')+
  ylab('Transcriptome similarity')+
  RotatedAxis()+
  theme(legend.position = 'none',
        axis.title = element_text(size = 30, color = 'black'),
        axis.text.x  = element_text(size = 30, color = 'black'),
        axis.text.y  = element_text(size = 30, color = 'black'))+
  stat_compare_means(paired = T, comparisons = list(c('PC3_PROSTATE','MET_OE2')))

ggsave('./pcdata/paper/plot/section5/engineered.pdf', width = 28, height = 12)


#Fig6d engineered cell line lineage analysis -------------------------------------------------


load(file = './pcdata/paper/data/section5/engineered.cell.line.RData')

engineered.cell.line.tpm   <- all.sample.tpm[,c(lncap.ko.id, lncap.enza.id, vcap.oe.id, 'PC3_PROSTATE')]


load( file = './pcdata/paper/data/section5/three.cell.type.seob.RData')

stem.cell.pesudobulk <- AggregateExpression(
  three.cell.type.seob,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "cell_type_new",
  add.ident = NULL,
  slot = "counts",
  verbose = TRUE
)

stem.cell.pesudobulk                     <- stem.cell.pesudobulk$RNA
stem.cell.pesudobulk.cpm                 <- apply(stem.cell.pesudobulk, 2, function(x) { x/sum(x)*1000000 })
stem.cell.pesudobulk.log2cpm             <- log2(stem.cell.pesudobulk.cpm +1)

three.cell.type.seob <- FindVariableFeatures(three.cell.type.seob, selection.method = "vst", nfeatures = 5000)

norm.5000.marker.gene <- head(VariableFeatures(three.cell.type.seob), 5000)


engineered.cell.line.tpm.symbol <- change_rownames_ensemble_to_symbol(engineered.cell.line.tpm)

nor.data                                                        <- engineered.cell.line.tpm.symbol

gg                 <- intersect(rownames(nor.data),norm.5000.marker.gene)
it.data            <- nor.data[gg,]  ##scrna expression matrix
norm.data          <- stem.cell.pesudobulk.log2cpm[gg,]
cor.matrix <- foreach(it = iterators::iter(obj = it.data,by = 'column',chunksize= 100) ,.combine = 'rbind') %dopar% {
  rs         <- cor(it,norm.data,method = 'spearman')
  rs  
}
cor.matrix  ## result get  one cell - one cell line  correlation
cor.max                             <- colnames(cor.matrix)[apply(cor.matrix, 1, function(x) which(x == max(x))[1])]
cor.anno.matrix                     <- cbind(cor.matrix, cor.max)
cor.anno.matrix <- as.data.frame(cor.anno.matrix)
cor.anno.matrix <- cor.anno.matrix[engineered.order,]
cor.anno.matrix$name <- c('PC3','VCaP-MET-OE1','VCaP-MET-OE2','LNCap-TP53/RB1-KO1','LNCap-TP53/RB1-KO2','LNCaP-enza1','LNCaP-enza2','LNCaP-enza3')

p1 <- pheatmap(cor.matrix, show_rownames = F, cluster_rows = F, cluster_cols = F)



####heatmap
mat <- cor.anno.matrix[,1:3]


row_df <- cor.anno.mamatrow_df <- cor.anno.matrix[,c('cor.max','name'), drop=F]
row_df <- as.data.frame(row_df)

type_colors    <- list(cor.max = c('LE'='#2CAD3F' , 'BE'='#D65813'),
                       name = c('PC3'="firebrick", 'VCap-MET-OE'='#FBCE6A', 'LNCap-TP53/RB1-KO'='#974F9F', 'LNCap-enza'='#354898'))
row_anno <- rowAnnotation(
  df = data.frame(cor.max = row_df$cor.max),
  
  col = type_colors,
  annotation_name_side = "top",
  show_annotation_name = TRUE
)
mat[] <- lapply(mat, as.numeric)
boxplot(mat)
mat <- as.matrix(mat)
colnames(mat) <- c('Basal', 'Luminal', 'Neuroendocrine')
col.fun    <- circlize::colorRamp2(seq(0.3,0.7, length.out = 101),colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(101))
rownames(mat) <- c('PC3','VCaP-MET-OE1','VCaP-MET-OE2','LNCaP-TP53/RB1-KO1','LNCaP-TP53/RB1-KO2','LNCaP-enza1','LNCaP-enza2','LNCaP-enza3')
p1 <- Heatmap(
  mat,
  name = "Correlation",
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = T,
  right_annotation = row_anno,
  row_names_side = "left",
  col = col.fun,
  column_names_rot = 315
  
)

p1
pdf('./pcdata/paper/plot/section5/lineage.annotation.pdf', width = 8, height = 8)
p1
dev.off()




#Fig6e 单个MSPC恶性细胞与organoid算相似性 -------------------------------------------------



load(file = './pcdata/paper/data/section5/stem.organoid.ccle.RData')
load(file = './pcdata/paper/data/section5/stem.organoid.RData')
## ccle stem bulk

stem.bulk.organoid.tpm
stem.id <- colnames(stem.bulk.organoid.tpm)
MSPC_seob     <- subset(CRPC.confident.seob, tumor_subtype =='MSPC')

stem.organoid.ccle.symbol <- change_rownames_ensemble_to_symbol(stem.organoid.ccle)

gene.id       <- intersect(rownames(MSPC_seob), rownames(stem.organoid.ccle.symbol))
gene.id       <- intersect(gene.id, CCLE.rna.seq.marker.gene.1000.symbol)

organoid.pc3.id <- c(colnames(stem.bulk.organoid.tpm),'PC3_PROSTATE')

MSPC.stem.cor  <- cor(MSPC_seob@assays$RNA@data[gene.id, ] %>% as.matrix() ,stem.organoid.ccle.symbol[gene.id,organoid.pc3.id], method = 'spearman')
organoid.order <- apply(MSPC.stem.cor, 2, median) %>% sort(decreasing = T) %>% names()

MSPC.stem.cor <- MSPC.stem.cor[,organoid.order]
head(MSPC.stem.cor)

MSPC.stem.cor <- MSPC.stem.cor %>% reshape2::melt()
head(MSPC.stem.cor)
colnames(MSPC.stem.cor)[2] <- 'organoid'

MSPC.stem.cor$organoid <- factor(MSPC.stem.cor$organoid)

# 获取 levels
organoid_levels <- levels(MSPC.stem.cor$organoid)

# 创建颜色映射：PC3_PROSTATE 为蓝色，其余为红色
color_map <- setNames(
  ifelse(organoid_levels == "PC3_PROSTATE", "firebrick", 'steelblue'),
  organoid_levels
)

organoid_median_order <- MSPC.stem.cor %>%
  group_by(organoid) %>%
  summarise(median_value = median(value)) %>%
  arrange(desc(median_value)) %>%
  pull(organoid)

# 按照中位数从高到低重新设定 factor 顺序
MSPC.stem.cor$organoid <- factor(MSPC.stem.cor$organoid, levels = organoid_median_order)



MSPC.stem.cor %>% ggplot(aes(x = organoid, y = value))+
  
  geom_violin(aes(fill = organoid, colour = organoid), alpha = 0.5) +
  # alpha控制不透明度
  geom_boxplot(aes(colour = organoid), width = 0.2)+
  scale_color_manual(values = color_map) +
  scale_fill_manual(values = color_map)+
  scale_x_discrete(labels = c('PC3_PROSTATE' = paste('PC3', sep = '\n', '(n = 2785)'),
                              'MSKPCa12' = paste('MSKPCa12', sep = '\n', '(n = 2785)'),
                              'MSKPCa9' = paste('MSKPCa9', sep = '\n', '(n = 2785)'),
                              'MSKPCa13'=paste('MSKPCa13', sep = '\n','(n= 2785)'),
                              'MSKPCa18'=paste('MSKPCa18', sep = '\n','(n= 2785)'),
                              'MSKPCa8'=paste('MSKPCa8',sep = '\n','(n = 2785)'),
                              'MSKPCa3'=paste('MSKPCa3',sep = '\n','(n = 2785)'),
                              'MSKPCa11'=paste('MSKPCa11',sep = '\n','(n = 2785)'),
                              'MSKPCa17'=paste('MSKPCa17',sep = '\n','(n = 2785)'),
                              'MSKPCa15'=paste('MSKPCa15',sep = '\n','(n = 2785)'),
                              'MSKPCa20'=paste('MSKPCa20',sep = '\n','(n = 2785)')))+
  
  theme_niwot()+
  xlab('')+
  ylab('Transcriptome similarity')+
  RotatedAxis()+
  theme(legend.position = 'none',
        axis.title = element_text(size = 30, color = 'black'),
        axis.text.x  = element_text(size = 30, color = 'black'),
        axis.text.y  = element_text(size = 30, color = 'black'))+
  stat_compare_means(paired = T, comparisons = list(c('PC3_PROSTATE','MSKPCa12'),
                                                    c('PC3_PROSTATE','MSKPCa9'),
                                                    c('PC3_PROSTATE','MSKPCa13'),
                                                    c('PC3_PROSTATE','MSKPCa18')))

ggsave('./pcdata/paper/plot/section5/organoid.pdf', width = 28, height = 12)





#Fig6f organoid lineage analysis -----------------------------------------------



load( file = './pcdata/paper/data/section5/three.cell.type.seob.RData')


stem.cell.pesudobulk <- AggregateExpression(
  three.cell.type.seob,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "cell_type_new",
  add.ident = NULL,
  slot = "counts",
  verbose = TRUE
)

stem.cell.pesudobulk                     <- stem.cell.pesudobulk$RNA
stem.cell.pesudobulk.cpm                 <- apply(stem.cell.pesudobulk, 2, function(x) { x/sum(x)*1000000 })
stem.cell.pesudobulk.log2cpm             <- log2(stem.cell.pesudobulk.cpm +1)

three.cell.type.seob <- FindVariableFeatures(three.cell.type.seob, selection.method = "vst", nfeatures = 5000)

norm.5000.marker.gene <- head(VariableFeatures(three.cell.type.seob), 5000)


stem.organoid.ccle.symbol <- change_rownames_ensemble_to_symbol(stem.organoid.ccle)
nor.data                                                        <- stem.organoid.ccle.symbol

gg                 <- intersect(rownames(nor.data),norm.5000.marker.gene)
it.data            <- nor.data[gg,c(colnames(stem.bulk.organoid.tpm), 'PC3_PROSTATE')]  ##scrna expression matrix
norm.data          <- stem.cell.pesudobulk.log2cpm[gg,]
cor.matrix <- foreach(it = iterators::iter(obj = it.data,by = 'column',chunksize= 100) ,.combine = 'rbind') %dopar% {
  rs         <- cor(it,norm.data,method = 'spearman')
  rs  
}
cor.matrix  ## result get  one cell - one cell line  correlation
cor.max                             <- colnames(cor.matrix)[apply(cor.matrix, 1, function(x) which(x == max(x))[1])]
cor.anno.matrix                     <- cbind(cor.matrix, cor.max)
cor.anno.matrix <- as.data.frame(cor.anno.matrix)
organoid_median_order <- levels(organoid_median_order) %>% c
cor.anno.matrix <- cor.anno.matrix[organoid_median_order,]




####heatmap
mat <- cor.anno.matrix[,1:3]


row_df <- cor.anno.mamatrow_df <- cor.anno.matrix[,c('cor.max'), drop=F]
row_df <- as.data.frame(row_df)

type_colors    <- list(cor.max = c('LE'='#2CAD3F' , 'BE'='#D65813')
)
row_anno <- rowAnnotation(
  df = data.frame(cor.max = row_df$cor.max),
  
  col = type_colors,
  annotation_name_side = "top",
  show_annotation_name = TRUE
)
mat[] <- lapply(mat, as.numeric)
boxplot(mat)
mat <- as.matrix(mat)
colnames(mat) <- c('Basal', 'Luminal', 'Neuroendocrine')
col.fun    <- circlize::colorRamp2(seq(0.3,0.7, length.out = 101),colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(101))

p1 <- Heatmap(
  mat,
  name = "Correlation",
  cluster_rows = F,
  cluster_columns = F,
  show_row_names = T,
  right_annotation = row_anno,
  row_names_side = "left",
  col = col.fun,
  column_names_rot = 315
  
)
p1
pdf('./pcdata/paper/plot/section5/lineage.annotation.pdf', width = 8, height = 8)
p1
dev.off()



# Fig6g-h pc3 organoid basal marker expression -------------------------------------------------


load(file = './pcdata/paper/data/section5/stem.organoid.RData')
## ccle stem bulk
CCLE.log2.tpm.matrix
stem.bulk.organoid.tpm
gene.id <- intersect(rownames(CCLE.log2.tpm.matrix), rownames(stem.bulk.organoid.tpm))
ccle.tpm <- CCLE.log2.tpm.matrix[gene.id,]
stem.tpm <- stem.bulk.organoid.tpm[gene.id,]
identical(rownames(ccle.tpm), rownames(stem.tpm))
ccle.stem         <- cbind(stem.tpm, ccle.tpm)
stem.pc3          <- ccle.stem[, c(colnames(stem.bulk.organoid.tpm), 'PC3_PROSTATE')]
stem.pc3.data     <- 2^stem.pc3-1
stem.pc3.data.cpm <- apply(stem.pc3.data, 2, function(x) { x/sum(x)*1000000 })
stem.pc3.data.log2cpm             <- log2(stem.pc3.data.cpm +1)

stem.pc3.data.log2cpm.symbol <- change_rownames_ensemble_to_symbol(stem.pc3.data.log2cpm)
df <- stem.pc3.data.log2cpm.symbol['TP63',] %>% sort() %>% as.data.frame()
colnames(df) <- 'expression'
df$sample <- rownames(df)

df <- df[order(df$expression, decreasing = T),]
df$sample <- factor(df$sample, levels = rev(df$sample))

ggplot(df, aes(x = sample, y = expression)) +
  geom_bar(stat = "identity", fill = "#D8CDEB") +
  
  coord_flip() +
  labs(title = "",
       x = "",
       y = "")+
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 20, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 20, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.25, "cm"),         # 增加刻度线长度
    
  )
ggsave('./pcdata/paper/plot/section5/TP63_expr.pdf', width = 8, height = 6)

df <- stem.pc3.data.log2cpm.symbol['KRT5',] %>% sort() %>% as.data.frame()
colnames(df) <- 'expression'
df$sample <- rownames(df)

df <- df[order(df$expression, decreasing = T),]
df$sample <- factor(df$sample, levels = rev(df$sample))

ggplot(df, aes(x = sample, y = expression)) +
  geom_bar(stat = "identity", fill = "#D8CDEB") +
  
  coord_flip() +
  labs(title = "",
       x = "",
       y = "")+
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 20, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 20, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.25, "cm"),         # 增加刻度线长度
    
  )
ggsave('./pcdata/paper/plot/section5/KRT5_expr.pdf', width = 8, height = 6)




