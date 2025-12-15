#Fig4a PCA CCLE ----------------------------------------------------------------
load('./processed.data/CCLE.transcriptome.RData')
hgnc.data                    <- read.delim("./raw.data/HGNC/hgnc_complete_set.txt", stringsAsFactors=FALSE)

flag                         <- grepl(colnames(CCLE.log2.tpm.matrix), pattern ='PROSTATE')
pc.cell.line.log2.tpm.matrix <- CCLE.log2.tpm.matrix[,flag]
pc.cell.line.name            <- colnames(pc.cell.line.log2.tpm.matrix)
flag                         <- colnames(pc.cell.line.log2.tpm.matrix) != 'PRECLH_PROSTATE' ## remove no cancer 
pc.cell.line.log2.tpm.matrix <- pc.cell.line.log2.tpm.matrix[,flag]
pca.rs                       <- prcomp(pc.cell.line.log2.tpm.matrix %>% t)    
plot(pca.rs$x[,1:2])
pc1.tail.gene                <- pca.rs$rotation[,'PC1'] %>% sort %>% tail(200) %>% names
pc1.head.gene                <- pca.rs$rotation[,'PC1'] %>% sort %>% head(200) %>% names
pca.rs$rotation[,'PC1'] %>% sort(decreasing = T)
pca                          <- pca.rs$x[,1:2]
pca                          <- as.data.frame(pca)
pca$cell.line                <- gsub(x = rownames(pca), pattern = '_PROSTATE',replacement = '')

p1 <- ggplot(pca,aes(x= PC1, y = PC2, color = cell.line))+
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
  ggplot.style+
  theme(legend.position = 'none')

pdf(file = './pcdata/paper/plot/section4/pca.pdf',width = 8, height = 6)
p1
dev.off()
#Fig4b pc1 tail 200 gene 富集分析 --------------------------------------------------
`MSigDB_Hallmark_2020_table.(2)` <- read.delim("~/met_pc_cell_line/pcdata/explore_pc3/input/MSigDB_Hallmark_2020_table (2).txt")
msigdb.enrich.result <- `MSigDB_Hallmark_2020_table.(2)`

msigdb.enrich.result <- msigdb.enrich.result %>% 
  separate(col = 'Overlap',into = c('count','gene_set_count'),sep = '/',remove = F)

msigdb.enrich.result$count <- as.numeric(msigdb.enrich.result$count)
msigdb.enrich.result <- msigdb.enrich.result %>% mutate(gene_ratio=msigdb.enrich.result$count/200)


msigdb.enrich.result <- msigdb.enrich.result[order(-msigdb.enrich.result$count), ]
msigdb.enrich.result <- msigdb.enrich.result %>% filter(Adjusted.P.value < 0.05)


ggplot(msigdb.enrich.result, aes(x = gene_ratio, y = reorder(Term, count,decreasing = F))) +
  geom_point(aes(size = count, fill = P.value), shape = 21) +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(
    
    size = "count",
    x = "Gene Ratio",
    y = "",
    title = "MsigDB Hallmark"
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 38))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color="black", linetype="solid", size = 1),
        axis.text.x = element_text( size = 14,color = 'black'),
        axis.text.y = element_text(size = 14,color = 'black'),
        legend.text = element_text(color = "black",size = 13),
        legend.title = element_text(color = "black",size = 13),
        axis.ticks.x = element_line(size=1),
        axis.ticks.y = element_line(size=1),
        plot.title = element_text(color = "black",face = "bold", size = 17,hjust = 0.5))


ggsave('./pcdata/paper/plot/section4/PC1.200.ORA.pdf', height = 6.5, width = 9)

#Fig.4c-e CD44 AR SYP expression cell line--------------------------------------------------
load(file = './pcdata/paper/data/section1/CCLE.log2.tpm.matrix.symbol.RData')

expr <- CCLE.log2.tpm.matrix.symbol[c('CD44','AR','SYP'), pc.cell.line.name] %>% t %>% as.data.frame()
expr$cell.line <- gsub(x = rownames(expr), pattern = '_PROSTATE', replacement = '')

p1 <- ggplot(expr, aes(x = reorder(cell.line, CD44, decreasing = T), y = CD44))+
  geom_col(width = 0.6, fill = "#D8CDEB")+
  scale_x_discrete(labels = c('NCIH660'='NCI-H660',
                              'VCAP'='VCaP',
                              'MDAPCA2B'='MDA-PCa-2b',
                              'LNCAPCLONEFGC'='LNCaP'))+
  theme_bw()+
  ggplot.style+
  RotatedAxis()+
  ylab('CD44 expression')+xlab('')

p2 <- ggplot(expr, aes(x = reorder(cell.line, AR, decreasing = T), y = AR))+
  geom_col(width = 0.6, fill = "#D8CDEB")+
  scale_x_discrete(labels = c('NCIH660'='NCI-H660',
                              'VCAP'='VCaP',
                              'MDAPCA2B'='MDA-PCa-2b',
                              'LNCAPCLONEFGC'='LNCaP'))+
  
  ggplot.style+
  RotatedAxis()+
  ylab('AR expression')+xlab('')


p3 <- ggplot(expr, aes(x = reorder(cell.line, SYP, decreasing = T), y = SYP))+
  geom_col(width = 0.6, fill = "#D8CDEB")+
  scale_x_discrete(labels = c('NCIH660'='NCI-H660',
                              'VCAP'='VCaP',
                              'MDAPCA2B'='MDA-PCa-2b',
                              'LNCAPCLONEFGC'='LNCaP'))+
  
  ggplot.style+
  RotatedAxis()+
  ylab('SYP expression')+xlab('')

pdf(file = './pcdata/paper/plot/section3/CD44.AR.SYP.pdf', width = 15, height = 5)
p1|p2|p3
dev.off()


#Fig.4f.4h CRPC AR SYP CD44 TP63 subtype marker experssion --------------------------------------------
load(file = './pcdata/paper/data/section4/CRPC.confident.seob.RData')
CRPC.tumor.seob <- subset(CRPC.confident.seob, tumor_subtype != 'normal')
p1 <- FeaturePlot(CRPC.tumor.seob, features = c('rna_AR'))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('AR')+
  NoLegend()
p2 <- FeaturePlot(CRPC.tumor.seob, features = c('SYP'))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('SYP')+
  NoLegend()




p3 <- FeaturePlot(CRPC.tumor.seob, features = c('rna_CD44'))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('CD44')+
  NoLegend()

p4 <- DimPlot(CRPC.tumor.seob, group.by = 'tumor_subtype',label = T, 
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

pdf('./pcdata/paper/plot/section4/basal.expression.pdf',width = 20, height = 5)
p4|p1|p2|p3
dev.off()


FeaturePlot(CRPC.tumor.seob, features = c('TP63'))+
  tidydr::theme_dr(xlength = 0.2, 
                   ylength = 0.2,
                   arrow = arrow(length = unit(0.2, "inches"),type = "closed"))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 2,hjust = 0.03),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('SYP')+
  NoLegend()
marker.expr <- FetchData(CRPC.tumor.seob, vars = c('rna_AR','SYP','rna_CD44', 'TP63'))


#Fig.4g MSPC correlation with ccle  -----------------------------------
hgnc.data                 <- read.delim("/home/liuxueying/met_pc_cell_line/external.data/HGNC/hgnc_complete_set.txt", stringsAsFactors=FALSE)
change_rownames_symbol_to_ensemble <- function(expression.matrix){
  match.pos          <- match(x= rownames(expression.matrix), table = hgnc.data$symbol)
  ense.id            <- hgnc.data[match.pos,]$ensembl_gene_id
  df                 <- data.frame(ense.id = ense.id, symbol.id = rownames(expression.matrix))
  df                 <- df[complete.cases(df),]
  expression.matrix  <- expression.matrix[df$symbol.id, ,drop=F]
  if(identical(rownames(expression.matrix), df$symbol.id)){
    rownames(expression.matrix) <- df$ense.id
  }else{
    expression.matrix
  }
  return(expression.matrix)
}

load ('./processed.data/CCLE.transcriptome.RData')  ##CCLE DATA
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
    
  )+ggtitle('MSPC')
ggsave('./pcdata/paper/plot/section4/MSPC.pdf', width = 9, height = 8)



#Fig.4i MSPC pheatmap annotation and boxplot -----------------------------------------------------
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

