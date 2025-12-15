# TCGA SU2C CCLE mutation burden Fig.2a------------------------------------------
load(file = './pcdata/paper/data/section1/mutation.burden.RData')
TCGA.df <- data.frame(x = rank(TCGA.mutation.burden), y = log2(TCGA.mutation.burden))
TCGA.df <- data.frame(id = names(TCGA.mutation.burden), mutation.burden = log2(TCGA.mutation.burden))
TCGA.df$type <- 'TCGA'


SU2C.df <- data.frame(x = rank(SU2C.mutation.burden), y = log2(SU2C.mutation.burden))
SU2C.df <- data.frame(id = names(SU2C.mutation.burden), mutation.burden = log2(SU2C.mutation.burden))
SU2C.df$type <- 'SU2C'


CCLE.df      <- data.frame(x = rank(CCLE.mutation.burden), y = log2(CCLE.mutation.burden))
CCLE.df      <- data.frame(id = names(CCLE.mutation.burden), mutation.burden = log2(CCLE.mutation.burden))
CCLE.df$type <- 'CCLE'
mutation.burden.df      <- rbind(TCGA.df, SU2C.df, CCLE.df)
mutation.burden.df$rank <- rank(mutation.burden.df$mutation.burden)



mutation.burden.df$type2 <- ifelse(mutation.burden.df$mutation.burden > 9 ,paste0(mutation.burden.df$type,'hyper'), 'other')
type_colors <- c(
  "CCLEhyper" = "#377EB9",     # 蓝色
  "SU2Chyper" = "#4DAE48",   # 绿色
  "TCGAhyper" = "#974F9F",     # 紫色
  'other'='black'
)


ggplot(mutation.burden.df, aes(x = rank, y = mutation.burden, color = type2))+
  geom_point(size = 3)+ 
  geom_text_repel(data = subset(mutation.burden.df, type =='CCLE' & mutation.burden>9),
                  aes(label = id),
                  nudge_y = 0.1)+
  xlab('Rank') + ylab('Mutation burden')+
  geom_hline(yintercept = 9, color = "black")+
  scale_color_manual(values = type_colors) +  # 这里指定颜色
  theme_bw()+
  theme(
    
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )
ggsave('./pcdata/paper/plot/section2/three.dataset.mutation.burden.pdf', width = 12,height = 10)


#Fig. 2b MMR mutation profile ---------------------------------------------------------
load(file = './pcdata/paper/data/section1/CCLE.SU2C.TCGA.mutation.matrix.RData')

MMR.mutation.matrix <- cbind(TCGA.mutation.matrix[MMR.gene,],SU2C.mutation.matrix[MMR.gene,],CCLE.mutation.matrix[MMR.gene,])
MMR.mut.sample      <- colnames(MMR.mutation.matrix)[apply(MMR.mutation.matrix,2,sum) >0 ]  ####MMR 修复基因突变的样本
MMR.gene            <- c('MLH1','MLH3','MSH2','MSH6','MSH3','PMS1','PMS2','POLE','POLD1') # https://www.nature.com/articles/ncomms15180/, Fig 1c
TCGA.MSI.sample     <- names(TCGA.mutation.burden)[log2(TCGA.mutation.burden) > 9]
SU2C.MSI.sample     <- names(SU2C.mutation.burden)[log2(SU2C.mutation.burden) > 9]
CCLE.MSI.sample     <- names(CCLE.mutation.burden)[log2(CCLE.mutation.burden) > 9]
MSI.sample          <- c(TCGA.MSI.sample,SU2C.MSI.sample,CCLE.MSI.sample)     ####### 超突变的样本
MSI.no.sample       <- setdiff(colnames(MMR.mutation.matrix),MSI.sample)  ####不是超突变的样本
MMR.no.mut.sample   <- setdiff(colnames(MMR.mutation.matrix),MMR.mut.sample)   ######修复基因没突变的样本



MSI.sample.id        <- data.frame(id = MSI.sample, type = rep('MSI',length(MSI.sample)))
MSI.no.sample.id     <- data.frame(id = MSI.no.sample, type = rep('MSI.no.sample', length(MSI.no.sample)))
MSI.meta             <- rbind(MSI.sample.id, MSI.no.sample.id)

MMR.mutation.matrix  <- cbind(TCGA.mutation.matrix[MMR.gene,],SU2C.mutation.matrix[MMR.gene,],CCLE.mutation.matrix[MMR.gene,])
MMR.mutation.matrix  <- MMR.mutation.matrix[,MSI.meta$id]





mutation.profile <- MMR.mutation.matrix
col.df.pro <- MSI.meta
col.annotation    <-  HeatmapAnnotation(type=col.df.pro$type, 
                                        col = list(type = c("MSI" =  "red", "MSI.no.sample" = "blue")) ,
                                        show_annotation_name = F,## 不显示type
                                        show_legend = F
)


# 定义 alter_fun
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h, gp = gpar(fill = "#CCCCCC", col = NA))
  },
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w - unit(0.5, "mm"), h * 0.8, gp = gpar(fill = "black", col = NA))
  }
)

# 定义颜色
col <- c("MUT" = "black", "background" = "#CCCCCC")
sample_order <- col.df.pro$id
mutation.profile <- mutation.profile[,col.df.pro$id]
identical(colnames(mutation.profile), col.df.pro$id)
# 绘制 OncoPrint
p1 <- oncoPrint(
  mutation.profile,
  get_type = function(x) {ifelse(x == 1, "MUT", "background")},
  row_names_side = "left",          # 行名放置在左侧
  show_heatmap_legend = F,      # 隐藏图例
  alter_fun = alter_fun,            # 突变样式函数
  col = col,                        # 颜色定义
  show_column_names = FALSE,        # 不显示列名
  show_row_names = T, # 不显示行名
  top_annotation = col.annotation,  # 顶部注释
  right_annotation = NULL,          # 确保右侧无条形图
  show_pct = FALSE,                 # 禁用百分比显示
  row_names_gp = gpar(fontsize = 40, fontface = 'bold'),
  column_order = sample_order,
  row_order = gene.order
)





MSI.sample.id        <- data.frame(id = MSI.sample, type = rep('MSI',length(MSI.sample)))
MSI.no.sample.id     <- data.frame(id = MSI.no.sample, type = rep('MSI.no.sample', length(MSI.no.sample)))




MSI.mutation.fequency      <- apply(MMR.mutation.matrix[,MSI.sample], 1, sum) / length(MSI.sample)


MSI.no.mutation.fequency   <- apply(MMR.mutation.matrix[,MSI.no.sample], 1, sum) / length(MSI.no.sample)
mutation.frequency         <- cbind(MSI.mutation.fequency, MSI.no.mutation.fequency)

col.df                  <- c('MSI','MSI.no.sample')

gene.order              <-  MSI.mutation.fequency %>% sort(decreasing = T) %>% names()

names(col.df)    <- colnames(mutation.frequency)
col.df           <- as.data.frame(col.df)
colnames(col.df) <- 'type'
col.ha <- HeatmapAnnotation(type = col.df$type,
                            col = list(type = c("MSI" =  "red", "MSI.no.sample" = "blue")),
                            show_legend = T,
                            show_annotation_name = F,### rm type 
                            annotation_name_gp = gpar(fontface = "bold"),
                            annotation_legend_param = list(
                              title_gp = gpar(fontsize = 20, fontface = 'bold'), ## 图例字体
                              labels_gp = gpar(fontsize = 20),
                              legend_height = unit(10, "cm"), ###调整col annotation 图例大小
                              legend_width = unit(5, "cm")
                            ))
p2 <- Heatmap(mutation.frequency,
              show_column_dend = FALSE, 
              col=colorRamp2(c(0,1),c('blue','red')),
              row_names_gp = gpar(fontsize = 40,fontface='bold'),
              width = unit(4, "cm"),
              show_heatmap_legend = T,##显示图例
              top_annotation=col.ha,
              show_column_names = F,
              cluster_columns = F,
              cluster_rows = F,
              column_order = c("MSI.mutation.fequency", "MSI.no.mutation.fequency"), # 强制指定列顺序
              row_order = gene.order,
              name = 'Mutation \nfrequency',
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 20, fontface = "bold"), # 图例标题字体
                labels_gp = gpar(fontsize = 20),                  # 图例标签字体
                legend_height = unit(10, "cm"),                    # 图例高度
                legend_width = unit(5, "cm")                      # 图例宽度
              ))





pdf('./pcdata/paper/plot/section2/mutation.profile.MMR.pdf', width = 30, height = 15)
p1+p2
dev.off()




# Fig.2c MSI de analysis volcano -------------------------------------------------



load('./TCGA/Prostate Adenocarcinoma.RData')  
load(file = './pcdata/paper/data/section1/mutation.burden.RData')
load('./processed.data/organize.mutation.data.RData')
TCGA.MSI.patient    <- names(TCGA.mutation.burden)[log2(TCGA.mutation.burden) >= 9]
TCGA.MSI.sample     <- PC.maf.data.list$TCGA$Tumor_Sample_Barcode[match(x= TCGA.MSI.patient,table = PC.maf.data.list$TCGA$dcast.id)]
TCGA.MSI.sample     <- intersect(TCGA.MSI.sample, colnames(log2.fpkm.matrix))  ###TCGA 超突变的原发癌样本
TCGA.non.MSI.sample <- setdiff(colnames(log2.fpkm.matrix),TCGA.MSI.sample )    #### TCGA 不是超突变的原发癌样本

source('./code/BioKLab.util.R')
de.res <- perform.DE.analysis.between.TRE.and.CON(CON.log2.read.count.matrix = log2.read.count.matrix[,TCGA.non.MSI.sample],
                                                  TRE.log2.read.count.matrix = log2.read.count.matrix[,TCGA.MSI.sample],
                                                  CON.log2.tpm.matrix        = log2.tpm.matrix[,TCGA.non.MSI.sample],
                                                  TRE.log2.tpm.matrix        = log2.tpm.matrix[,TCGA.MSI.sample]
)   ####得到超突变样本和非超突变样本差异表达基因

up.gene <- rownames(de.res)[de.res$log2FoldChange > 1  & de.res$padj < 0.05] # DAVIDE analysis shows cell cycle - related process elevated, and T cell abundance elevated 
dn.gene <- rownames(de.res)[de.res$log2FoldChange < -1 & de.res$padj < 0.05] # Gamma-delta T cell marker TRGV9 down-regulated, interesting! 
de.res  <- change_rownames_ensemble_to_symbol(de.res)
de.res$gene <- rownames(de.res)
de.res$sig <- ifelse(de.res$padj<0.05 & de.res$log2FoldChange > 1, 'up',
                     ifelse(de.res$padj <0.05 & de.res$log2FoldChange < -1, 'dn', 'no'))

de_genes <- c('CD3D', 'CD3E', 'CD8A','GZMA','GZMB','PRF1','PDCD1','CTLA4')

ggplot(de.res, aes(x = log2FoldChange, y = -log10(padj),color=sig))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "black")+
  geom_point()+
  scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
  scale_size_continuous(range = c(0,1))+
  theme(legend.title = element_blank() )+
  geom_label_repel(data = subset(de.res, gene %in% de_genes),
                   max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                   aes(label = gene),
                   size = 8, 
                   color = 'black',fontface = "bold.italic",
                   segment.color = "black",   
                   segment.linetype = "solid", 
                   segment.size = 0.5 ,
                   box.padding = 3,        
                   point.padding = 0.2 ) +
  xlab("log2FC")+
  ylab("-log10 p-adj")+
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 30),
    axis.title = element_text(size = 30, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 30, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )
ggsave('./pcdata/paper/plot/section2/de.hyper.non.volcano.pdf', width = 10, height = 10)


# Fig.2d hallmark ----------------------------------------------------
hall <- read.delim(file = './pcdata/paper/data/section2/MSigDB_Hallmark_2020_table (3).txt') 
de.enrich.result <- hall %>% 
  separate(col = 'Overlap', into = c('count', 'gene_set_count'), sep = '/', remove =F) 

de.enrich.result$count <- as.numeric(de.enrich.result$count)
de.enrich.result <- de.enrich.result %>% mutate(gene_ratio = de.enrich.result$count/125)

de.enrich.result <- de.enrich.result %>% filter(Adjusted.P.value < 0.05)

de.enrich.result$Term <- gsub("\\(GO:\\d+\\)", "", de.enrich.result$Term)
de.enrich.result <- de.enrich.result[order(de.enrich.result$Adjusted.P.value),]
#write.csv(hall, file = "hall.de.enrich.result.csv", row.names = FALSE)


ggplot(de.enrich.result, aes(x = gene_ratio, y = reorder(Term, count,decreasing = F))) +
  geom_point(aes(size = count, fill = P.value), shape =21) +
  scale_fill_gradient(low = "red", high = "blue") +
  
  labs(
    
    size = "Count",
    x = "Gene Ratio",
    y = "",
    title = "MSigDB Hallmark"
  ) +
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, color = 'black'),
    axis.title.x = element_text(size = 25,  color = 'black'),
    axis.text.x = element_text(size = 22,  color = 'black'),
    axis.text.y = element_text(size = 25, color = 'black'),
    panel.border = element_rect(color = "black", linewidth = 1.2, fill = NA),
    
    axis.ticks.length = unit(0.2, "cm"),  # Increase tick mark length
    axis.ticks = element_line(size = 1)  # Increase tick mark thickness
  )

ggsave('./pcdata/paper/plot/section2/DE.MSI.non.MSI.ORA.pdf', width = 16, height = 12)












