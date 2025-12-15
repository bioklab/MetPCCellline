
#Supp Fig.2a binominal --------------------------------------------------------
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
ggplot(mutation.burden.df, aes(x = mutation.burden)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 50,
    fill = "#7FB3D5",
    color = "black"
    
  ) +
  geom_density(
    linewidth = 0.5,
    color = "red"
  ) +
  xlab("Mutation burden") +
  ylab("Density") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13, face = "bold")
  )




#Supp Fig. 2b validate hypermutated  ------------------------------------------


load(file = './pcdata/paper/data/section1/validate.hyper.mutaed.sample.RData')

fisher.test(x)   #######修复基因突变和修复基因不突变的样本  超突变的频率不一样
rownames(x) <- c('MMR.mut.sample','MMR.no.mut.sample')
colnames(x) <- c('MSI.sample', 'MSI.no.sample')

p.val <- fisher.test(x)$p.value

data <- data.frame(
  Type = c("MMR.mut.sample", "MMR.no.mut.sample"),
  MSI.sample = c(15, 3),
  MSI.no.sample = c(40, 758)
)


data_long <- reshape2::melt(data, id.vars = "Type", variable.name = "Sample", value.name = "Count")
data_normalized <- data_long %>%
  dplyr::group_by(Sample) %>%
  mutate(Percentage = Count / sum(Count)*100 )




# 绘制归一化堆积柱状图
data_normalized$Sample <- factor(data_normalized$Sample, levels = c('MSI.no.sample','MSI.sample'))
type_colors <- c(
  "MMR.mut.sample" = "#F48892",     # 绿色
  "MMR.no.mut.sample" = '#91CAE8'      # 紫色
)
ggplot(data_normalized, aes(x = Sample, y = Percentage, fill = Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  labs(
    title = " ",
    x = "",
    y = "Percentage (%)",
    fill = "Type"
  ) +
  scale_x_discrete(labels=c('MSI.sample'='Hypermutated', 'MSI.no.sample'='Non-hypermutated')) +
  scale_fill_manual(values = type_colors)+
  annotate(
    'text',
    x = 0.5,            # x 轴位置（在两个柱中间）
    y = 105,           # y 轴位置（高于图顶端）
    label = paste0("p = ", format(p.val, scientific = TRUE)), # p 值标签
    size = 5,           # 字体大小
    color = "black",      # 字体颜色
    hjust = 0           # 水平对齐方式
  )+ggplot.style



ggsave('./pcdata/paper/plot/section2/vaidate.MSI.pdf', width = 15, height = 10)


# supp Fig. 2c ------------------------------------------------------------



imm.expr                  <- TCGA.tpm[c('CD3D', 'CD3E', 'CD8A', 'GZMA', 'GZMB', 'PRF1','PDCD1','CTLA4'),]
imm.expr                  <- imm.expr %>% t %>% as.data.frame()

imm.expr                  <- imm.expr[TCGA.meta$id,]
identical(rownames(imm.expr), TCGA.meta$id)
imm.expr$type             <- TCGA.meta$type
imm.expr                  <- imm.expr %>% melt() %>% as.data.frame()
imm.expr$type          <- factor(imm.expr$type, levels = c("non.MSI", "MSI"))

imm.expr %>% ggplot(aes(x= variable, y = value, color = type))+
  # boxplot 分组并靠近
  geom_boxplot(position = position_dodge(0.6), outlier.shape = NA, width = 0.5, lwd = 0.6) +
  
  # jitter 散点叠加
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
              size = 1.2) +
  scale_color_manual(values =c("non.MSI" =  '#8FC3DF', 'MSI' = "#E8848E")
  )+
  stat_compare_means(aes(group = type), 
                     method = "wilcox.test",
                     label = 'p.format')+
  ylab('Expression')+
  xlab('')+
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 30),
    axis.title = element_text(size = 30, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 30, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )


ggsave('./pcdata/paper/plot/section2/imm.exp.pdf', width = 15, height = 10)

#supp 2d xCell -------------------------------------------------------------------



library(xCell)
load('./TCGA/Prostate Adenocarcinoma.RData') 
log2.tpm.matrix[1:4,1:4]
TCGA.tpm <- change_rownames_ensemble_to_symbol(log2.tpm.matrix)
TCGA.tpm[1:4,1:4]
TCGA.MSI.patient    <- names(TCGA.mutation.burden)[log2(TCGA.mutation.burden) >= 9]
TCGA.MSI.sample     <- PC.maf.data.list$TCGA$Tumor_Sample_Barcode[match(x= TCGA.MSI.patient,table = PC.maf.data.list$TCGA$dcast.id)]
TCGA.MSI.sample     <- intersect(TCGA.MSI.sample, colnames(log2.tpm.matrix))  ###TCGA 超突变的原发癌样本
TCGA.non.MSI.sample <- setdiff(colnames(log2.tpm.matrix),TCGA.MSI.sample )    #### TCGA 不是超突变的原发癌样本


TCGA.meta  <- data.frame(id = c(TCGA.MSI.sample, TCGA.non.MSI.sample), 
                         type = c(rep('MSI', length(TCGA.MSI.sample)), rep('non.MSI',length(TCGA.non.MSI.sample))))




save(list = c('TCGA.meta', ''))
xCell  <-  xCellAnalysis(TCGA.tpm, 
                         rnaseq = TRUE, # RNA-seq
                         parallel.sz = 10, #线程数 ，window建议改为1
                         parallel.type = "FORK") 

scores      <- xCell
TCGA.scores <- t(scores)


TCGA.scores               <- as.data.frame(TCGA.scores)
TCGA.scores               <- TCGA.scores[TCGA.meta$id,]
identical(rownames(TCGA.scores),TCGA.meta$id)
TCGA.scores               <- cbind(TCGA.scores, TCGA.meta$type)
colnames(TCGA.scores)[68] <- 'type'

TCGA.scores$type          <- factor(TCGA.scores$type, levels = c("non.MSI", "MSI"))

p1 <- ggplot(TCGA.scores, aes(x = type, y = `CD8+ T-cells`, color = type))+
  # boxplot 分组并靠近
  geom_boxplot(position = position_dodge(0.6), outlier.shape = NA, width = 0.5, lwd = 0.6) +
  
  # jitter 散点叠加
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
              size = 1.2)+
  scale_color_manual(values =c("non.MSI" =  '#8FC3DF', 'MSI' = "#E8848E")
  )+
  scale_x_discrete(labels = c('non.MSI'=paste('Non-hypermutated',sep = '\n','(n = 490)'),
                              'MSI' = paste('Hypermutated', sep = '\n','(n = 4)')))+
  stat_compare_means()+
  ylab('CD8+ T cell score')+
  xlab('')+
  ggtitle('')+
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 30),
    axis.title = element_text(size = 30, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 30, color = 'black'),      # 坐标轴刻度加粗
    #panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    panel.grid = element_blank(),
    
    panel.border = element_blank(),  # 移除所有边框
    
    axis.line.x.bottom = element_line(color = "black"), # 添加下方框线
    
    axis.line.y.left = element_line(color = "black"), # 添加左侧框线
    axis.line.x.top = element_blank(), # 隐藏顶部框线
    axis.line.y.right = element_blank(), # 隐藏右侧框线
    
  )





#supp 2d TIDE --------------------------------------------------------------------


### CD8A, CD8B, GZMA, GZMB and PRF1

TCGA.tpm[1:4,1:4]
ctl.exp <- TCGA.tpm[c('CD8A','CD8B','GZMA','GZMB','PRF1'),]
ctl.exp.ave <- colMeans(ctl.exp)
ctl.exp.ave <- as.data.frame(ctl.exp.ave)
ctl.exp.ave <- ctl.exp.ave[TCGA.meta$id,, drop=F]
identical(rownames(ctl.exp.ave), TCGA.meta$id)
ctl.exp.ave$type <- TCGA.meta$type

ctl.exp.ave$type <- factor(ctl.exp.ave$type, levels = c('non.MSI', 'MSI'))

p2 <- ctl.exp.ave %>% ggplot(aes(x = type, y =ctl.exp.ave, color = type))+
  # boxplot 分组并靠近
  geom_boxplot(position = position_dodge(0.6), outlier.shape = NA, width = 0.5, lwd = 0.6) +
  
  # jitter 散点叠加
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6),
              size = 1.2)+
  scale_color_manual(values =c("non.MSI" =  '#8FC3DF', 'MSI' = "#E8848E")
  )+
  scale_x_discrete(labels = c('non.MSI'=paste('Non-hypermutated',sep = '\n','(n = 490)'),
                              'MSI' = paste('Hypermutated', sep = '\n','(n = 4)')))+
  stat_compare_means()+
  ylab('CTL level')+
  xlab('')+
  ggtitle('')+
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 30),
    axis.title = element_text(size = 30, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 30, color = 'black'),      # 坐标轴刻度加粗
    #panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    panel.grid = element_blank(),
    
    panel.border = element_blank(),  # 移除所有边框
    
    axis.line.x.bottom = element_line(color = "black"), # 添加下方框线
    
    axis.line.y.left = element_line(color = "black"), # 添加左侧框线
    axis.line.x.top = element_blank(), # 隐藏顶部框线
    axis.line.y.right = element_blank(), # 隐藏右侧框线
    
  )




pdf('./pcdata/paper/plot/section2/immu.an.pdf', width = 20, height = 8)
p1|p2 
dev.off()







