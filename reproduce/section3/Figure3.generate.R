
# Figure.3a ---------------------------------------------------------------

#pc.cell.line.name
flag                         <- grepl(colnames(CCLE.log2.tpm.matrix), pattern ='PROSTATE')
pc.cell.line.log2.tpm.matrix <- CCLE.log2.tpm.matrix[,flag]
flag                         <- colnames(pc.cell.line.log2.tpm.matrix) != 'PRECLH_PROSTATE' ## remove no cancer 
pc.cell.line.log2.tpm.matrix <- pc.cell.line.log2.tpm.matrix[,flag]
pc.cell.line.name            <- colnames(pc.cell.line.log2.tpm.matrix) 

# TC rs function 
source('./code/BioKLab.util.R')
source('./pcdata/my_function.R')
load("~/met_pc_cell_line/pcdata/CCLE.rna.seq.marker.gene.1000.symbol.RData")
load ('./processed.data/CCLE.transcriptome.RData')
CCLE.log2.tpm.matrix.symbol <- change_rownames_ensemble_to_symbol(CCLE.log2.tpm.matrix)
TC.rs.list <- function(seob){
  TC.rs             <- pick.out.cell.line(expr.of.samples = seob@assays$RNA$data %>% as.matrix(), expr.of.cell.lines = CCLE.log2.tpm.matrix.symbol,marker.gene = CCLE.rna.seq.marker.gene.1000.symbol)
  median.cor        <- TC.rs$cell.line.median.cor %>% as.data.frame()
  pc.cell.line.name <- grep(pattern = '_PROSTATE', x = rownames(median.cor), value = T)
  pc.median.cor     <- median.cor[pc.cell.line.name,, drop=F]
  TC.rs.list        <-  list(median.cor = median.cor,
                             pc.median.cor = pc.median.cor,
                             correlation.matrix = TC.rs[["correlation.matrix"]])
}
# bulk Adenocarcinoma TC analysis 
#MET500 prostate adenocarcinoma
load('./processed.data/BioklabData/liuke/KeData/MET500/MET500.RData')
flag                                     <- grepl(x=MET500.sample.meta$cancer.type,pattern='Prostate Adenocarcinoma') 
MET500.sample                            <- MET500.sample.meta$Run[flag]
tmp                                      <- MET500.sample.meta[MET500.sample,]

pc.MET.bone.polyA.sample  <- tmp$Run[tmp$biopsy.site == 'BONE'       & tmp$LibrarySelection == 'PolyA']
pc.MET.liver.polyA.sample <- tmp$Run[tmp$biopsy.site == 'LIVER'      & tmp$LibrarySelection == 'PolyA']
pc.MET.lymph.polyA.sample <- tmp$Run[tmp$biopsy.site == 'LYMPH_NODE' & tmp$LibrarySelection == 'PolyA']

pc.MET.bone.hs.sample  <- tmp$Run[tmp$biopsy.site == 'BONE'       & tmp$LibrarySelection == 'Hybrid Selection']
pc.MET.liver.hs.sample <- tmp$Run[tmp$biopsy.site == 'LIVER'      & tmp$LibrarySelection == 'Hybrid Selection']
pc.MET.lymph.hs.sample <- tmp$Run[tmp$biopsy.site == 'LYMPH_NODE' & tmp$LibrarySelection == 'Hybrid Selection']


pc.MET.sample.list <- list(pc.MET.bone.polyA.sample= pc.MET.bone.polyA.sample,
                           pc.MET.liver.polyA.sample = pc.MET.liver.polyA.sample,
                           pc.MET.lymph.polyA.sample = pc.MET.lymph.polyA.sample,
                           pc.MET.bone.hs.sample=pc.MET.bone.hs.sample,
                           pc.MET.liver.hs.sample = pc.MET.liver.hs.sample,
                           pc.MET.lymph.hs.sample = pc.MET.lymph.hs.sample
)



#  Compute PCa cell line correlation with all MET500 samples
source('./code/BioKLab.util.R')
MET500.TC.rs   <- pick.out.cell.line(expr.of.samples = MET500.log2.tpm.matrix[,pc.MET.sample.list %>% unlist],expr.of.cell.lines = CCLE.log2.tpm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
MET500.TC.rs$cell.line.median.cor %>% as.data.frame() %>% View()
df                    <- MET500.TC.rs$cell.line.median.cor %>% as.data.frame()
df$name               <- rownames(df)
df$prostate           <- ifelse(rownames(df) %in% pc.cell.line.name,'prostate','others')
colnames(df)[1]       <- 'correlation'
df_ordered            <- df[order(df$correlation, decreasing = T),]
df$rank               <- rank(df$correlation) 
df_ranked             <- df[order(df$rank),]

p1 <- ggplot(df_ranked, aes(x = rank, y = correlation, color = prostate)) +
  geom_point(size = 3) +  # 绘制点，设置大小
  labs(title = "", x = "Rank", y = "Transcriptome similarity") +
  scale_color_manual(values = c("prostate" = "red", "others" = "black")) +  
  geom_point(data = subset(df, name == "PC3_PROSTATE"), color = "red", size = 4) +
  geom_text_repel(data = subset(df, name == "PC3_PROSTATE"| name == 'LNCAPCLONEFGC_PROSTATE' | name == 'VCAP_PROSTATE' | name == 'MDAPCA2B_PROSTATE'), 
                  aes(label = gsub(x= name, pattern = '_PROSTATE','')), 
                  nudge_y = 0.02,  # 标签在纵坐标上的偏移量
                  color = "black",
                  size = 4,
                  arrow = arrow(length = unit(0.02, "npc")))+
  ggplot.style+
  theme(legend.position = 'none')

ggsave('./pcdata/paper/plot/section2/MET500.adenocarcinoma.pdf', plot = p1,width = 15, height = 10)


#Fig.3b bulk RNA site-specific ---------------------------------------------------
load('./processed.data/BioklabData/liuke/KeData/MET500/MET500.RData')
flag                                     <- grepl(x=MET500.sample.meta$cancer.type,pattern='Prostate Adenocarcinoma') 
MET500.sample                            <- MET500.sample.meta$Run[flag]
tmp                                      <- MET500.sample.meta[MET500.sample,]

pc.MET.bone.polyA.sample  <- tmp$Run[tmp$biopsy.site == 'BONE'       & tmp$LibrarySelection == 'PolyA']
pc.MET.liver.polyA.sample <- tmp$Run[tmp$biopsy.site == 'LIVER'      & tmp$LibrarySelection == 'PolyA']
pc.MET.lymph.polyA.sample <- tmp$Run[tmp$biopsy.site == 'LYMPH_NODE' & tmp$LibrarySelection == 'PolyA']

pc.MET.bone.hs.sample  <- tmp$Run[tmp$biopsy.site == 'BONE'       & tmp$LibrarySelection == 'Hybrid Selection']
pc.MET.liver.hs.sample <- tmp$Run[tmp$biopsy.site == 'LIVER'      & tmp$LibrarySelection == 'Hybrid Selection']
pc.MET.lymph.hs.sample <- tmp$Run[tmp$biopsy.site == 'LYMPH_NODE' & tmp$LibrarySelection == 'Hybrid Selection']


pc.MET.sample.list <- list(pc.MET.bone.polyA.sample= pc.MET.bone.polyA.sample,
                           pc.MET.liver.polyA.sample = pc.MET.liver.polyA.sample,
                           pc.MET.lymph.polyA.sample = pc.MET.lymph.polyA.sample,
                           pc.MET.bone.hs.sample=pc.MET.bone.hs.sample,
                           pc.MET.liver.hs.sample = pc.MET.liver.hs.sample,
                           pc.MET.lymph.hs.sample = pc.MET.lymph.hs.sample
)


MET500.site.TC.rs.list <- foreach (item = pc.MET.sample.list) %do% {
  rs   <- pick.out.cell.line(expr.of.samples = MET500.log2.tpm.matrix[,item],expr.of.cell.lines = CCLE.log2.tpm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
  rs
}
names(MET500.site.TC.rs.list)      <-  names(pc.MET.sample.list)
all.cell.line                      <-  colnames(CCLE.log2.tpm.matrix)
pc.cell.line                       <-  colnames(CCLE.log2.tpm.matrix)[colnames(CCLE.log2.tpm.matrix) %>% grepl(pattern='PROSTATE')]

MET500.site.cell.line.correlation.matrix <- foreach(item = MET500.site.TC.rs.list,.combine='cbind')  %do% {
  item$cell.line.median.cor[all.cell.line]
}
colnames(MET500.site.cell.line.correlation.matrix)     <- names(MET500.site.TC.rs.list)

MET500.site.cell.line.pc.correlation.matrix <- MET500.site.cell.line.correlation.matrix[pc.cell.line.name,]
flag                                                   <- rownames(MET500.site.cell.line.pc.correlation.matrix) != 'PRECLH_PROSTATE'
MET500.site.cell.line.pc.correlation.matrix            <- MET500.site.cell.line.pc.correlation.matrix[flag,]

cor1 <- cor(MET500.site.cell.line.pc.correlation.matrix[,1],MET500.site.cell.line.pc.correlation.matrix[,2], method = 'spearman')
df <- MET500.site.cell.line.pc.correlation.matrix
site.specific.cor <- MET500.site.cell.line.pc.correlation.matrix[,c(1,2,3)]
colnames(site.specific.cor) <- c('Bone','Liver','Lymph')
mat <- site.specific.cor

source('./pcdata/paper/code/chart.correlation.R')
pdf(file  = './pcdata/paper/plot/section3/bulk.site.specific.pdf',width=20,height=20)
chart.Correlation(mat, histogram=FALSE, pch=19,method = 'spearman',cex.cor.tuning = 0.7,cor.range = c(-1,1))
dev.off()

my_palette <- colorRampPalette(c("blue", "red"))(n = 101)
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}
pdf(file  = './pcdata/paper/plot/section3/color.bar.pdf',width=20,height=20)
color.bar(my_palette,min = -1,max=1) # draw the color bar
dev.off()




#Fig3c scRNA adenocacinoma -----------------------------------------------------------
source('./code/BioKLab.util.R')

TC.rs.list <- function(seob){
  TC.rs             <- pick.out.cell.line(expr.of.samples = seob@assays$RNA@data %>% as.matrix(), expr.of.cell.lines = CCLE.log2.tpm.matrix.symbol,marker.gene = CCLE.rna.seq.marker.gene.1000.symbol)
  median.cor        <- TC.rs$cell.line.median.cor %>% as.data.frame()
  pc.cell.line.name <- grep(pattern = '_PROSTATE', x = rownames(median.cor), value = T)
  pc.median.cor     <- median.cor[pc.cell.line.name,, drop=F]
  TC.rs.list        <-  list(median.cor = median.cor,
                             pc.median.cor = pc.median.cor,
                             correlation.matrix = TC.rs[["correlation.matrix"]])
}
#14 mCRPC
load('./pcdata/SCP1244_14mCRPC/scp1224.seob.RData')

scp1224.malignant.seob <- subset(scp1224.seob, cell_type == 'prostate cancer cell')


adeno.seob <- subset(scp1224.malignant.seob, biosample_id != '09171135') # remove small cell cacinoma

adeno.res <- TC.rs.list(adeno.seob)
adeno.res$median.cor %>% View()

# plot --------------------------------------------------------------------

df                    <- adeno.res$median.cor %>% as.data.frame()
df$name               <- rownames(df)
df$prostate           <- ifelse(rownames(df) %in% pc.cell.line.name, 'prostate', 'others')
colnames(df)[1]       <- 'correlation'
df_ordered            <- df[order(df$correlation, decreasing = T),]
df$rank               <- rank(df$correlation) 
df_ranked             <- df[order(df$rank),]

df_highlight <- subset(df, name %in% c("PC3_PROSTATE", "LNCAPCLONEFGC_PROSTATE", "VCAP_PROSTATE", "MDAPCA2B_PROSTATE"))
ggplot(df_ranked, aes(x = rank, y = correlation, color = prostate)) +
  geom_point(size = 2) +  # 绘制所有点
  labs(title = "Adenocarcinoma scRNA-seq", x = "Rank", y = "Transcriptome similarity") +
  scale_color_manual(values = c("prostate" = "red", "others" = "black")) +  
  geom_point(data = subset(df, prostate == "prostate"), color = "red", size = 3)+
  ggplot.style +
  theme(legend.position = 'none')+
  
  # 添加箭头，从标签尾部指向数据点
  geom_segment(data = df_highlight, 
               aes(x = rank - 50, xend = rank, y = correlation, yend = correlation),  # x 起点是标签，xend 是数据点
               arrow = arrow(length = unit(0.01, "npc")),  # 添加箭头
               color = "black") +
  
  # 添加标签
  geom_text(data = df_highlight, 
            aes(x = rank - 50, y = correlation, label = gsub("_PROSTATE", "", name)),  # 标签位于箭头起点
            hjust = 1,  # 右对齐标签
            color = "black", 
            size = 4)

figure.folder <- "./pcdata/paper/plot/section2/"
ggsave(sprintf("%s/%s",figure.folder,'scRNA-seq.Adenocarcinoma.pdf'),width = 20,height=15 )

#scRNA-seq site-specific correlation ---------------------------------------------------------- 
# bone
adeno.seob$organ__ontology_label %>% table()
bone.seob <- subset(adeno.seob, organ__ontology_label == 'bone tissue')
bone.res <- TC.rs.list(bone.seob)
bone.res$median.cor %>% View()
# lymph node
ln.seob   <- subset(adeno.seob, organ__ontology_label == 'lymph node')
ln.res    <- TC.rs.list(ln.seob)
ln.res$median.cor %>% View()



# bone ln scRNA spearman correlation
bone.median.cor                    <- bone.res$pc.median.cor %>% as.data.frame()
colnames(bone.median.cor)          <- 'bone.cor.median'
flag                               <- rownames(bone.median.cor) != 'PRECLH_PROSTATE'
bone.pc.median.cor                 <- bone.median.cor[flag,,drop = F ]

ln.median.cor                      <- ln.res$pc.median.cor %>% as.data.frame()
colnames(ln.median.cor)            <- 'ln.cor.median'
flag                <- rownames(ln.median.cor) !='PRECLH_PROSTATE'
ln.pc.median.cor    <- ln.median.cor[flag,,drop = F]
ln.pc.median.cor    <- ln.pc.median.cor[rownames(bone.pc.median.cor),,drop = F]
c.median.cor        <- cbind(ln.pc.median.cor, bone.pc.median.cor)
bone.ln.cor         <- cor(c.median.cor[,1],c.median.cor[,2], method = 'spearman')
df                  <- c.median.cor

ggplot(df, aes(x= ln.cor.median, y=bone.cor.median))+
  geom_point(size=8)+
  theme_bw()+
  xlab('Lymph node')+
  ylab('Bone')+
  ggtitle('')+
  xlim(c(0,0.3))+ylim(c(0,0.3))+geom_abline(slope = 1,intercept = 0)+
  annotate(
    'text',
    x = 0.2,            # x 轴位置（在两个柱中间）
    y = 0.15,           # y 轴位置（高于图顶端）
    label = paste0("Spearman-rank correlation = ", round(bone.ln.cor, digits = 2)), # p 值标签
    size = 5,           # 字体大小
    color = "black",      # 字体颜色
    hjust = 0           # 水平对齐方式
  )+theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )

