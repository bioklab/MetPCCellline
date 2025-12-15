#Supp Fig1a-c dm mutation site-specific -----------------------------------------------
load(file = './output/compare.somatic.mutation.R.output/compare.somatic.mutation.RData')

load(file = './pcdata/paper/data/section1/somatic.site.dm.gene.RData')
bone.dm.rs$sig <- ifelse(bone.dm.rs$fdr < 0.05,'yes','no'
)
de_genes <- MET.dm.gene
ggplot(bone.dm.rs, aes(x = delta, y = -log10(fdr),color=sig))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  #geom_vline(xintercept = c(-0.3,0.3), linetype = "dashed", color = "black")+
  geom_point()+
  scale_color_manual(values=c("#BBBBBB","#2f5688"))+
  scale_size_continuous(range = c(0,1))+
  theme(legend.title = element_blank() )+
  geom_label_repel(data = subset(bone.dm.rs, gene %in% de_genes),
                   max.overlaps = Inf,
                   aes(label = gene),
                   size = 9, 
                   color = 'black',fontface = "bold.italic",
                   segment.color = "black",   
                   segment.linetype = "solid", 
                   segment.size = 0.5 ,
                   box.padding = 3,        
                   point.padding = 0.2 ) +
  xlab("Delta")+
  ylab("-log10 p-adj")+
  ggtitle('Bone')+
  ggplot.style+
  theme(legend.position = 'none')

ggsave('./pcdata/paper/plot/section1/dm.somatic.mutation.bone.pdf', width = 15, height = 15)




liver.dm.rs$sig <- ifelse(liver.dm.rs$fdr < 0.05,'yes','no'
)
de_genes <- c('AR','TP53','MUC4','FAM189B')
ggplot(liver.dm.rs, aes(x = delta, y = -log10(fdr),color=sig))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  #geom_vline(xintercept = c(-0.3,0.3), linetype = "dashed", color = "black")+
  geom_point()+
  scale_color_manual(values=c("#BBBBBB","#2f5688"))+
  scale_size_continuous(range = c(0,1))+
  theme(legend.title = element_blank() )+
  geom_label_repel(data = subset(liver.dm.rs, gene %in% de_genes),
                   max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                   aes(label = gene),
                   size = 9, 
                   color = 'black',fontface = "bold.italic",
                   segment.color = "black",   
                   segment.linetype = "solid", 
                   segment.size = 0.5 ,
                   box.padding = 3,        
                   point.padding = 0.2 ) +
  xlab("Delta")+
  ylab("-log10 p-adj")+
  ggtitle('Liver')+
  ggplot.style+
  theme(legend.position = 'none')

ggsave('./pcdata/paper/plot/section1/dm.somatic.mutation.liver.pdf', width = 15, height = 15)




ln.dm.rs$sig <- ifelse(ln.dm.rs$fdr < 0.05,'yes','no'
)
de_genes <- MET.dm.gene
ggplot(ln.dm.rs, aes(x = delta, y = -log10(fdr),color=sig))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  #geom_vline(xintercept = c(-0.3,0.3), linetype = "dashed", color = "black")+
  geom_point()+
  scale_color_manual(values=c("#BBBBBB","#2f5688"))+
  scale_size_continuous(range = c(0,1))+
  theme(legend.title = element_blank() )+
  geom_label_repel(data = subset(ln.dm.rs, gene %in% de_genes),
                   max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                   aes(label = gene),
                   size = 9, 
                   color = 'black',fontface = "bold.italic",
                   segment.color = "black",   
                   segment.linetype = "solid", 
                   segment.size = 0.5 ,
                   box.padding = 3,        
                   point.padding = 0.2 ) +
  xlab("Delta")+
  ylab("-log10 p-adj")+
  ggtitle('Lymph node')+
  ggplot.style+
  theme(legend.position = 'none')

ggsave('./pcdata/paper/plot/section1/dm.somatic.mutation.ln.pdf', width = 15, height = 15)



#Supp1d hotspot barplot ---------------------------------------------------------

load(file = './pcdata/paper/data/section1/supplementary.table1.RData')
supplementary.table1 %>% View()

supplementary.table1 <- na.omit(supplementary.table1)
supplementary.table1$HGNC.symbol %>% table() %>% sort(decreasing = T) %>% length()
View(supplementary.table1)
df <- supplementary.table1 %>%
  group_by(HGNC.symbol, mutation.code) %>%
  tally(name = "count") %>%
  filter(count >= 3) %>%
  arrange(desc(count))
gene.count <- df$HGNC.symbol %>% table() %>% sort(decreasing = T) %>% as.data.frame()
colnames(gene.count)[1] <- 'gene'
ggplot(data = gene.count )+geom_point(aes(x = gene, y = Freq))
df <- gene.count
df$color <- ifelse(df$Freq >= 3, "High", "Low")

ggplot(df, aes(y = reorder(gene, Freq), x = Freq, fill = color)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("High" = "firebrick", "Low" = "grey70")) +
  coord_flip() +
  labs(x = "Frequency", y = "", title = "") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),  # 横坐标标签倾斜
    axis.text.y = element_text(size = 14, face = "bold"),
    panel.border = element_rect(),
    axis.line.x.bottom = element_line(color = "black"), # 添加下方框线
    
    axis.line.y.left = element_line(color = "black"), # 添加左侧框线
    axis.line.x.top = element_blank(), # 隐藏顶部框线
    axis.line.y.right = element_blank(), # 隐藏右侧框线
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    axis.ticks = element_line(size = 1.6)        # 增加刻度线粗细
  )


ggsave('./pcdata/paper/plot/section1/hotspot.count.pdf', width = 12, height = 10)










# Supp Figure1 e-g --------------------------------------------------------

load(file = './pcdata/paper/data/section1/dm.cnv.three.site.RData')
load(file = '/home/liuxueying/met_pc_cell_line/pcdata/paper/data/section1/su2cid.RData')

# volcano plot ------------------------------------------------------------
rs.df.LN$sig <- ifelse(rs.df.LN$p.adj<0.05 & rs.df.LN$delta > 0.3, 'up',
                       ifelse(rs.df.LN$p.adj <0.05 & rs.df.LN$delta < -0.3, 'dn', 'no'))
rs.df.LN$log10.p.adj <- -log10(rs.df.LN$p.adj)
de_genes <- c('AR')

##Lymph node
p1 <- ggplot(rs.df.LN, aes(x = delta, y = -log10(p.adj),color=sig))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-0.3,0.3), linetype = "dashed", color = "black")+
  geom_point()+
  scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
  scale_size_continuous(range = c(0,1))+
  theme(legend.title = element_blank() )+
  geom_label_repel(data = subset(rs.df.LN, gene %in% de_genes),
                   max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                   aes(label = gene),
                   size = 5, 
                   color = 'black',fontface = "bold.italic",
                   segment.color = "black",   
                   segment.linetype = "solid", 
                   segment.size = 0.5 ,
                   box.padding = 3,        
                   point.padding = 0.2 ) +
  xlab("delta")+
  ylab("-log10 p-adj")+
  theme_bw()+
  ggtitle('Lymph node')+
  geom_hline(yintercept = 50)+
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),     # 坐标轴标题加粗
    axis.text = element_text(size = 14, face = "bold"),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA) # 加粗边框
  )

## bone
rs.df.bone$sig <- ifelse(rs.df.bone$p.adj<0.05 & rs.df.bone$delta > 0.3, 'up',
                         ifelse(rs.df.bone$p.adj <0.05 & rs.df.bone$delta < -0.3, 'dn', 'no'))



p2 <- ggplot(rs.df.bone, aes(x = delta, y = -log10(p.adj),color=sig))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-0.3,0.3), linetype = "dashed", color = "black")+
  geom_point()+
  scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
  scale_size_continuous(range = c(0,1))+
  theme(legend.title = element_blank() )+
  geom_label_repel(data = subset(rs.df.bone, gene %in% de_genes),
                   max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                   aes(label = gene),
                   size = 5, 
                   color = 'black',fontface = "bold.italic",
                   segment.color = "black",   
                   segment.linetype = "solid", 
                   segment.size = 0.5 ,
                   box.padding = 3,        
                   point.padding = 0.2 ) +
  xlab("delta")+
  ylab("-log10 p-adj")+
  theme_bw()+
  ggtitle('Bone')+
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),     # 坐标轴标题加粗
    axis.text = element_text(size = 14, face = "bold"),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA) # 加粗边框
  )


## liver
rs.df.liver$sig <- ifelse(rs.df.liver$p.adj<0.05 & rs.df.liver$delta > 0.3, 'up',
                          ifelse(rs.df.liver$p.adj <0.05 & rs.df.liver$delta < -0.3, 'dn', 'no'))



p3 <- ggplot(rs.df.liver, aes(x = delta, y = -log10(p.adj),color=sig))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-0.3,0.3), linetype = "dashed", color = "black")+
  geom_point()+
  geom_hline(yintercept = 50)+
  geom_vline(xintercept = 0.33)+
  scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
  scale_size_continuous(range = c(0,1))+
  theme(legend.title = element_blank() )+
  geom_label_repel(data = subset(rs.df.liver, gene %in% de_genes),
                   max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                   aes(label = gene),
                   size = 5, 
                   color = 'black',fontface = "bold.italic",
                   segment.color = "black",   
                   segment.linetype = "solid", 
                   segment.size = 0.5 ,
                   box.padding = 3,        
                   point.padding = 0.2 ) +
  geom_label_repel(
    data = subset(rs.df.liver, delta > 0.3 & -log10(p.adj) > 50),
    aes(label = gene),
    size = 5,
    color = 'black',
    
    segment.size = 0.5,
    box.padding = 3,
    point.padding = 0.2,
    max.overlaps = Inf
  )+
  xlab("delta")+
  ylab("-log10 p-adj")+
  theme_bw()+
  ggtitle('Liver')+
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),     # 坐标轴标题加粗
    axis.text = element_text(size = 14, face = "bold"),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA) # 加粗边框
  )
ggsave('./pcdata/paper/plot/section1/aa.pdf', width = 29,height = 15)
pdf('./pcdata/paper/plot/section1/dm.cnv.supp.pdf', width = 15, height = 6)
p1|p2|p3             
dev.off()


