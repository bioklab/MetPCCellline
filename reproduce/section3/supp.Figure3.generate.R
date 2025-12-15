#Supp 3a bulk bone TC ------------------------------------------------------------
load(file = './pcdata/paper/data/section3/MET500.site.cell.line.correlation.matrix.RData')

# bulk bone
df                    <- MET500.site.cell.line.correlation.matrix %>% as.data.frame()
df <- df[,1,drop = F]
df$name               <- rownames(df)
df$prostate           <- ifelse(rownames(df) %in% pc.cell.line.name,'prostate','others')

df_ordered            <- df[order(df$pc.MET.bone.polyA.sample, decreasing = T),]
df$rank               <- rank(df$pc.MET.bone.polyA.sample) 
df_ranked             <- df[order(df$rank),]

ggplot(df_ranked, aes(x = rank, y = pc.MET.bone.polyA.sample, color = prostate)) +
  geom_point(size = 3) +  # 绘制点，设置大小
  labs(title = "Bone", x = "Rank", y = "Transcriptome similarity") +
  scale_color_manual(values = c("prostate" = "red", "others" = "black")) +  
  theme_bw()+
  geom_point(data = subset(df, name == "PC3_PROSTATE"), color = "red", size = 4) +
  geom_text_repel(data = subset(df, name == "PC3_PROSTATE"| name == 'LNCAPCLONEFGC_PROSTATE' | name == 'VCAP_PROSTATE' | name == 'MDAPCA2B_PROSTATE'), 
                  aes(label = gsub(pattern = '_PROSTATE',x = name, replacement = '')), 
                  nudge_y = 0.02,  # 标签在纵坐标上的偏移量
                  color = "black",
                  size = 4,
                  arrow = arrow(length = unit(0.02, "npc")))+
  ggplot.style+
  theme(legend.position = 'none')


ggsave('./pcdata/paper/plot/section3/Bulk.bone.pdf',width = 12,height=10 )


#supp3b bulk liver TC -----------------------------------------------------------

#bulk liver
df                    <- MET500.site.cell.line.correlation.matrix %>% as.data.frame()
df <- df[,2,drop = F]
df$name               <- rownames(df)
df$prostate           <- ifelse(rownames(df) %in% pc.cell.line.name,'prostate','others')

df_ordered            <- df[order(df$pc.MET.liver.polyA.sample, decreasing = T),]
df$rank               <- rank(df$pc.MET.liver.polyA.sample) 
df_ranked             <- df[order(df$rank),]

ggplot(df_ranked, aes(x = rank, y = pc.MET.liver.polyA.sample, color = prostate)) +
  geom_point(size = 3) +  # 绘制点，设置大小
  labs(title = "Liver", x = "Rank", y = "Transcriptome similarity") +
  scale_color_manual(values = c("prostate" = "red", "others" = "black")) +  
  theme_bw()+
  geom_point(data = subset(df, name == "PC3_PROSTATE"), color = "red", size = 4) +
  geom_text_repel(data = subset(df, name == "PC3_PROSTATE"| name == 'LNCAPCLONEFGC_PROSTATE' | name == 'VCAP_PROSTATE' | name == 'MDAPCA2B_PROSTATE'), 
                  aes(label = gsub(pattern = '_PROSTATE',x = name, replacement = '')), 
                  nudge_y = 0.02,  # 标签在纵坐标上的偏移量
                  color = "black",
                  size = 4,
                  arrow = arrow(length = unit(0.02, "npc")))+
  ggplot.style+
  theme(legend.position = 'none')


ggsave('./pcdata/paper/plot/section3/Bulk.liver.pdf',width = 12,height=10 )







#Supp3c bulk lymph --------------------------------------------------------------


#bulk lymph
df                    <- MET500.site.cell.line.correlation.matrix %>% as.data.frame()
df <- df[,3,drop = F]
df$name               <- rownames(df)
df$prostate           <- ifelse(rownames(df) %in% pc.cell.line.name,'prostate','others')

df_ordered            <- df[order(df$pc.MET.lymph.polyA.sample, decreasing = T),]
df$rank               <- rank(df$pc.MET.lymph.polyA.sample) 
df_ranked             <- df[order(df$rank),]

ggplot(df_ranked, aes(x = rank, y = pc.MET.lymph.polyA.sample, color = prostate)) +
  geom_point(size = 3) +  # 绘制点，设置大小
  labs(title = "Lymph node", x = "Rank", y = "Transcriptome similarity") +
  scale_color_manual(values = c("prostate" = "red", "others" = "black")) +  
  theme_bw()+
  geom_point(data = subset(df, name == "PC3_PROSTATE"), color = "red", size = 4) +
  geom_text_repel(data = subset(df, name == "PC3_PROSTATE"| name == 'LNCAPCLONEFGC_PROSTATE' | name == 'VCAP_PROSTATE' | name == 'MDAPCA2B_PROSTATE'), 
                  aes(label = gsub(pattern = '_PROSTATE',x = name, replacement = '')), 
                  nudge_y = 0.02,  # 标签在纵坐标上的偏移量
                  color = "black",
                  size = 4,
                  arrow = arrow(length = unit(0.02, "npc")))+
  
  theme(legend.position = 'none')


ggsave('./pcdata/paper/plot/section3/Lymph.node.pdf',width = 12,height=10)
#Supp3d-e scRNA-seq site specific ---------------------------------------------------------------

load(file = './pcdata/paper/data/section2/TC.res.RData')

ln.res$median.cor %>% View()
ln.median.cor                      <- ln.res$pc.median.cor %>% as.data.frame()
bone.res$median.cor %>% View()


#Supp 3d scRNA bone
df                    <- bone.res$median.cor %>% as.data.frame()
df$name               <- rownames(df)
df$prostate           <- ifelse(rownames(df)%in% pc.cell.line.name,'prostate','others')
colnames(df)[1] <- 'correlation'
df_ordered            <- df[order(df$correlation, decreasing = T),]
df$rank               <- rank(df$correlation) 
df_ranked             <- df[order(df$rank),]

ggplot(df_ranked, aes(x = rank, y = correlation, color = prostate)) +
  geom_point(size = 3) +  # 绘制点，设置大小
  labs(title = "scRNA bone", x = "Rank", y = "Transcriptome similarity") +
  scale_color_manual(values = c("prostate" = "red", "others" = "black")) +  
  geom_point(data = subset(df, name == "PC3_PROSTATE"), color = "red", size = 4) +
  geom_text_repel(data = subset(df, name == "PC3_PROSTATE"| name == 'LNCAPCLONEFGC_PROSTATE' | name == 'VCAP_PROSTATE' | name == 'MDAPCA2B_PROSTATE'|name=='22RV1_PROSTATE'), 
                  aes(label = gsub(pattern = '_PROSTATE',x = name, replacement = '')), 
                  nudge_y = 0.02,  # 标签在纵坐标上的偏移量
                  color = "black",
                  size = 4,
                  arrow = arrow(length = unit(0.02, "npc")))+
  ggtitle('Bone')+
  
  theme(legend.position = 'none')
ggsave(sprintf("%s/%s",figure.folder,'scRNA.bone.pdf'),width = 12,height=10)

#Supp 3e scRNA lymph
df                    <- ln.res$median.cor %>% as.data.frame()
df$name               <- rownames(df)
df$prostate           <- ifelse(rownames(df)%in% pc.cell.line.name,'prostate','others')
colnames(df)[1] <- 'correlation'
df_ordered            <- df[order(df$correlation, decreasing = T),]
df$rank               <- rank(df$correlation) 
df_ranked             <- df[order(df$rank),]

ggplot(df_ranked, aes(x = rank, y = correlation, color = prostate)) +
  geom_point(size = 3) +  # 绘制点，设置大小
  labs(title = "scRNA lymph", x = "Rank", y = "Transcriptome similarity") +
  scale_color_manual(values = c("prostate" = "red", "others" = "black")) +  
  geom_point(data = subset(df, name == "PC3_PROSTATE"), color = "red", size = 4) +
  geom_text_repel(data = subset(df, name == "PC3_PROSTATE"| name == 'LNCAPCLONEFGC_PROSTATE' | name == 'VCAP_PROSTATE' | name == 'MDAPCA2B_PROSTATE'), 
                  aes(label = gsub(pattern = '_PROSTATE',x = name, replacement = '')), 
                  nudge_y = 0.02,  # 标签在纵坐标上的偏移量
                  color = "black",
                  size = 4,
                  arrow = arrow(length = unit(0.02, "npc")))+
  ggtitle('Lymph node')+
  
  
  theme(legend.position = 'none')
ggsave(sprintf("%s/%s",figure.folder,'scRNA.lymph.pdf'),width = 12,height=10)







#SuppFig3f pubmed citation plot ----------------------------------------------------
Supplementary_Table_2 <- read_csv("./pcdata/paper/data/section3/pubmed.citation.csv")
colnames(Supplementary_Table_2) <- c('cell.line.name', 'pubmed.citation.number')


Supplementary_Table_2$pubmed.citation.number <- as.numeric(Supplementary_Table_2$pubmed.citation.number)
Supplementary_Table_2$percent <- Supplementary_Table_2$pubmed.citation.number/sum(Supplementary_Table_2$pubmed.citation.number)*100


# ggplot(Supplementary_Table_2, aes(x = reorder(cell.line.name, percent, decreasing=T), y = percent))+geom_point(size=6)+xlab('')+
#   ylab('Percent')+
#   theme_bw()+
#   RotatedAxis()
ggplot(Supplementary_Table_2, aes(x = reorder(cell.line.name, percent, decreasing = TRUE), y = percent)) +
  geom_col(width = 0.7, fill = "steelblue") +  # 使用geom_col生成柱状图，设定柱宽和颜色
  xlab('') +
  ylab('Percentage(%)') +
  theme_bw() +
  RotatedAxis() +  # 如果这是自定义函数用于旋转x轴标签，可以保留
  theme(
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 16, face = "bold"),
    panel.border = element_rect(color = "black", linewidth = 1.2, fill = NA)
  )

ggsave('./pcdata/paper/plot/section3/pubmed.citation.pdf', width = 8, height = 6)

