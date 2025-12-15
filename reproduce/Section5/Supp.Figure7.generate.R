# mspc pc3 up 2126 富集分析 --------------------------------------------------
mspc.pc3.up <- read.delim("./pcdata/paper/data/section4/MSPC_PC3_up_hallmark.txt")


mspc.pc3.up <- mspc.pc3.up %>% 
  separate(col = 'Overlap',into = c('count','gene_set_count'),sep = '/',remove = F)

mspc.pc3.up$count <- as.numeric(mspc.pc3.up$count)
mspc.pc3.up <- mspc.pc3.up %>% mutate(gene_ratio=mspc.pc3.up$count/2126)


mspc.pc3.up <- mspc.pc3.up[order(-mspc.pc3.up$count), ]
mspc.pc3.up <- mspc.pc3.up %>% filter(Adjusted.P.value < 0.05)
mspc.pc3.up <- mspc.pc3.up[order(mspc.pc3.up$P.value), ] %>% head(n=10)

ggplot(mspc.pc3.up, aes(x = gene_ratio, y = reorder(Term, count,decreasing = F))) +
  geom_point(aes(size = count, fill = P.value), shape = 21) +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(
    
    size = "Count",
    x = "Gene Ratio",
    y = "",
    title = "PC3-vs-MSPC down-regulated "
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 38))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color="black", linetype="solid", size = 1),
        axis.text.x = element_text( size = 20,color = 'black'),
        axis.text.y = element_text(size = 20,color = 'black'),
        axis.title.x = element_text(size = 20,color = 'black'),
        legend.text = element_text(color = "black",size = 15),
        legend.title = element_text(color = "black",size = 15),
        axis.ticks.x = element_line(size=1),
        axis.ticks.y = element_line(size=1),
        plot.title = element_text(color = "black",face = "bold", size = 17,hjust = 0.5))
ggsave('./pcdata/paper/plot/section4/PC3-vs-MSPC down-regulated.pdf', width = 10, height = 8)

# mspc pc3 dn 967 富集分析 --------------------------------------------------
mspc.pc3.dn <- read.delim("./pcdata/paper/data/section4/MSPC_PC3_dn_hallmark.txt")


mspc.pc3.dn <- mspc.pc3.dn %>% 
  separate(col = 'Overlap',into = c('count','gene_set_count'),sep = '/',remove = F)

mspc.pc3.dn$count <- as.numeric(mspc.pc3.dn$count)
mspc.pc3.dn <- mspc.pc3.dn %>% mutate(gene_ratio=mspc.pc3.dn$count/967)


mspc.pc3.dn <- mspc.pc3.dn[order(-mspc.pc3.dn$count), ]
mspc.pc3.dn <- mspc.pc3.dn %>% filter(Adjusted.P.value < 0.05)
mspc.pc3.dn <- mspc.pc3.dn[order(mspc.pc3.dn$P.value), ] %>% head(n=10)

ggplot(mspc.pc3.dn, aes(x = gene_ratio, y = reorder(Term, count,decreasing = F))) +
  geom_point(aes(size = count, fill = P.value), shape = 21) +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(
    
    size = "Count",
    x = "Gene Ratio",
    y = "",
    title = "PC3-vs-MSPC up-regulated"
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 38))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color="black", linetype="solid", size = 1),
        axis.text.x = element_text( size = 20,color = 'black'),
        axis.text.y = element_text(size = 20,color = 'black'),
        axis.title.x = element_text(size = 20,color = 'black'),
        legend.text = element_text(color = "black",size = 15),
        legend.title = element_text(color = "black",size = 15),
        axis.ticks.x = element_line(size=1),
        axis.ticks.y = element_line(size=1),
        plot.title = element_text(color = "black",face = "bold", size = 17,hjust = 0.5))
ggsave('./pcdata/paper/plot/section4/PC3-vs-MSPC up-regulated.pdf', width = 10, height = 8)
# nepc.ncih660.up 1852     富集分析 --------------------------------------------------
nepc.ncih660.up <- read.delim("./pcdata/paper/data/section4/NEPC_NCIH660_up_hallmark.txt")


nepc.ncih660.up <- nepc.ncih660.up %>% 
  separate(col = 'Overlap',into = c('count','gene_set_count'),sep = '/',remove = F)

nepc.ncih660.up$count <- as.numeric(nepc.ncih660.up$count)
nepc.ncih660.up <- nepc.ncih660.up %>% mutate(gene_ratio=nepc.ncih660.up$count/1852)


nepc.ncih660.up <- nepc.ncih660.up[order(-nepc.ncih660.up$count), ]
nepc.ncih660.up <- nepc.ncih660.up %>% filter(Adjusted.P.value < 0.05)
nepc.ncih660.up <- nepc.ncih660.up[order(nepc.ncih660.up$P.value), ] %>% head(n=10)

ggplot(nepc.ncih660.up, aes(x = gene_ratio, y = reorder(Term, count,decreasing = F))) +
  geom_point(aes(size = count, fill = P.value), shape = 21) +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(
    
    size = "Count",
    x = "Gene Ratio",
    y = "",
    title = "NCI-H660-vs-NEPC down-regulated"
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 38))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color="black", linetype="solid", size = 1),
        axis.text.x = element_text( size = 20,color = 'black'),
        axis.text.y = element_text(size = 20,color = 'black'),
        axis.title.x = element_text(size = 20,color = 'black'),
        legend.text = element_text(color = "black",size = 15),
        legend.title = element_text(color = "black",size = 15),
        axis.ticks.x = element_line(size=1),
        axis.ticks.y = element_line(size=1),
        plot.title = element_text(color = "black",face = "bold", size = 17,hjust = 0.5))

ggsave('./pcdata/paper/plot/section4/NCI-H660-vs-NEPC down-regulated.pdf', width = 10, height = 8)
# nepc.ncih660.dn 2650 富集分析 --------------------------------------------------
nepc.ncih660.dn <- read.delim("./pcdata/paper/data/section4/NEPC_NCIH660_dn_hallmark.txt")


nepc.ncih660.dn <- nepc.ncih660.dn %>% 
  separate(col = 'Overlap',into = c('count','gene_set_count'),sep = '/',remove = F)

nepc.ncih660.dn$count <- as.numeric(nepc.ncih660.dn$count)
nepc.ncih660.dn <- nepc.ncih660.dn %>% mutate(gene_ratio=nepc.ncih660.dn$count/2650)


nepc.ncih660.dn <- nepc.ncih660.dn[order(-nepc.ncih660.dn$count), ]
nepc.ncih660.dn <- nepc.ncih660.dn %>% filter(Adjusted.P.value < 0.05)
nepc.ncih660.dn <- nepc.ncih660.dn[order(nepc.ncih660.dn$P.value), ] %>% head(n=10)

ggplot(nepc.ncih660.dn, aes(x = gene_ratio, y = reorder(Term, count,decreasing = F))) +
  geom_point(aes(size = count, fill = P.value), shape = 21) +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(
    
    size = "Count",
    x = "Gene Ratio",
    y = "",
    title = "NCI-H660-vs-NEPC up-regulated"
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 38))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color="black", linetype="solid", size = 1),
        axis.text.x = element_text( size = 20,color = 'black'),
        axis.text.y = element_text(size = 20,color = 'black'),
        axis.title.x = element_text(size = 20,color = 'black'),
        legend.text = element_text(color = "black",size = 15),
        legend.title = element_text(color = "black",size = 15),
        axis.ticks.x = element_line(size=1),
        axis.ticks.y = element_line(size=1),
        plot.title = element_text(color = "black",face = "bold", size = 17,hjust = 0.5))

ggsave('./pcdata/paper/plot/section4/NCI-H660-vs-NEPC up-regulated.pdf', width = 10, height = 8)

#arpc.vcap.up 3855 富集分析 --------------------------------------------------
arpc.vcap.up <- read.delim("./pcdata/paper/data/section4/ARPC_VCAP_up_hallmark.txt")


arpc.vcap.up <- arpc.vcap.up %>% 
  separate(col = 'Overlap',into = c('count','gene_set_count'),sep = '/',remove = F)

arpc.vcap.up$count <- as.numeric(arpc.vcap.up$count)
arpc.vcap.up <- arpc.vcap.up %>% mutate(gene_ratio=arpc.vcap.up$count/3855)


#arpc.vcap.up <- arpc.vcap.up[order(-arpc.vcap.up$count), ]
arpc.vcap.up <- arpc.vcap.up %>% filter(Adjusted.P.value < 0.05)
arpc.vcap.up <- arpc.vcap.up[order(arpc.vcap.up$P.value), ] %>% head(n=10)

ggplot(arpc.vcap.up, aes(x = gene_ratio, y = reorder(Term, count,decreasing = F))) +
  geom_point(aes(size = count, fill = P.value), shape = 21) +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(
    
    size = "Count",
    x = "Gene Ratio",
    y = "",
    title = "VCaP-vs-ARPC down-regulated"
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 38))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color="black", linetype="solid", size = 1),
        axis.text.x = element_text( size = 20,color = 'black'),
        axis.text.y = element_text(size = 20,color = 'black'),
        axis.title.x = element_text(size = 20,color = 'black'),
        legend.text = element_text(color = "black",size = 15),
        legend.title = element_text(color = "black",size = 15),
        axis.ticks.x = element_line(size=1),
        axis.ticks.y = element_line(size=1),
        plot.title = element_text(color = "black",face = "bold", size = 17,hjust = 0.5))
ggsave('./pcdata/paper/plot/section4/VCaP-vs-ARPC down-regulated.pdf', width = 10, height = 8)

#arpc.vcap.dn 2210 富集分析 --------------------------------------------------
arpc.vcap.dn <- read.delim("./pcdata/paper/data/section4/ARPC_VCAP_dn_hallmark.txt")


arpc.vcap.dn <- arpc.vcap.dn %>% 
  separate(col = 'Overlap',into = c('count','gene_set_count'),sep = '/',remove = F)

arpc.vcap.dn$count <- as.numeric(arpc.vcap.dn$count)
arpc.vcap.dn <- arpc.vcap.dn %>% mutate(gene_ratio=arpc.vcap.dn$count/2210)


#arpc.vcap.dn <- arpc.vcap.dn[order(-arpc.vcap.dn$count), ]
arpc.vcap.dn <- arpc.vcap.dn %>% filter(Adjusted.P.value < 0.05)
arpc.vcap.dn <- arpc.vcap.dn[order(arpc.vcap.dn$P.value), ] %>% head(n=10)

ggplot(arpc.vcap.dn, aes(x = gene_ratio, y = reorder(Term, count,decreasing = F))) +
  geom_point(aes(size = count, fill = P.value), shape = 21) +
  scale_fill_gradient(low = "red", high = "blue") +
  labs(
    
    size = "Count",
    x = "Gene Ratio",
    y = "",
    title = "VCaP-vs-ARPC up-regulated"
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 38))+
  theme_bw()+
  theme(panel.border = element_rect(fill = NA, color="black", linetype="solid", size = 1),
        axis.text.x = element_text( size = 20,color = 'black'),
        axis.text.y = element_text(size = 20,color = 'black'),
        axis.title.x = element_text(size = 20,color = 'black'),
        legend.text = element_text(color = "black",size = 15),
        legend.title = element_text(color = "black",size = 15),
        axis.ticks.x = element_line(size=1),
        axis.ticks.y = element_line(size=1),
        plot.title = element_text(color = "black",face = "bold", size = 17,hjust = 0.5))
ggsave('./pcdata/paper/plot/section4/VCaP-vs-ARPC up-regulated.pdf', width = 9, height = 8)
