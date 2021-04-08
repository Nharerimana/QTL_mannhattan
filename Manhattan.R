library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(grid)
library(gridExtra)
library(gtable)
library(readxl)
library(stringi)
library(openxlsx)
library(tibble)

# Make the Manhattan plot with ggplot
df<- read.csv("df.csv")
# Prepare the dataset
don <- df %>% 
    # Compute chromosome size
  group_by(SNP_chr) %>% 
  summarise(chr_len=max(SNP_pos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(sigeQTLs, ., by=c("SNP_chr"="SNP_chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(SNP_chr, SNP_pos) %>%
  mutate( BPcum= SNP_pos+tot) %>%

  # Add highlight and annotation information
  mutate( is_annotate=ifelse(-log10(pvalue)>8, "yes", "no")) 

# Prepare X axis
axisdf <- don %>% group_by(SNP_chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
# Make the plot
PanelA <- ggplot(don, aes(x=BPcum, y=-log10(pvalue))) +
  labs(x = "Chromosome", y = expression(paste('−log10(P'[eQTL],')')))+
  # Show all points
  geom_point( aes(color=as.factor(SNP_chr)), alpha=0.8, size=1.7) +
  scale_color_manual(values = rep(c("green4"), 22 )) +
  # custom X axis:
  scale_x_continuous( label = axisdf$SNP_chr, breaks= axisdf$center ) +
  # Add label using ggrepel to avoid overlapping
  geom_text_repel(data=subset(don, is_annotate=="yes"), aes(label=hgnc_symbol), size= 4,angle = 0, 
                    min.segment.length = 0 ,arrow = arrow(length = unit(0.015, "npc"))) +
  # Custom the theme:
  theme_bw() +theme(legend.position="none",
    axis.line.x=element_line(),
    axis.line.y=element_line(),
    axis.text.y=element_text(size=12),
    axis.text.x=element_text(size=12),
    axis.title=element_text(size=12),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank())