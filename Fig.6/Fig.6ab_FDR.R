rm(list=ls())
setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.6/")
df <- read.csv("Relative changes in edaphic factors.csv")
df <- df[1:21, 1:2]
df$sample <- factor(df$sample, levels = c(df$sample))
bar_colors <- ifelse(df$EFC < 0, "#375DA2", "#EE7154")
p <- ggplot(data = df, aes(x = sample, y = EFC * 100)) +
  geom_col(fill = bar_colors, color = NA) +
  labs(x = '', y = "Environment factor changes (%)") +
  ggtitle("") +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    plot.title = element_text(size = 14, hjust = 0.5, color = "black", 
                              face = "bold", margin = margin(b = 10)),  
    axis.text.x = element_text(angle = 90, vjust = 0.5, face = "bold",
                               color = "black", size = 6),
    plot.margin = margin(10, 10, 10, 10, unit = "pt"),
    legend.position = "none",
    text = element_text(size = 6)
  )

p
ggsave("Fig_6a.pdf",plot = p,device = "pdf",dpi = 200,width =70 ,height = 50,units = "mm")
###################################################
library(tidyverse)
library(ggplot2)
library(reshape2)
library(psych)
df <- read.csv("The difference in soil factors between RZS and MFS.csv", header = T, row.names = 1)

mydata <- df[3:22]
mydata$pH <- NULL
mydata <- scale(mydata)
head(mydata)
Zcore <- (mydata - mean(mydata)) / sd(mydata)
Zcore <- as.data.frame(Zcore)
Zcore$deltaEMF <- rowSums(Zcore) / 19
Zcore <- Zcore[rownames(df), ]
df$deltaEMF <- Zcore$deltaEMF

df$EMF <- NULL

# ----- 计算 Pearson 相关 -----
r_result <- corr.test(df$Nitrogen, df[3:23], method = 'pearson')

# ----- 对 p 值进行 FDR 校正 -----
# r_result$p 是矩阵，先转成向量校正再变回矩阵
p_adj <- p.adjust(as.vector(r_result$p), method = "fdr")
dim(p_adj) <- dim(r_result$p)
dimnames(p_adj) <- dimnames(r_result$p)

r <- 
  r_result$r %>% 
  melt() %>% 
  set_names(c('Genes', 'Index', 'r'))

p <- 
  p_adj %>%   # 使用校正后的 p 值矩阵
  melt() %>% 
  set_names(c('Genes', 'Index', 'P_value')) %>% 
  mutate(P_value_sig = case_when(P_value > 0.05 ~ " ",
                                 P_value <= 0.05 & P_value > 0.01 ~ "*",
                                 P_value <= 0.01 & P_value > 0.001 ~ "**",
                                 P_value <= 0.001 ~ "***",
                                 TRUE ~ NA_character_))

data <- cbind(r, p) %>% select(-(4:5))
head(data)

df <- data[order(-data$r), ]

pos <- df[df$r > 0, ]
neg <- df[df$r < 0, ]
df$Index <- factor(df$Index, levels = df$Index)
df$max <- max(df$r)
df$min <- min(df$r)
pos <- df[df$r >= 0, ]
neg <- df[df$r < 0, ]
pos$label <- pos$r + pos$max * 0.2
neg$label <- neg$r + neg$min * 0.2
df <- rbind(pos, neg)
df$Index <- factor(df$Index, levels = rev(df$Index))

bar_colors <- ifelse(df$r < 0, "#375DA2", "#EE7154")
p <- ggplot(data = df, aes(x = Index, y = r)) +
  geom_col(fill = bar_colors, color = NA) +
  labs(x = '', y = "Pearson's r") +
  ggtitle("") +
  geom_text(aes(label = P_value_sig, y = label), size = 2) +
  theme_classic() +  
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),  
        plot.title = element_text(size = 14, hjust = 0.5, color = "black", 
                                  face = "bold", margin = margin(b = 10)),  
        axis.text.x = element_text(angle = 90, vjust = 0.5, face = "bold",
                                   color = "black", size = 6),
        plot.margin = margin(10, 10, 10, 10, unit = "pt"),
        legend.position = "none") +
  theme(text = element_text(size = 6))

p
# 输出文件名添加 _FDR 后缀，以便区分
ggsave("Fig_6b_FDR.pdf", plot = p, device = "pdf", dpi = 200, width = 70, height = 50, units = "mm")
