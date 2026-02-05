setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.5/")
library(vegan)
gene<-read.csv("Nitrogen gene abundance.csv",
               header = T,row.names = 1)
group<-read.csv("group.csv",header = T,row.names = 1)

gene2 <- merge(gene,group,by = "row.names")


gene3<-gene2[c(1,4:18)]
gene3$Anammox_Bacteria <- NULL
row.names(gene3) <- gene3$Row.names
gene3$Row.names <- NULL
#
variables <- gene3[1:12]
group_factors <- gene3[13:14]
set.seed(123) 

df <- gene3
permanova_nitrogen <- adonis2(variables ~ treatment, data = df, permutations = 9999)
print(permanova_nitrogen)

permanova_niche <- adonis2(variables ~ Niche, data = df, permutations = 9999)
print(permanova_niche)

r_value <- c(permanova_nitrogen$R2[1], permanova_niche$R2[1])
p_value <- c(permanova_nitrogen$`Pr(>F)`[1], permanova_niche$`Pr(>F)`[1])
sample <- c("Nitrogen", "Niche")

data <- data.frame(sample, r_value, p_value)

print(data)
head(data)

data$signif <- ifelse(data$p_value < 0.001, "***",
                      ifelse(data$p_value < 0.01, "**",
                             ifelse(data$p_value < 0.05, "*", "ns")))
print(data)

library(ggplot2)
plot <- ggplot(data, aes(x = sample, y = r_value, fill = sample)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  labs(x = "", y = "PERMANOVA R2", title = "Bar Plot of R Values by Sample") +
  geom_text(data = data, aes(x = sample, label = signif, y = max(r_value)*1.05), size = 4) +
  scale_fill_manual(values = c("Nitrogen" = "#093793","Niche" = "#DB2F2F")) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title = element_text(size = 14, hjust = 0.5, color = "#333333", face = "bold", margin = margin(b = 10)),
    axis.text.x = element_text(angle = 0, vjust = 0.5, face = "bold", color = "#333333"),
    plot.margin = margin(10, 10, 10, 10, unit = "pt"),
    legend.position = "none",
    text = element_text(family = "serif"),
    panel.background = element_blank(),
    plot.background = element_blank()
  )

plot
ggsave(plot = plot,filename = "Niche and Nitrogen variations in gene abundance.pdf",device = "pdf",dpi = 300,width = 3,height = 6)
