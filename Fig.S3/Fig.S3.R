setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.S3")
#### Load data, import two essential packages
library(vegan)
library(ggplot2)
############## Rhizobiales #########################################################################################################
group <- read.csv("group.csv")
colnames(group)[1] <- "names"
group <- group[grepl("N0|N2|N4|N6", group$treatment), ]
a1 <- read.csv("Rhizobia_ASV.csv", row.names = 1, header = T)
a1 <- a1[, grepl("N0|N2|N4|N6", colnames(a1))]
a1 <- t(a1)

ly.rare <- function(a, n) {
  b <- rowSums(a)
  
  i = 1
  {
    n.i <- as.vector(c(seq(0, b[i], n), b[i]))
    rare.i <- as.integer(rarefy(a[i, ], n.i))
    rare <- cbind(rep(names(b)[i], length(n.i)), n.i, rare.i)
  }
  
  for (i in 2:length(b)) {
    n.i <- as.vector(c(seq(0, b[i], n), b[i]))
    rare.i <- as.integer(rarefy(a[i, ], n.i))
    rare.ni <- cbind(rep(names(b)[i], length(n.i)), n.i, rare.i)
    rare <- rbind(rare, rare.ni)
  }
  
  colnames(rare)[1:3] <- c("names", "reads", "value")
  rare <- as.data.frame(rare)
  rare$reads <- as.numeric(as.character(rare$reads))
  rare$value <- as.numeric(as.character(rare$value))
  rare
}
# Usage: ly.rare(a, n), where 'a' is the OTU abundance table and 'n' is the step size for read sampling.
# Run the above function, input data (using non-normalized data as an example) and set the sampling step size. Here, set every 2000 reads as a sampling interval.
rare <- ly.rare(a1, 2000)

## Add grouping information to the new data frame
rare.new <- merge(rare, group, by = "names")
## Generate plot using ggplot, a basic ggplot plotting function. Fine-tune as needed.
rare.new$treatment <- factor(rare.new$treatment, levels = c("N0", "N2", "N4", "N6"))
rare.new$Niche <- factor(rare.new$Niche, levels = c("Original", "Rootzone", "Rhizosphere", "Root", "Nodule"),
                         labels = c("MFS", "RZS", "RS", "RE", "NE"))
p <- ggplot(rare.new, aes(reads, value, group = names, colour = Niche)) +
  geom_line(lwd = 0.3) +
  labs(x = "Number of sequences", y = "Number of OTUs") +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black')) +
  guides(colour = guide_legend(title = NULL)) +
  theme(legend.key = element_rect(fill = "white")) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.text = element_text(size = 8))
p
ggsave(plot = p, filename = "Rhizobia_Rarefaction_curves(Niche).pdf", device = "pdf", dpi = 300, width = 8, height = 4)
p <- ggplot(rare.new, aes(reads, value, group = names, colour = treatment)) +
  geom_line(lwd = 0.3) +
  labs(x = "Number of sequences", y = "Number of OTUs") +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black')) +
  guides(colour = guide_legend(title = NULL)) +
  theme(legend.key = element_rect(fill = "white")) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.text = element_text(size = 8))
p
ggsave(plot = p, filename = "Rhizobia_Rarefaction_curves(Treatment).pdf", device = "pdf", dpi = 300, width = 8, height = 4)
########### Bacteria ############################################################################################################
group <- read.csv("group.csv")
colnames(group)[1] <- "names"
group <- group[grepl("N0|N2|N4|N6", group$treatment), ]
a1 <- read.csv("Bacteria_ASV.csv", row.names = 1, header = T)
a1 <- a1[, grepl("N0|N2|N4|N6", colnames(a1))]
a1 <- t(a1)

ly.rare <- function(a, n) {
  b <- rowSums(a)
  
  i = 1
  {
    n.i <- as.vector(c(seq(0, b[i], n), b[i]))
    rare.i <- as.integer(rarefy(a[i, ], n.i))
    rare <- cbind(rep(names(b)[i], length(n.i)), n.i, rare.i)
  }
  
  for (i in 2:length(b)) {
    n.i <- as.vector(c(seq(0, b[i], n), b[i]))
    rare.i <- as.integer(rarefy(a[i, ], n.i))
    rare.ni <- cbind(rep(names(b)[i], length(n.i)), n.i, rare.i)
    rare <- rbind(rare, rare.ni)
  }
  
  colnames(rare)[1:3] <- c("names", "reads", "value")
  rare <- as.data.frame(rare)
  rare$reads <- as.numeric(as.character(rare$reads))
  rare$value <- as.numeric(as.character(rare$value))
  rare
}

rare <- ly.rare(a1, 2000)

rare.new <- merge(rare, group, by = "names")

rare.new$treatment <- factor(rare.new$treatment, levels = c("N0", "N2", "N4", "N6"))
rare.new$Niche <- factor(rare.new$Niche, levels = c("Original", "Rootzone", "Rhizosphere", "Root", "Nodule"),
                         labels = c("MFS", "RZS", "RS", "RE", "NE"))
p <- ggplot(rare.new, aes(reads, value, group = names, colour = Niche)) +
  geom_line(lwd = 0.3) +
  labs(x = "Number of sequences", y = "Number of OTUs") +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black')) +
  guides(colour = guide_legend(title = NULL)) +
  theme(legend.key = element_rect(fill = "white")) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.text = element_text(size = 8))
p
ggsave(plot = p, filename = "Bacteria_Rarefaction_curves(Niche).pdf", device = "pdf", dpi = 300, width = 8, height = 4)
p <- ggplot(rare.new, aes(reads, value, group = names, colour = treatment)) +
  geom_line(lwd = 0.3) +
  labs(x = "Number of sequences", y = "Number of OTUs") +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black')) +
  guides(colour = guide_legend(title = NULL)) +
  theme(legend.key = element_rect(fill = "white")) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.text = element_text(size = 8))
p
ggsave(plot = p, filename = "Bacteria_Rarefaction_curves(Treatment).pdf", device = "pdf", dpi = 300, width = 8, height = 4)
############## AMF #########################################################################################################
group <- read.csv("group.csv")
colnames(group)[1] <- "names"
group <- group[grepl("N0|N2|N4|N6", group$treatment), ]
a1 <- read.csv("AMF_ASV.csv", row.names = 1, header = T)
a1 <- a1[, grepl("N0|N2|N4|N6", colnames(a1))]
a1 <- t(a1)


ly.rare <- function(a, n) {
  b <- rowSums(a)
  
  i = 1
  {
    n.i <- as.vector(c(seq(0, b[i], n), b[i]))
    rare.i <- as.integer(rarefy(a[i, ], n.i))
    rare <- cbind(rep(names(b)[i], length(n.i)), n.i, rare.i)
  }
  
  for (i in 2:length(b)) {
    n.i <- as.vector(c(seq(0, b[i], n), b[i]))
    rare.i <- as.integer(rarefy(a[i, ], n.i))
    rare.ni <- cbind(rep(names(b)[i], length(n.i)), n.i, rare.i)
    rare <- rbind(rare, rare.ni)
  }
  
  colnames(rare)[1:3] <- c("names", "reads", "value")
  rare <- as.data.frame(rare)
  rare$reads <- as.numeric(as.character(rare$reads))
  rare$value <- as.numeric(as.character(rare$value))
  rare
}
# Usage: ly.rare(a, n), where 'a' is the OTU abundance table and 'n' is the step size for read sampling.
# Run the above function, input data (using non-normalized data as an example) and set the sampling step size. Here, set every 2000 reads as a sampling interval.
rare <- ly.rare(a1, 2000)

## Add grouping information to the new data frame
rare.new <- merge(rare, group, by = "names")
## Generate plot using ggplot, a basic ggplot plotting function. Fine-tune as needed.
rare.new$treatment <- factor(rare.new$treatment, levels = c("N0", "N2", "N4", "N6"))
rare.new$Niche <- factor(rare.new$Niche, levels = c("Original", "Rootzone", "Rhizosphere", "Root", "Nodule"),
                         labels = c("MFS", "RZS", "RS", "RE", "NE"))
p <- ggplot(rare.new, aes(reads, value, group = names, colour = Niche)) +
  geom_line(lwd = 0.3) +
  labs(x = "Number of sequences", y = "Number of OTUs") +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black')) +
  guides(colour = guide_legend(title = NULL)) +
  theme(legend.key = element_rect(fill = "white")) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.text = element_text(size = 8))
p
ggsave(plot = p, filename = "AMF_Rarefaction_curves(Niche).pdf", device = "pdf", dpi = 300, width = 8, height = 4)
p <- ggplot(rare.new, aes(reads, value, group = names, colour = treatment)) +
  geom_line(lwd = 0.3) +
  labs(x = "Number of sequences", y = "Number of OTUs") +
  theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black')) +
  guides(colour = guide_legend(title = NULL)) +
  theme(legend.key = element_rect(fill = "white")) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(size = 10)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 12)) +
  theme(legend.text = element_text(size = 8))
p
ggsave(plot = p, filename = "AMF_Rarefaction_curves(Treatment).pdf", device = "pdf", dpi = 300, width = 8, height = 4)