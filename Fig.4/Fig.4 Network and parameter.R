path = "C:/Users/ADMIN/Desktop/Longterm N code/Fig.4/"
setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.4/")
library(ggplot2)
library(ggpubr)
library(ggalt)
library(phyloseq)
library(igraph)
library(network)
library(sna)
library(tidyverse)
library(ggClusterNet)
metadata<-read.csv("metadata.csv",header = T,row.names = 1)
metadata <- metadata[!grepl("Original",metadata$Location),]

AMF<-read.csv("AMF_ASV_rarefaction.csv",header = T,row.names = 1)
AMF<-AMF[grepl("N0|N2|N4|N6",colnames(AMF)),]
AMF<-AMF[rowSums(AMF)!=0,]
RHI<-read.csv("Rhizobia_ASV_rarefaction.csv",header = T,row.names = 1)
RHI<-RHI[grepl("N0|N2|N4|N6",colnames(RHI)),]
RHI<-RHI[rowSums(RHI)!=0,]
RHI<-RHI[,rownames(metadata)]#调整列顺序一致
AMF<-AMF[,rownames(metadata)]#调整列顺序一致
AMF_TAX<-read.csv("AMF_tax.csv",header = T,row.names = 1)
RHI_TAX<-read.csv("Rhizobia_tax.csv",header = T,row.names = 1)
AMF_TAX<-AMF_TAX[rownames(AMF),]
RHI_TAX<-RHI_TAX[rownames(RHI),]
AMF_TAX<-AMF_TAX[,colnames(RHI_TAX)]


ps16 = phyloseq(sample_data(metadata),
                otu_table(as.matrix(RHI), taxa_are_rows=TRUE),
                tax_table(as.matrix(RHI_TAX)))
psIT=phyloseq(sample_data(metadata),
              otu_table(as.matrix(AMF), taxa_are_rows=TRUE),
              tax_table(as.matrix(AMF_TAX)))
str(ps16)

ps.merge <- ggClusterNet::merge16S_ITS(ps16s = ps16,
                                       psITS= psIT,
                                       NITS = 500,
                                       N16s = 500)
map =  phyloseq::sample_data(ps.merge)
head(map)
phyloseq::sample_data(ps.merge) <- map
#################Niche Networks#########################
library(sna)
library(ggClusterNet)
library(igraph)

result <- ggClusterNet::corBionetwork(ps = ps.merge,
                                      N = 0,
                                      r.threshold = 0.8, 
                                      p.threshold = 0.05,
                                      big = T,
                                      group = "Location",
                                      path = paste0(path,"output_Niche"),
                                      fill = "Phylum", 
                                      size = "igraph.degree", 
                                      scale = TRUE, 
                                      bio = TRUE, 
                                      zipi = T, 
                                      step = 100, 
                                      width = 18,
                                      label = TRUE,
                                      height = 10,
                                      layout_net = "model_maptree2"
)
p = result[[1]]
p
data = result[[2]]
write.csv(data,paste0(path,"output_Niche","/co-occurrence_Grobel_net.csv"))
#################Nitrogen Networks#########################
library(sna)
library(ggClusterNet)
library(igraph)
result <- ggClusterNet::corBionetwork(ps = ps.merge,
                                      N = 0,
                                      r.threshold = 0.8, 
                                      p.threshold = 0.05,
                                      big = T,
                                      group = "Nitrogen",
                                      path = paste0(path,"output_Nitrogen"),
                                      fill = "Phylum", 
                                      size = "igraph.degree", 
                                      scale = TRUE, 
                                      bio = TRUE, 
                                      zipi = T, 
                                      step = 100, 
                                      width = 18,
                                      label = TRUE,
                                      height = 10,
                                      layout_net = "model_maptree2"
)
p = result[[1]]
p
data = result[[2]]
write.csv(data,paste0(path,"output_Nitrogen","/co-occurrence_Grobel_net.csv"))
#############################Network_parameter_plot###################################################
#######Niche############
setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.4/output_Niche/")
a <- read.csv("co-occurrence_Grobel_net.csv",header = T,row.names = 1)
head(a)
parameters <- c("average.degree","average.path.length","connectance")
a2 <- as.data.frame(t(a[parameters,]))
rownames(a2) <- c("RZS","RS","RE","NE")
a2$network <- rownames(a2)
a2$group <- c("Niche")
a2$network <- factor(a2$network,levels = c("RZS","RS","RE","NE"))
phy.cols<-c("#85BA8F")
p_Nice_1<-ggplot(data=a2,aes(x=network,y=average.degree,
                       group=group,colour=group))+
  geom_point(size=3)+
  labs(x="",y="Degree")+
  geom_xspline(spline_shape=-0.5)+
  scale_x_discrete(limits=c("RZS","RS","RE","NE"))+
  scale_colour_manual(values=phy.cols)+
  theme_bw()+
  theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8))+
  theme(axis.title.y=element_text(size=12))+theme(axis.title.x=element_text(size=12))+
  theme(legend.title=element_text(size=5),legend.text=element_text(size=5))+theme(legend.position="bottom")
p_Nice_1

p_Nice_2<-ggplot(data=a2,aes(x=network,y=average.path.length,
                             group=group,colour=group))+
  geom_point(size=3)+
  labs(x="",y="Average path length")+
  geom_xspline(spline_shape=-0.5)+
  scale_x_discrete(limits=c("RZS","RS","RE","NE"))+
  scale_colour_manual(values=phy.cols)+
  theme_bw()+
  theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8))+
  theme(axis.title.y=element_text(size=12))+theme(axis.title.x=element_text(size=12))+
  theme(legend.title=element_text(size=5),legend.text=element_text(size=5))+theme(legend.position="bottom")
p_Nice_2

p_Nice_3<-ggplot(data=a2,aes(x=network,y=connectance,
                             group=group,colour=group))+
  geom_point(size=3)+
  labs(x="",y="Connectance")+
  geom_xspline(spline_shape=-0.5)+
  scale_x_discrete(limits=c("RZS","RS","RE","NE"))+
  scale_colour_manual(values=phy.cols)+
  theme_bw()+
  theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8))+
  theme(axis.title.y=element_text(size=12))+theme(axis.title.x=element_text(size=12))+
  theme(legend.title=element_text(size=5),legend.text=element_text(size=5))+theme(legend.position="bottom")
p_Nice_3
########Nitrogen#########
setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.4/output_Nitrogen/")
a <- read.csv("co-occurrence_Grobel_net.csv",header = T,row.names = 1)
head(a)
parameters <- c("average.degree","average.path.length","connectance")
a2 <- as.data.frame(t(a[parameters,]))
a2$network <- rownames(a2)
a2$group <- c("Nitrogen")
a2$network <- factor(a2$network,levels = c("N0","N2","N4","N6"))
phy.cols<-c("#A3C8DC")
p_Nitrogen_1<-ggplot(data=a2,aes(x=network,y=average.degree,
                             group=group,colour=group))+
  geom_point(size=3)+
  labs(x="",y="Degree")+
  geom_xspline(spline_shape=-0.5)+
  scale_x_discrete(limits=c("N0","N2","N4","N6"))+
  scale_colour_manual(values=phy.cols)+
  theme_bw()+
  theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8))+
  theme(axis.title.y=element_text(size=12))+theme(axis.title.x=element_text(size=12))+
  theme(legend.title=element_text(size=5),legend.text=element_text(size=5))+theme(legend.position="bottom")
p_Nitrogen_1

p_Nitrogen_2<-ggplot(data=a2,aes(x=network,y=average.path.length,
                             group=group,colour=group))+
  geom_point(size=3)+
  labs(x="",y="Average path length")+
  geom_xspline(spline_shape=-0.5)+
  scale_x_discrete(limits=c("N0","N2","N4","N6"))+
  scale_colour_manual(values=phy.cols)+
  theme_bw()+
  theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8))+
  theme(axis.title.y=element_text(size=12))+theme(axis.title.x=element_text(size=12))+
  theme(legend.title=element_text(size=5),legend.text=element_text(size=5))+theme(legend.position="bottom")
p_Nitrogen_2

p_Nitrogen_3<-ggplot(data=a2,aes(x=network,y=connectance,
                             group=group,colour=group))+
  geom_point(size=3)+
  labs(x="",y="Connectance")+
  geom_xspline(spline_shape=-0.5)+
  scale_x_discrete(limits=c("N0","N2","N4","N6"))+
  scale_colour_manual(values=phy.cols)+
  theme_bw()+
  theme(axis.text.x=element_text(size=8),axis.text.y=element_text(size=8))+
  theme(axis.title.y=element_text(size=12))+theme(axis.title.x=element_text(size=12))+
  theme(legend.title=element_text(size=5),legend.text=element_text(size=5))+theme(legend.position="bottom")
p_Nitrogen_3

library(patchwork)
combined_plot <- (p_Nice_1 / p_Nice_2 / p_Nice_3) | 
  (p_Nitrogen_1 / p_Nitrogen_2 / p_Nitrogen_3)
print(combined_plot)
setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.4/")
ggsave(filename = "Network_parameter_plot.pdf", plot = combined_plot,device = "pdf", width = 5, height = 8, dpi = 300)
