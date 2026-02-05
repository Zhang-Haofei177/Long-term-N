#######################Niche#######################################
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(reshape2)
library(psych)
setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.5/")
gene<-read.csv("Nitrogen gene abundance.csv",
               header = T,row.names = 1)
env<-read.csv("env.csv",
              header = T,row.names = 1)
env$Niches[env$Niches == "Root zone"] <- "1"
env$Niches[env$Niches == "Rhizosphere"] <- "2"
env$Niches[env$Niches == "Root"] <- "3"
env$Niches[env$Niches == "Nodules"] <- "4"
group<-read.csv("group.csv",header = T,row.names = 1)
gene2<-gene[3:ncol(gene)]
env2<-env[2:ncol(env)]
env2$Niches <- as.numeric(env2$Niches)

e_rz<-env2[env2$Nitrogen==0,]
g_rz<-gene2[rownames(e_rz),]

r_result <- corr.test(g_rz,e_rz, method = 'spearman')
r <- 
  r_result$r %>% 
  melt() %>% 
  set_names(c('Genes', 'Index', 'r'))

p <- 
  r_result$p %>% 
  melt() %>% 
  set_names(c('Genes', 'Index', 'P_value')) %>% 
  mutate(P_value_sig = case_when(P_value > 0.05 ~ " ",
                                 P_value <= 0.05 & P_value > 0.01 ~ "*",
                                 P_value <= 0.01 & P_value > 0.001 ~ "**",
                                 P_value <= 0.001 ~ "***",
                                 TRUE ~ NA_character_))

data <- cbind(r,p) %>% select(-(4:5))
head(data)

data_r<-data[1:3]
data_p<-data[c(1,2,5)]
N0_wide_r <- spread(data_r, Index, r)
N0_wide_r$Nitrogen <- "N0"
N0_wide_p <- spread(data_p, Index, P_value_sig)
N0_wide_p$Nitrogen <- "N0"

e_rh<-env2[env2$Nitrogen==200,]
g_rh<-gene2[rownames(e_rh),]

r_result <- corr.test(g_rh,e_rh, method = 'spearman')
r <- 
  r_result$r %>% 
  melt() %>% 
  set_names(c('Genes', 'Index', 'r'))

p <- 
  r_result$p %>% 
  melt() %>% 
  set_names(c('Genes', 'Index', 'P_value')) %>% 
  mutate(P_value_sig = case_when(P_value > 0.05 ~ " ",
                                 P_value <= 0.05 & P_value > 0.01 ~ "*",
                                 P_value <= 0.01 & P_value > 0.001 ~ "**",
                                 P_value <= 0.001 ~ "***",
                                 TRUE ~ NA_character_))

data <- cbind(r,p) %>% select(-(4:5))
head(data)

data_r<-data[1:3]
data_p<-data[c(1,2,5)]
N2_wide_r <- spread(data_r, Index, r)
N2_wide_r$Nitrogen <- "N2"
N2_wide_p <- spread(data_p, Index, P_value_sig)
N2_wide_p$Nitrogen <- "N2"

e_rt<-env2[env2$Nitrogen==400,]
g_rt<-gene2[rownames(e_rt),]

r_result <- corr.test(g_rt,e_rt, method = 'spearman')
r <- 
  r_result$r %>% 
  melt() %>% 
  set_names(c('Genes', 'Index', 'r'))

p <- 
  r_result$p %>% 
  melt() %>% 
  set_names(c('Genes', 'Index', 'P_value')) %>% 
  mutate(P_value_sig = case_when(P_value > 0.05 ~ " ",
                                 P_value <= 0.05 & P_value > 0.01 ~ "*",
                                 P_value <= 0.01 & P_value > 0.001 ~ "**",
                                 P_value <= 0.001 ~ "***",
                                 TRUE ~ NA_character_))

data <- cbind(r,p) %>% select(-(4:5))
head(data)

data_r<-data[1:3]
data_p<-data[c(1,2,5)]
N4_wide_r <- spread(data_r, Index, r)
N4_wide_r$Nitrogen <- "N4"
N4_wide_p <- spread(data_p, Index, P_value_sig)
N4_wide_p$Nitrogen <- "N4"

e_nod<-env2[env2$Nitrogen==600,]
g_nod<-gene2[rownames(e_nod),]

r_result <- corr.test(g_nod,e_nod, method = 'spearman')
r <- 
  r_result$r %>% 
  melt() %>% 
  set_names(c('Genes', 'Index', 'r'))

p <- 
  r_result$p %>% 
  melt() %>% 
  set_names(c('Genes', 'Index', 'P_value')) %>% 
  mutate(P_value_sig = case_when(P_value > 0.05 ~ " ",
                                 P_value <= 0.05 & P_value > 0.01 ~ "*",
                                 P_value <= 0.01 & P_value > 0.001 ~ "**",
                                 P_value <= 0.001 ~ "***",
                                 TRUE ~ NA_character_))

data <- cbind(r,p) %>% select(-(4:5))
head(data)

data_r<-data[1:3]
data_p<-data[c(1,2,5)]
N6_wide_r <- spread(data_r, Index, r)
N6_wide_r$Nitrogen <- "N6"
N6_wide_p <- spread(data_p, Index, P_value_sig)
N6_wide_p$Nitrogen <- "N6"

total <- rbind(N0_wide_r,N2_wide_r,N4_wide_r,N6_wide_r)
total <- total[1:3]
total_p <- rbind(N0_wide_p,N2_wide_p,N4_wide_p,N6_wide_p)
total_p <- total_p[1:3]
write.csv(total,"Gene Abundance_Niche_Spearman r.csv")
write.csv(total_p,"Gene Abundance_Niche_Spearman p.csv")
#####################Nitrogen##########################################
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(reshape2)
library(psych)
setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.5/")
gene<-read.csv("Nitrogen gene abundance.csv",
               header = T,row.names = 1)
env<-read.csv("env.csv",
              header = T,row.names = 1)

group<-read.csv("group.csv",header = T,row.names = 1)
gene2<-gene[3:ncol(gene)]
env2<-env[2:ncol(env)]
env2$Nitrogen <- as.numeric(env2$Nitrogen)

e_rz<-env2[env2$Niches=="Root zone",]
g_rz<-gene2[rownames(e_rz),]

r_result <- corr.test(g_rz,e_rz[2:23], method = 'spearman')
r <- 
  r_result$r %>% 
  melt() %>% 
  set_names(c('Genes', 'Index', 'r'))

p <- 
  r_result$p %>% 
  melt() %>% 
  set_names(c('Genes', 'Index', 'P_value')) %>% 
  mutate(P_value_sig = case_when(P_value > 0.05 ~ " ",
                                 P_value <= 0.05 & P_value > 0.01 ~ "*",
                                 P_value <= 0.01 & P_value > 0.001 ~ "**",
                                 P_value <= 0.001 ~ "***",
                                 TRUE ~ NA_character_))

data <- cbind(r,p) %>% select(-(4:5))
head(data)

data_r<-data[1:3]
data_p<-data[c(1,2,5)]
RZS_wide_r <- spread(data_r, Index, r)
RZS_wide_r$Niche <- "RZS"
RZS_wide_p <- spread(data_p, Index, P_value_sig)
RZS_wide_p$Niche <- "RZS"

e_rh<-env2[env2$Niches=="Rhizosphere",]
g_rh<-gene2[rownames(e_rh),]

r_result <- corr.test(g_rh,e_rh[2:23], method = 'spearman')
r <- 
  r_result$r %>% 
  melt() %>% 
  set_names(c('Genes', 'Index', 'r'))

p <- 
  r_result$p %>% 
  melt() %>% 
  set_names(c('Genes', 'Index', 'P_value')) %>% 
  mutate(P_value_sig = case_when(P_value > 0.05 ~ " ",
                                 P_value <= 0.05 & P_value > 0.01 ~ "*",
                                 P_value <= 0.01 & P_value > 0.001 ~ "**",
                                 P_value <= 0.001 ~ "***",
                                 TRUE ~ NA_character_))

data <- cbind(r,p) %>% select(-(4:5))
head(data)

data_r<-data[1:3]
data_p<-data[c(1,2,5)]
RS_wide_r <- spread(data_r, Index, r)
RS_wide_r$Niche <- "RS"
RS_wide_p <- spread(data_p, Index, P_value_sig)
RS_wide_p$Niche <- "RS"

e_rt<-env2[env2$Niches=="Root",]
g_rt<-gene2[rownames(e_rt),]

r_result <- corr.test(g_rt,e_rt[2:23], method = 'spearman')
r <- 
  r_result$r %>% 
  melt() %>% 
  set_names(c('Genes', 'Index', 'r'))

p <- 
  r_result$p %>% 
  melt() %>% 
  set_names(c('Genes', 'Index', 'P_value')) %>% 
  mutate(P_value_sig = case_when(P_value > 0.05 ~ " ",
                                 P_value <= 0.05 & P_value > 0.01 ~ "*",
                                 P_value <= 0.01 & P_value > 0.001 ~ "**",
                                 P_value <= 0.001 ~ "***",
                                 TRUE ~ NA_character_))

data <- cbind(r,p) %>% select(-(4:5))
head(data)

data_r<-data[1:3]
data_p<-data[c(1,2,5)]
RE_wide_r <- spread(data_r, Index, r)
RE_wide_r$Niche <- "RE"
RE_wide_p <- spread(data_p, Index, P_value_sig)
RE_wide_p$Niche <- "RE"

e_nod<-env2[env2$Niches=="Nodules",]
g_nod<-gene2[rownames(e_nod),]

r_result <- corr.test(g_nod,e_nod[2:23], method = 'spearman')
r <- 
  r_result$r %>% 
  melt() %>% 
  set_names(c('Genes', 'Index', 'r'))

p <- 
  r_result$p %>% 
  melt() %>% 
  set_names(c('Genes', 'Index', 'P_value')) %>% 
  mutate(P_value_sig = case_when(P_value > 0.05 ~ " ",
                                 P_value <= 0.05 & P_value > 0.01 ~ "*",
                                 P_value <= 0.01 & P_value > 0.001 ~ "**",
                                 P_value <= 0.001 ~ "***",
                                 TRUE ~ NA_character_))

data <- cbind(r,p) %>% select(-(4:5))
head(data)

data_r<-data[1:3]
data_p<-data[c(1,2,5)]
NE_wide_r <- spread(data_r, Index, r)
NE_wide_r$Niche <- "NE"
NE_wide_p <- spread(data_p, Index, P_value_sig)
NE_wide_p$Niche <- "NE"

total <- rbind(RZS_wide_r,RS_wide_r,RE_wide_r,NE_wide_r)
total <- total[c("Genes","Nitrogen","Niche")]
total_p <- rbind(RZS_wide_p,RS_wide_p,RE_wide_p,NE_wide_p)
total_p <- total_p[c("Genes","Nitrogen","Niche")]
write.csv(total,"Gene Abundance_Nitrogen_Spearman r.csv")
write.csv(total_p,"Gene Abundance_Nitrogen_Spearman p.csv")
