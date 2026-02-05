#######################Niche LMM#######################################
rm(list=ls())
setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.5/")
library(tidyverse)
library(magrittr)
library(ape)
library(base)
library(vegan)
library(phytools)
library(scales)
library(lme4)
library(Matrix)
library(dplyr)
setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.5/")
gene<-read.csv("Nitrogen gene abundance.csv",
               header = T,row.names = 1)
env<-read.csv("env.csv",
              header = T,row.names = 1)
env <- env %>%
  mutate(Niche2 = case_when(
    Niches == "Root zone" ~ 1,
    Niches == "Rhizosphere" ~ 2,
    Niches == "Root" ~ 3,
    Niches == "Nodules" ~ 4,
    TRUE ~ NA_real_  
  ))
group<-read.csv("group.csv",header = T,row.names = 1)
gene2<-gene[3:ncol(gene)]
env2<-env[3:ncol(env)]

treatused <- env2
divindex <- gene2
divindex <- scale(divindex)

divs2 <- do.call(cbind, lapply(1:ncol(divindex), function(j){
  message("Now j=", j, " in ", ncol(divindex), ". ", date())
  
  if (length(unique(divindex[, j])) < 3){
    result <- rep(NA, 38)
  } else {
    div <- data.frame(divtest = divindex[, j], treatused)
    
    fm1 <- lmer(divtest ~ Niche2 + (1 | Nitrogen), data = div)
    
    presult <- car::Anova(fm1, type=2)
    coefs <- coef(summary(fm1))[ , "Estimate"]
    names(coefs) <- paste0(names(coefs), ".mean")
    
    SEvalues <- coef(summary(fm1))[ , "Std. Error"]
    names(SEvalues) <- paste0(names(SEvalues), ".se")
    
    tvalues <- coef(summary(fm1))[ , "t value"]
    names(tvalues) <- paste0(names(tvalues), ".t")
    
    chisqP <- c(presult[,1], presult[,3])
    names(chisqP) <- c(paste0(row.names(presult), ".chisq"), paste0(row.names(presult), ".P"))
    
    result <- c(coefs, tvalues, SEvalues, chisqP)
    
    print(length(result))
    
    if(length(result) < 8) {
      result <- c(result, rep(NA, 8 - length(result)))
    }
  }
  
  result
}))


colnames(divs2) <- colnames(divindex)
divs2<-divs2[1:8,]
write.csv(divs2,file = "Niche_LMM result.csv")
#################Nitrogen LMM##############################
setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.5/")
gene<-read.csv("Nitrogen gene abundance.csv",
               header = T,row.names = 1)
env<-read.csv("env.csv",
              header = T,row.names = 1)
env <- env %>%
  mutate(Niche2 = case_when(
    Niches == "Root zone" ~ 1,
    Niches == "Rhizosphere" ~ 2,
    Niches == "Root" ~ 3,
    Niches == "Nodules" ~ 4,
    TRUE ~ NA_real_  
  ))
group<-read.csv("group.csv",header = T,row.names = 1)
gene2<-gene[3:ncol(gene)]
env2<-env[3:ncol(env)]

treatused <- env2
divindex <- gene2
divindex <- scale(divindex)

divs2 <- do.call(cbind, lapply(1:ncol(divindex), function(j){
  message("Now j=", j, " in ", ncol(divindex), ". ", date())
  
  if (length(unique(divindex[, j])) < 3){
    result <- rep(NA, 38)
  } else {
    div <- data.frame(divtest = divindex[, j], treatused)
    
    fm1 <- lmer(divtest ~ Nitrogen + (1 | Niche2), data = div)
    
    presult <- car::Anova(fm1, type=2)
    coefs <- coef(summary(fm1))[ , "Estimate"]
    names(coefs) <- paste0(names(coefs), ".mean")
    
    SEvalues <- coef(summary(fm1))[ , "Std. Error"]
    names(SEvalues) <- paste0(names(SEvalues), ".se")
    
    tvalues <- coef(summary(fm1))[ , "t value"]
    names(tvalues) <- paste0(names(tvalues), ".t")
    
    chisqP <- c(presult[,1], presult[,3])
    names(chisqP) <- c(paste0(row.names(presult), ".chisq"), paste0(row.names(presult), ".P"))
    
    result <- c(coefs, tvalues, SEvalues, chisqP)
    

    print(length(result))
    

    if(length(result) < 8) {
      result <- c(result, rep(NA, 8 - length(result)))
    }
  }
  
  result
}))

colnames(divs2) <- colnames(divindex)
divs2<-divs2[1:8,]
write.csv(divs2,file = "Nitrogen_LMM result.csv")
