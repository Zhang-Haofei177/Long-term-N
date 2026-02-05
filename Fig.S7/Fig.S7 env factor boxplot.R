setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.S7/")
rm(list = ls())
a <- read.csv("env.csv", header = TRUE)
a$Treatment <- NULL
library(reshape2)
library(multcomp)
library(ggplot2)
library(tidyverse)
library(dplyr)

# Extract Nitrogen treatment samples
Nitrogen <- a[grepl("N0|N2|N4|N6", a$sample), ]

### Original samples
Nitrogen_O <- subset(Nitrogen, Nitrogen$Niche == "Original")
Nitrogen_z <- subset(Nitrogen, Nitrogen$Niche == "Root_zone")
Nitrogen2 <- Nitrogen_O
data <- Nitrogen2
values <- colnames(Nitrogen2)[4:24]
group <- c("group")

stat_anova <- NULL
data_list <- NULL

for (i in group) {
  for (j in values) {
    dat <- data[c(i, j)]
    names(dat) <- c('class', 'var')
    dat$class <- factor(dat$class)
    
    fit <- aov(var ~ class, dat)
    p_value <- summary(fit)[[1]][1, 5]
    
    # Determine significance level
    if (p_value < 0.001) {
      sig <- '***'
    } else if (p_value >= 0.001 & p_value < 0.01) {
      sig <- '**'
    } else if (p_value >= 0.01 & p_value < 0.05) {
      sig <- '*'
    } else {
      sig <- ''
    }
    
    stat_anova <- rbind(stat_anova, c(i, j, p_value, sig))
    
    # Tukey HSD test (multcomp package), multiple comparisons
    tuk <- cld(glht(fit, alternative = 'two.sided', linfct = mcp(class = 'Tukey')), sig = p, decreasing = TRUE)
    
    # Extract significance letters "abc" from cld() output
    sig <- data.frame(tuk$mcletters$Letters, stringsAsFactors = FALSE)
    names(sig) <- 'sig'
    sig$class <- rownames(sig)
    sig$var <- j
    
    dat_new <- dat %>% group_by(class) %>% mutate(mean = mean(var), sd = sd(var))
    dat_new2 <- unique(dat_new[, c(1, 3, 4)])
    dat_new2$class <- as.character(dat_new2$class)
    dat_new2$group <- i
    
    dat_new2 <- merge(dat_new2, sig, by = 'class')
    dat_new2 <- dat_new2[c(4, 6, 1, 2, 3, 5)]
    
    data_list <- rbind(data_list, dat_new2)
  }
}

stat_anova
data_list
ori <- data_list

library(rstatix)
Nitrogen2 <- Nitrogen_O
Nitrogen_l <- melt(Nitrogen2, id.vars = c("sample", "group", "Niche"), variable.name = "Variable", value.name = "Value")

library(dplyr)
test <- data_list
test <- test[, c(2:6)]
colnames(test)[1] <- c("Variable")
colnames(test)[2] <- c("group")

new_df <- Nitrogen2 %>% group_by(group) %>% summarize_all(max)
Nitrogen_m <- melt(new_df, id.vars = c("group", "sample", "Niche"), variable.name = "Variable", value.name = "max")
test$max <- Nitrogen_m$max[match(paste(test$Variable, test$group, sep = "_"), paste(Nitrogen_m$Variable, Nitrogen_m$group, sep = "_"))]

new_df2 <- Nitrogen2 %>% summarize_all(max)
Nitrogen_x <- melt(new_df2, id.vars = c("group", "sample", "Niche"), variable.name = "Variable", value.name = "max")
test$x <- Nitrogen_x$max[match(test$Variable, Nitrogen_x$Variable)]

colors <- c("#88CCEE", "#44AA99", "#117733", "#999933",
            "#BFDDEE", "#78C1D6", "#459F80", "#C5C297")

ori_l <- Nitrogen_l
ori_test <- test

######## RZS samples ####
Nitrogen2 <- Nitrogen_z
data <- Nitrogen2
values <- colnames(Nitrogen2)[4:24]
group <- c("group")

stat_anova <- NULL
data_list <- NULL

for (i in group) {
  for (j in values) {
    dat <- data[c(i, j)]
    names(dat) <- c('class', 'var')
    dat$class <- factor(dat$class)
    
    fit <- aov(var ~ class, dat)
    p_value <- summary(fit)[[1]][1, 5]
    
    # Determine significance level
    if (p_value < 0.001) {
      sig <- '***'
    } else if (p_value >= 0.001 & p_value < 0.01) {
      sig <- '**'
    } else if (p_value >= 0.01 & p_value < 0.05) {
      sig <- '*'
    } else {
      sig <- ''
    }
    
    stat_anova <- rbind(stat_anova, c(i, j, p_value, sig))
    
    # Tukey HSD test (multcomp package), multiple comparisons
    tuk <- cld(glht(fit, alternative = 'two.sided', linfct = mcp(class = 'Tukey')), sig = p, decreasing = TRUE)
    
    # Extract significance letters "abc" from cld() output
    sig <- data.frame(tuk$mcletters$Letters, stringsAsFactors = FALSE)
    names(sig) <- 'sig'
    sig$class <- rownames(sig)
    sig$var <- j
    
    dat_new <- dat %>% group_by(class) %>% mutate(mean = mean(var), sd = sd(var))
    dat_new2 <- unique(dat_new[, c(1, 3, 4)])
    dat_new2$class <- as.character(dat_new2$class)
    dat_new2$group <- i
    
    dat_new2 <- merge(dat_new2, sig, by = 'class')
    dat_new2 <- dat_new2[c(4, 6, 1, 2, 3, 5)]
    
    data_list <- rbind(data_list, dat_new2)
  }
}

stat_anova
data_list
rz <- data_list

library(rstatix)
Nitrogen2 <- Nitrogen_z
Nitrogen_l <- melt(Nitrogen2, id.vars = c("sample", "group", "Niche"), variable.name = "Variable", value.name = "Value")

library(dplyr)
test <- data_list
test <- test[, c(2:6)]
colnames(test)[1] <- c("Variable")
colnames(test)[2] <- c("group")

new_df <- Nitrogen2 %>% group_by(group) %>% summarize_all(max)
Nitrogen_m <- melt(new_df, id.vars = c("group", "sample", "Niche"), variable.name = "Variable", value.name = "max")
test$max <- Nitrogen_m$max[match(paste(test$Variable, test$group, sep = "_"), paste(Nitrogen_m$Variable, Nitrogen_m$group, sep = "_"))]

new_df2 <- Nitrogen2 %>% summarize_all(max)
Nitrogen_x <- melt(new_df2, id.vars = c("group", "sample", "Niche"), variable.name = "Variable", value.name = "max")
test$x <- Nitrogen_x$max[match(test$Variable, Nitrogen_x$Variable)]

colors <- c("#88CCEE", "#44AA99", "#117733", "#999933",
            "#BFDDEE", "#78C1D6", "#459F80", "#C5C297")

rz_l <- Nitrogen_l
rz_test <- test

#### Combine boxplots ####
total_l <- rbind(ori_l, rz_l)
total_test <- rbind(ori_test, rz_test)

total_l$group <- factor(total_l$group, levels = c("N0.O", "N2.O", "N4.O", "N6.O",
                                                  "N0.Z", "N2.Z", "N4.Z", "N6.Z"))

total_test$group <- factor(total_test$group, levels = c("N0.O", "N2.O", "N4.O", "N6.O",
                                                        "N0.Z", "N2.Z", "N4.Z", "N6.Z"))

# Note: 'var' variable not defined in provided code, using values instead
total_l$Variable <- factor(total_l$Variable, levels = values)
total_test$Variable <- factor(total_test$Variable, levels = values)

colors <- rep(c("#3162A0", "#3FBDA7", "#F08080", "#F0C986"), 2)

library(dplyr)
total_test <- total_test %>%
  group_by(Variable) %>%
  mutate(ymax = max(max))

p_nitrogen <- ggplot(total_l, aes(group, Value)) +
  geom_rect(aes(xmin = "N0.Z", xmax = "N6.Z", ymin = -Inf, ymax = Inf), fill = "#ECECEC", alpha = 0.3) +
  geom_rect(aes(xmin = 4.5, xmax = 8.5, ymin = -Inf, ymax = Inf), fill = "#ECECEC", alpha = 0.3) +
  geom_boxplot(aes(group, Value, color = group)) +
  facet_wrap(~Variable, ncol = 4, scales = "free") +
  geom_jitter(aes(x = group, y = Value, color = group), width = 0.01) +
  geom_text(data = total_test, aes(x = group, label = sig, y = max + ymax * 0.2), size = 4) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(x = NULL, y = "") +
  ggtitle("") +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title = element_text(size = 14, hjust = 0.5, color = "#333333", face = "bold", margin = margin(b = 10)),
    axis.text.x = element_text(angle = 0, vjust = 0.5, face = "bold", color = "#333333"),
    plot.margin = margin(10, 10, 10, 10, unit = "pt"),
    legend.position = "none",
    text = element_text(family = "serif")
  )

p_nitrogen
ggsave("Fig.S7 env factor boxplot.pdf", plot = p_nitrogen, device = "pdf", dpi = 200, width = 11, height = 10)