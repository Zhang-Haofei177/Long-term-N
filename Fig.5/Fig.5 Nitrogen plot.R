#########Fig.5 plot##########################
rm(list=ls())
library(tidyr)
library(ComplexHeatmap)
library(circlize)

#########data preparation####################
setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.5/")
#####Nitrogen######
spearman_r <- read.csv("Gene Abundance_Nitrogen_Spearman r.csv",header = T,row.names = 1)
spearman_p <- read.csv("Gene Abundance_Nitrogen_Spearman p.csv",header = T,row.names = 1)

spearman_r_wide <- spearman_r %>%
  pivot_wider(
    names_from = Niche,
    values_from = Nitrogen
  )

spearman_p_wide <- spearman_p %>%
  pivot_wider(
    names_from = Niche,
    values_from = Nitrogen
  )

LMM <- read.csv("Nitrogen_LMM result.csv",header = T,row.names = 1)
LMM <- LMM[c("Nitrogen.mean","Nitrogen.P"),]
LMM2 <- as.data.frame(t(LMM))
LMM2$Genes <- rownames(LMM2)

df <- merge(LMM2[c(1,3)],spearman_r_wide,by = "Genes")
pvalue <- merge(LMM2[c(2,3)],spearman_p_wide,by = "Genes")

colnames(df)[2] <- c("LM_effect")
colnames(pvalue)[2] <- c("LM_effect")

#############plot########################
pvalue$LM_effect <- ifelse(pvalue$LM_effect < 0.001, "***", 
                           ifelse(pvalue$LM_effect < 0.01, "**", 
                                  ifelse(pvalue$LM_effect < 0.05, "*", "")))

df <- df[!grepl("Anammox_Bacteria",df$Genes),]
pvalue <- pvalue[!grepl("Anammox_Bacteria",pvalue$Genes),]

heatmapdata <- df[-2]
heatmapdata <- as.data.frame(heatmapdata)
rownames(heatmapdata) <- heatmapdata$Genes
heatmapdata$Genes <- NULL

barplot_p <- pvalue[1:2]
pvalue <- as.data.frame(pvalue[-2])
rownames(pvalue) <- pvalue$Genes
pvalue$Genes <- NULL

# Fill NA values in pvalue with empty strings
pvalue <- pvalue %>%
  mutate_all(~ ifelse(is.na(.), "", .))

# Set heatmap related parameters
col_fun = colorRamp2(c(-1, 0, 1), c("#4979b6", "white", "#d9352a")) # Set scale range and colors
barcol <- read.csv("Gene Functions.csv", header = T, row.names = 1) # Group bars by gene functions
barcol <- barcol[!grepl("Anammox_Bacteria",rownames(barcol)), ] # Remove Anammox_Bacteria

# Extract effect values as bardata column for bar plot
bardata <- as.data.frame(df[1:2])
rownames(bardata) <- bardata$Genes

# Sort bardata in ascending order by LM_effect
bardata <- bardata[order(bardata$LM_effect), ]

# Reorder all data according to bardata row order
# Get sorted order
gene_order <- rownames(bardata)

# Reorder barcol according to same order
barcol <- barcol[gene_order, ]
bar_color <- barcol$colors
names(bar_color) <- rownames(barcol)

# Reorder heatmapdata according to same order
heatmapdata <- heatmapdata[gene_order, ]

# Reorder pvalue according to same order
pvalue <- pvalue[gene_order, ]

# Reorder barplot_p according to same order
barplot_p <- barplot_p[match(gene_order, barplot_p$Genes), ]

# Prepare bar plot data
bardata$Genes <- NULL

# Prepare text annotation data
lm_effect_values <- barplot_p$LM_effect

# Create bar plot annotation
bar_anno <- anno_barplot(
  bardata,
  which = "row",  # Specify as row annotation
  width = unit(3, "cm"),  # Bar plot width
  bar_width = 0.8,        # Bar width
  outline = FALSE,
  gp = gpar(fill = bar_color, lwd = 1)  # Bar plot colors and border width
)

# Create text annotation
text_anno <- anno_text(
  lm_effect_values,  # Text values
  which = "row",  # Specify as row annotation
  location = 0,  # 0 means left side, 1 means right side
  just = "left",  # Text alignment
  gp = gpar(fontsize = 8),  # Font size
  width = unit(1.5, "cm")  # Text annotation width
)

# Combine bar plot and text annotations
right_anno <- rowAnnotation(
  'Nitrogen effect' = bar_anno,
  'Value' = text_anno,
  annotation_width = unit(c(3, 1.5), "cm")  # Set width for both annotations
)

# Column annotations
group <- colnames(df)[3:6]
group_df <- data.frame(Niche = group)

topanno <- HeatmapAnnotation(
  Niche = group,  # Directly use named character vector
  col = list(
    Niche = c(
      'RZS' = '#7e7098',
      'RS' = '#416eb2',
      'RE' = '#82a29a',
      'NE' = '#e67f7e'
    )
  ),
  border = TRUE,
  show_annotation_name = TRUE
)

# Set row order - use already sorted order
x <- gene_order

# Draw heatmap
p <- Heatmap(
  heatmapdata,
  name = "Spearman's r",
  width = unit(4 * 0.8, "cm"),
  height = unit(13 * 0.8, "cm"),
  col = col_fun,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  row_names_side = 'left',  # Row annotations placed on left side
  top_annotation = topanno,
  right_annotation = right_anno,
  show_heatmap_legend = TRUE,
  show_column_names = TRUE,
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 1),  # Heatmap cell border color and width
  row_names_gp = gpar(fontsize = 10),
  border_gp = gpar(lwd = 1),
  row_order = x,  # Set row order
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%s", pvalue[i, j]), x, y, gp = gpar(fontsize = 10))
  }
)

p

pdf("Nitrogen_Spearman_LMM.pdf", width = 6, height = 8)
draw(p)
dev.off()