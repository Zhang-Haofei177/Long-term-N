setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.2andFig.3/iCAMP/Draw plot/")
library(reshape2)
library(ggplot2)

######################## AMF Process Importance Plot ###########
######################## Nitrogen ###########
a <- read.csv("AMF_1.ProcessImportance_EachGroup.csv", header = T)
group <- read.csv("group.csv", header = T, row.names = 1)
colnames(group)[3] <- "Group"
treat_name <- c("N0", "N2", "N4", "N6")
a2 <- a[a$Group %in% treat_name, ]
a2 <- a2[3:8]

a2_l <- melt(a2, id.vars = c("Group"), variable.name = "Assembly Process")
a2_l$`Assembly Process` <- factor(a2_l$`Assembly Process`, levels = c("HoS", "HeS", "DL", "HD", "DR"))
a2_l$Group <- factor(a2_l$Group, levels = c("N0", "N2", "N4", "N6"))
unique(a2_l$Group)

# Set color palette for each assembly process
color_assembly_Process <- c(
  "HeS" = "#DD6564",
  "HoS" = "#1B907E",
  "DL" = "#0C478D",
  "HD" = "#855396",
  "DR" = "#E2AB2E"
)

p <- ggplot(a2_l, aes(Group, value * 100, fill = `Assembly Process`)) +
  geom_rect(aes(xmin = 0, xmax = 5, ymin = -Inf, ymax = Inf), fill = "#f0f7ff", alpha = 1) +
  geom_col(position = 'stack', width = 0.6) +
  scale_fill_manual(values = color_assembly_Process) +
  labs(x = '', y = "Relative importance (%)", fill = 'Assembly Process') +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title = element_text(size = 14, hjust = 0.5, color = "#333333", face = "bold", margin = margin(b = 10)),
    axis.text.x = element_text(angle = 0, vjust = 0.5, face = "bold", color = "#333333"),
    plot.margin = margin(10, 10, 10, 10, unit = "pt"),
    legend.position = "none",
    text = element_text(family = "serif")
  )
p
ggsave(plot = p, filename = "AMF_ProcessImportance_EachNitrogen.pdf", device = "pdf", dpi = 300, width = 2, height = 3)

######################## Niche ###########
a <- read.csv("AMF_1.ProcessImportance_EachGroup.csv", header = T)
group <- read.csv("group.csv", header = T, row.names = 1)
colnames(group)[3] <- "Group"
treat_name <- c("Rootzone", "Rhizosphere", "Root", "Nodule")
a2 <- a[a$Group %in% treat_name, ]
a2 <- a2[3:8]

a2_l <- melt(a2, id.vars = c("Group"), variable.name = "Assembly Process")
a2_l$`Assembly Process` <- factor(a2_l$`Assembly Process`, levels = c("HoS", "HeS", "DL", "HD", "DR"))
a2_l$Group <- factor(a2_l$Group, levels = c("Rootzone", "Rhizosphere", "Root", "Nodule"),labels = c("RZS","RS","RE","NE"))
unique(a2_l$Group)

# Set color palette for each assembly process
color_assembly_Process <- c(
  "HeS" = "#DD6564",
  "HoS" = "#1B907E",
  "DL" = "#0C478D",
  "HD" = "#855396",
  "DR" = "#E2AB2E"
)

p <- ggplot(a2_l, aes(Group, value * 100, fill = `Assembly Process`)) +
  geom_rect(aes(xmin = 0, xmax = 5, ymin = -Inf, ymax = Inf), fill = "#f0f7ff", alpha = 1) +
  geom_col(position = 'stack', width = 0.6) +
  scale_fill_manual(values = color_assembly_Process) +
  labs(x = '', y = "Relative importance (%)", fill = 'Assembly Process') +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title = element_text(size = 14, hjust = 0.5, color = "#333333", face = "bold", margin = margin(b = 10)),
    axis.text.x = element_text(angle = 0, vjust = 0.5, face = "bold", color = "#333333"),
    plot.margin = margin(10, 10, 10, 10, unit = "pt"),
    legend.position = "none",
    text = element_text(family = "serif")
  )
p
ggsave(plot = p, filename = "AMF_ProcessImportance_EachNiche.pdf", device = "pdf", dpi = 300, width = 2, height = 3)

######################## Rhizobia Process Importance Plot ###########
######################## Nitrogen ###########
a <- read.csv("Rhi_1.ProcessImportance_EachGroup.csv", header = T)
group <- read.csv("group.csv", header = T, row.names = 1)
colnames(group)[3] <- "Group"
treat_name <- c("N0", "N2", "N4", "N6")
a2 <- a[a$Group %in% treat_name, ]
a2 <- a2[3:8]

a2_l <- melt(a2, id.vars = c("Group"), variable.name = "Assembly Process")
a2_l$`Assembly Process` <- factor(a2_l$`Assembly Process`, levels = c("HoS", "HeS", "DL", "HD", "DR"))
a2_l$Group <- factor(a2_l$Group, levels = c("N0", "N2", "N4", "N6"))
unique(a2_l$Group)

# Set color palette for each assembly process
color_assembly_Process <- c(
  "HeS" = "#DD6564",
  "HoS" = "#1B907E",
  "DL" = "#0C478D",
  "HD" = "#855396",
  "DR" = "#E2AB2E"
)

p <- ggplot(a2_l, aes(Group, value * 100, fill = `Assembly Process`)) +
  geom_rect(aes(xmin = 0, xmax = 5, ymin = -Inf, ymax = Inf), fill = "#f0f7ff", alpha = 1) +
  geom_col(position = 'stack', width = 0.6) +
  scale_fill_manual(values = color_assembly_Process) +
  labs(x = '', y = "Relative importance (%)", fill = 'Assembly Process') +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title = element_text(size = 14, hjust = 0.5, color = "#333333", face = "bold", margin = margin(b = 10)),
    axis.text.x = element_text(angle = 0, vjust = 0.5, face = "bold", color = "#333333"),
    plot.margin = margin(10, 10, 10, 10, unit = "pt"),
    legend.position = "none",
    text = element_text(family = "serif")
  )
p
ggsave(plot = p, filename = "Rhi_ProcessImportance_EachNitrogen.pdf", device = "pdf", dpi = 300, width = 2, height = 3)

######################## Niche ###########
a <- read.csv("Rhi_1.ProcessImportance_EachGroup.csv", header = T)
group <- read.csv("group.csv", header = T, row.names = 1)
colnames(group)[3] <- "Group"
treat_name <- c("Rootzone", "Rhizosphere", "Root", "Nodule")
a2 <- a[a$Group %in% treat_name, ]
a2 <- a2[3:8]

a2_l <- melt(a2, id.vars = c("Group"), variable.name = "Assembly Process")
a2_l$`Assembly Process` <- factor(a2_l$`Assembly Process`, levels = c("HoS", "HeS", "DL", "HD", "DR"))
a2_l$Group <- factor(a2_l$Group, levels = c("Rootzone", "Rhizosphere", "Root", "Nodule"),labels = c("RZS","RS","RE","NE"))
unique(a2_l$Group)

# Set color palette for each assembly process
color_assembly_Process <- c(
  "HeS" = "#DD6564",
  "HoS" = "#1B907E",
  "DL" = "#0C478D",
  "HD" = "#855396",
  "DR" = "#E2AB2E"
)

p <- ggplot(a2_l, aes(Group, value * 100, fill = `Assembly Process`)) +
  geom_rect(aes(xmin = 0, xmax = 5, ymin = -Inf, ymax = Inf), fill = "#f0f7ff", alpha = 1) +
  geom_col(position = 'stack', width = 0.6) +
  scale_fill_manual(values = color_assembly_Process) +
  labs(x = '', y = "Relative importance (%)", fill = 'Assembly Process') +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title = element_text(size = 14, hjust = 0.5, color = "#333333", face = "bold", margin = margin(b = 10)),
    axis.text.x = element_text(angle = 0, vjust = 0.5, face = "bold", color = "#333333"),
    plot.margin = margin(10, 10, 10, 10, unit = "pt"),
    legend.position = "none",
    text = element_text(family = "serif")
  )
p
ggsave(plot = p, filename = "Rhi_ProcessImportance_EachNiche.pdf", device = "pdf", dpi = 300, width = 2, height = 3)

######################### Stacked Bar Charts #####################################
############### AMF Stacked Bar Charts ##################################################
setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.2andFig.3/iCAMP/Draw plot/")
a <- read.csv("AMF_1.ProcessImportance_EachBin_EachGroup.csv", header = T)
niche_nitrogen <- subset(a, Group %in% c("Rootzone", "Rhizosphere", "Root", "Nodule", "N0", "N2", "N4", "N6") & Index %in% c("DominantProcess", "DominantProcessImportance", "DominantProcessPvalue"))
niche_nitrogen$Method <- NULL
niche_nitrogen$GroupBasedOn <- NULL
niche_nitrogen$Index <- gsub("DominantProcess", "Category",
                             gsub("DominantProcessImportance", "Percentage",
                                  gsub("DominantProcessPvalue", "p",
                                       niche_nitrogen$Index)))
new_rownames <- paste(niche_nitrogen$Group, niche_nitrogen$Index, sep = "_")
rownames(niche_nitrogen) <- new_rownames
niche_nitrogen <- niche_nitrogen[, !(colnames(niche_nitrogen) %in% c("Group", "Index"))]
niche_t <- t(niche_nitrogen[grepl("Rootzone|Rhizosphere|Root|Nodule", rownames(niche_nitrogen)), ])
nitrogen_t <- t(niche_nitrogen[grepl("N0|N2|N4|N6", rownames(niche_nitrogen)), ])
niche_t <- as.data.frame(niche_t)
nitrogen_t <- as.data.frame(nitrogen_t)
niche_RZS <- niche_t[grepl("Rootzone", colnames(niche_t))]
niche_RZS$Rootzone_p <- NULL
niche_RZS$Niche <- "RZS"
niche_RZS$bin <- rownames(niche_RZS)
colnames(niche_RZS) <- c("Category", "Percentage", "Niche", "bin")

niche_RS <- niche_t[grepl("Rhizosphere", colnames(niche_t))]
niche_RS$Rhizosphere_p <- NULL
niche_RS$Niche <- "RS"
niche_RS$bin <- rownames(niche_RS)
colnames(niche_RS) <- c("Category", "Percentage", "Niche", "bin")

niche_RE <- niche_t[grepl("Root", colnames(niche_t))]
niche_RE <- niche_RE[!grepl("Rootzone", colnames(niche_RE))]
niche_RE$Root_p <- NULL
niche_RE$Niche <- "RE"
niche_RE$bin <- rownames(niche_RE)
colnames(niche_RE) <- c("Category", "Percentage", "Niche", "bin")

niche_NE <- niche_t[grepl("Nodule", colnames(niche_t))]
niche_NE$Nodule_p <- NULL
niche_NE$Niche <- "NE"
niche_NE$bin <- rownames(niche_NE)
colnames(niche_NE) <- c("Category", "Percentage", "Niche", "bin")
niche_total <- rbind(niche_RZS, niche_RS, niche_RE, niche_NE)
niche_total <- niche_total[, c("bin", "Category", "Percentage", "Niche")]
rownames(niche_total) <- NULL

nitrogen_N0 <- nitrogen_t[grepl("N0", colnames(nitrogen_t))]
nitrogen_N0$N0_p <- NULL
nitrogen_N0$nitrogen <- "N0"
nitrogen_N0$bin <- rownames(nitrogen_N0)
colnames(nitrogen_N0) <- c("Category", "Percentage", "nitrogen", "bin")

nitrogen_N2 <- nitrogen_t[grepl("N2", colnames(nitrogen_t))]
nitrogen_N2$N2_p <- NULL
nitrogen_N2$nitrogen <- "N2"
nitrogen_N2$bin <- rownames(nitrogen_N2)
colnames(nitrogen_N2) <- c("Category", "Percentage", "nitrogen", "bin")

nitrogen_N4 <- nitrogen_t[grepl("N4", colnames(nitrogen_t))]
nitrogen_N4$N4_p <- NULL
nitrogen_N4$nitrogen <- "N4"
nitrogen_N4$bin <- rownames(nitrogen_N4)
colnames(nitrogen_N4) <- c("Category", "Percentage", "nitrogen", "bin")

nitrogen_N6 <- nitrogen_t[grepl("N6", colnames(nitrogen_t))]
nitrogen_N6$N6_p <- NULL
nitrogen_N6$nitrogen <- "N6"
nitrogen_N6$bin <- rownames(nitrogen_N6)
colnames(nitrogen_N6) <- c("Category", "Percentage", "nitrogen", "bin")
nitrogen_total <- rbind(nitrogen_N0, nitrogen_N2, nitrogen_N4, nitrogen_N6)
nitrogen_total <- nitrogen_total[, c("bin", "Category", "Percentage", "nitrogen")]
rownames(nitrogen_total) <- NULL

########################## AMF_Niche ########################################################
df <- niche_total
df$bin <- sub("^b", "B", df$bin)
colnames(df)[1] <- "Bin"
niche <- df[1:4]
tax <- read.csv("AMF_1.Bin_TopTaxon.csv", header = T)
tax <- tax[c(1:2, 10)]

niche_tax <- merge(niche, tax, by = "Bin")

df <- niche_tax

# Define niche order: RZS -> RS -> RE -> NE (increasing host selection)
niche_order <- c("RZS", "RS", "RE", "NE")

# 1. Get assembly process sequences for each bin across four niches (only complete data)
bin_sequences <- df %>%
  filter(Niche %in% niche_order) %>%
  # Ensure each bin has data for all four niches
  group_by(Bin) %>%
  filter(n_distinct(Niche) == 4) %>%  # Keep only bins with data in all four niches
  arrange(Bin, factor(Niche, levels = niche_order)) %>%
  summarise(
    # Assembly process sequence, e.g., "DR->DL->DL->DR"
    process_sequence = paste(Category, collapse = "->"),
    # Get genus annotation
    TopTaxon.Genus = first(TopTaxon.Genus),
    # Get relative abundance of the bin (each bin has only one abundance value)
    BinRA = first(BinRA)
  ) %>%
  ungroup()

# 2. Statistics by process sequence type and genus
sequence_genus_summary <- bin_sequences %>%
  group_by(process_sequence, TopTaxon.Genus) %>%
  summarise(
    bin_count = n(),
    total_RA = sum(BinRA, na.rm = TRUE),  # Add abundance only once per bin
    bins = paste(Bin, collapse = ", "),
    .groups = 'drop'
  ) %>%
  # Sort by sequence type and abundance
  arrange(process_sequence, desc(total_RA))

# 3. Prepare data for stacked chart: select main genera, merge others as "Other"
# Calculate total abundance per genus
genus_total_RA <- sequence_genus_summary %>%
  group_by(TopTaxon.Genus) %>%
  summarise(total_genus_RA = sum(total_RA)) %>%
  arrange(desc(total_genus_RA))

# Select top N main genera (e.g., top 10)
top_N <- 5
top_genera <- head(genus_total_RA$TopTaxon.Genus, top_N)

# Merge non-main genera as "Other"
stacked_data <- sequence_genus_summary %>%
  mutate(
    Genus_Group = ifelse(TopTaxon.Genus %in% top_genera, TopTaxon.Genus, "Other")
  ) %>%
  group_by(process_sequence, Genus_Group) %>%
  summarise(
    total_RA = sum(total_RA),
    bin_count = sum(bin_count),
    .groups = 'drop'
  )

# 4. Calculate total abundance per sequence type (for sorting)
sequence_totals <- stacked_data %>%
  group_by(process_sequence) %>%
  summarise(total_sequence_RA = sum(total_RA)) %>%
  arrange(desc(total_sequence_RA))

# Order sequence types by total abundance
stacked_data <- stacked_data %>%
  mutate(
    process_sequence = factor(process_sequence,
                              levels = sequence_totals$process_sequence)
  )

# 5. Create stacked bar chart
library(ggplot2)
library(viridis)
library(RColorBrewer)

# Choose color scheme - assign colors to main genera
n_colors <- length(unique(stacked_data$Genus_Group))
if (n_colors <= 12) {
  # Use Set3 palette, suitable for categorical data
  color_palette <- brewer.pal(n_colors, "Set3")
} else {
  # If not enough colors, use viridis
  color_palette <- viridis(n_colors)
}

# Create stacked bar chart
stacked_plot <- ggplot(stacked_data,
                       aes(x = process_sequence, y = total_RA, fill = Genus_Group)) +
  geom_col(alpha = 0.9, color = "black", linewidth = 0.3) +
  # Add total bin count above bars
  geom_text(data = sequence_totals,
            aes(x = process_sequence, y = total_sequence_RA,
                label = paste0("n=", sapply(process_sequence, function(seq) {
                  sum(stacked_data$bin_count[stacked_data$process_sequence == seq])
                }))),
            vjust = -0.5, size = 3.0, fontface = "bold", family = "Arial",
            inherit.aes = FALSE) +  # Add inherit.aes = FALSE
  # Add total abundance value inside bars
  geom_text(data = sequence_totals,
            aes(x = process_sequence, y = total_sequence_RA,
                label = sprintf("%.4f", total_sequence_RA)),
            vjust = 1.5, size = 2.5, color = "white", fontface = "bold", family = "Arial",
            inherit.aes = FALSE) +  # Add inherit.aes = FALSE
  labs(
    title = "Community Assembly Process Sequences by Genus",
    subtitle = "RZS → RS → RE → NE (Stacked by Genus)",
    x = "Assembly Process Sequence",
    y = "Cumulative Relative Abundance",
    fill = "Genus"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold", color = "black", margin = margin(t = 10)),
    axis.title.y = element_text(size = 11, face = "bold", color = "black", margin = margin(r = 10)),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "black", margin = margin(b = 10)),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(15, 15, 15, 15),
    plot.background = element_rect(color = "black", linewidth = 1, fill = "white")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = color_palette)

print(stacked_plot)

# 6. If there are too many sequence types, show only top N
top_sequences_count <- 15
if (nrow(sequence_totals) > top_sequences_count) {
  top_sequences <- head(sequence_totals$process_sequence, top_sequences_count)
  stacked_data_top <- stacked_data %>%
    filter(process_sequence %in% top_sequences)
  
  # Pre-calculate bin count per sequence type
  sequence_totals_top <- sequence_totals %>%
    filter(process_sequence %in% top_sequences) %>%
    mutate(
      bin_count = sapply(process_sequence, function(seq) {
        sum(stacked_data_top$bin_count[stacked_data_top$process_sequence == seq])
      })
    )
  stacked_plot_top <- ggplot(stacked_data_top,
                             aes(x = process_sequence, y = total_RA, fill = Genus_Group)) +
    geom_col(alpha = 0.9, color = "black", linewidth = 0.3) +
    geom_text(data = sequence_totals_top,
              aes(x = process_sequence, y = total_sequence_RA,
                  label = paste0("n=", bin_count)),
              vjust = -0.5, size = 3.0, fontface = "bold", family = "Arial",
              inherit.aes = FALSE) +  # Add inherit.aes = FALSE
    geom_text(data = sequence_totals_top,
              aes(x = process_sequence, y = total_sequence_RA,
                  label = sprintf("%.4f", total_sequence_RA)),
              vjust = 1.5, size = 2.5, color = "white", fontface = "bold", family = "Arial",
              inherit.aes = FALSE) +  # Add inherit.aes = FALSE
    labs(
      title = "Top Community Assembly Process Sequences by Genus",
      subtitle = paste("Top", top_sequences_count, "sequences by abundance"),
      x = "Assembly Process Sequence",
      y = "Cumulative Relative Abundance",
      fill = "Genus"
    ) +
    theme_classic() +
    theme(
      text = element_text(family = "Arial", color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.title.x = element_text(size = 11, face = "bold", color = "black", margin = margin(t = 10)),
      axis.title.y = element_text(size = 11, face = "bold", color = "black", margin = margin(r = 10)),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5, color = "black"),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "black", margin = margin(b = 10)),
      legend.position = "right",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.margin = margin(15, 15, 15, 15),
      plot.background = element_rect(color = "black", linewidth = 1, fill = "white")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    scale_fill_manual(values = color_palette)
  
  print(stacked_plot_top)
}

# 7. Save results
if (exists("stacked_plot_top")) {
  ggsave("AMF_Niche_Stacked_Top_Assembly_Process_Sequences_by_Genus.png", stacked_plot_top,
         width = 10, height = 5, dpi = 300)
  ggsave("AMF_Niche_Stacked_Top_Assembly_Process_Sequences_by_Genus.pdf", stacked_plot_top,
         width = 10, height = 5, device = cairo_pdf)
}

################# Pie chart ##########################
################# AMF Niche Pie chart ################
if (exists("sequence_totals") && exists("stacked_data")) {
  total_bins_all <- sum(stacked_data$bin_count, na.rm = TRUE)
  stable_sequences <- sequence_totals %>%
    mutate(
      processes = strsplit(process_sequence, "->"),
      is_stable = sapply(processes, function(x) length(unique(x)) == 1)
    ) %>%
    filter(is_stable)
  total_bin_stable <- 0
  stable_seq_list <- list()
  if (nrow(stable_sequences) > 0) {
    for (i in 1:nrow(stable_sequences)) {
      seq_type <- stable_sequences$process_sequence[i]
      seq_data <- stacked_data %>% filter(process_sequence == seq_type)
      bin_count_seq <- sum(seq_data$bin_count, na.rm = TRUE)
      total_bin_stable <- total_bin_stable + bin_count_seq
      
      stable_seq_list[[i]] <- data.frame(
        sequence = seq_type,
        type = "Stable",
        bin_count = bin_count_seq,
        relative_abundance = stable_sequences$total_sequence_RA[i]
      )
    }
  }
  changed_sequences <- sequence_totals %>%
    mutate(
      processes = strsplit(process_sequence, "->"),
      is_stable = sapply(processes, function(x) length(unique(x)) == 1)
    ) %>%
    filter(!is_stable)
  total_bin_changed <- 0
  changed_seq_list <- list()
  
  if (nrow(changed_sequences) > 0) {
    for (i in 1:nrow(changed_sequences)) {
      seq_type <- changed_sequences$process_sequence[i]
      seq_data <- stacked_data %>% filter(process_sequence == seq_type)
      bin_count_seq <- sum(seq_data$bin_count, na.rm = TRUE)
      total_bin_changed <- total_bin_changed + bin_count_seq
      
      changed_seq_list[[i]] <- data.frame(
        sequence = seq_type,
        type = "Changed",
        bin_count = bin_count_seq,
        relative_abundance = changed_sequences$total_sequence_RA[i]
      )
    }
  }
  
  all_sequences_info <- do.call(rbind, c(stable_seq_list, changed_seq_list))
  
  total_stable <- sum(all_sequences_info$bin_count[all_sequences_info$type == "Stable"], na.rm = TRUE)
  total_changed <- sum(all_sequences_info$bin_count[all_sequences_info$type == "Changed"], na.rm = TRUE)
  
  total_all <- total_stable + total_changed
  stable_ratio <- total_stable / total_all * 100
  changed_ratio <- total_changed / total_all * 100
  
  piechart_data <- data.frame(
    type = c("Stable", "Changed"),
    bin_count = c(total_stable, total_changed),
    percentage = c(stable_ratio, changed_ratio),
    label = c(
      paste0("Stable\n", total_stable, " bins\n", round(stable_ratio, 1), "%"),
      paste0("Changed\n", total_changed, " bins\n", round(changed_ratio, 1), "%")
    )
  )
  
  piechart_details <- list(
    summary = piechart_data,
    sequences = all_sequences_info,
    stats = data.frame(
      metric = c(
        "Total bins (stable)",
        "Total bins (changed)",
        "Total bins (all)",
        "Percentage (stable)",
        "Percentage (changed)"
      ),
      value = c(
        total_stable,
        total_changed,
        total_all,
        paste0(round(stable_ratio, 2), "%"),
        paste0(round(changed_ratio, 2), "%")
      )
    )
  )
  
} else {
  cat("Sequence type data not available\n\n")
  piechart_data <- NULL
  piechart_details <- NULL
}

if (!is.null(piechart_data)) {
  library(ggplot2)
  
  pie_chart <- ggplot(piechart_data, aes(x = "", y = bin_count, fill = type)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5),
              size = 5, color = "white", fontface = "bold") +
    scale_fill_manual(values = c("#4E79A7", "#F28E2B")) +
    labs(title = "Stable vs Changed bins",
         subtitle = paste("Based on", sum(piechart_data$bin_count), "total bins")
    ) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          legend.position = "right")
  print(pie_chart)
  # Save pie chart
  ggsave("AMF_Niche_piechart.pdf", pie_chart, width = 8, height = 8, dpi = 300)
}

##################### AMF_Nitrogen ##########################
df <- nitrogen_total
colnames(df) <- c("Bin", "Category", "Percentage", "Nitrogen")
df$Bin <- sub("^b", "B", df$Bin)
Nitrogen <- df
tax <- read.csv("AMF_1.Bin_TopTaxon.csv", header = T)
tax <- tax[c(1:2, 10)]
Nitrogen_tax <- merge(Nitrogen, tax, by = "Bin")

df <- Nitrogen_tax

# Define nitrogen order: N0 -> N2 -> N4 -> N6 (increasing nitrogen application)
Nitrogen_order <- c("N0", "N2", "N4", "N6")

# 1. Get assembly process sequences for each bin across four nitrogen levels (only complete data)
bin_sequences <- df %>%
  filter(Nitrogen %in% Nitrogen_order) %>%
  # Ensure each bin has data for all four nitrogen levels
  group_by(Bin) %>%
  filter(n_distinct(Nitrogen) == 4) %>%  # Keep only bins with data in all four nitrogen levels
  arrange(Bin, factor(Nitrogen, levels = Nitrogen_order)) %>%
  summarise(
    # Assembly process sequence, e.g., "DR->DL->DL->DR"
    process_sequence = paste(Category, collapse = "->"),
    # Get genus annotation
    TopTaxon.Genus = first(TopTaxon.Genus),
    # Get relative abundance of the bin (each bin has only one abundance value)
    BinRA = first(BinRA)
  ) %>%
  ungroup()

# 2. Statistics by process sequence type and genus
sequence_genus_summary <- bin_sequences %>%
  group_by(process_sequence, TopTaxon.Genus) %>%
  summarise(
    bin_count = n(),
    total_RA = sum(BinRA, na.rm = TRUE),  # Add abundance only once per bin
    bins = paste(Bin, collapse = ", "),
    .groups = 'drop'
  ) %>%
  # Sort by sequence type and abundance
  arrange(process_sequence, desc(total_RA))

# 3. Prepare data for stacked chart: select main genera, merge others as "Other"
# Calculate total abundance per genus
genus_total_RA <- sequence_genus_summary %>%
  group_by(TopTaxon.Genus) %>%
  summarise(total_genus_RA = sum(total_RA)) %>%
  arrange(desc(total_genus_RA))

# Select top N main genera (e.g., top 10)
top_N <- 5
top_genera <- head(genus_total_RA$TopTaxon.Genus, top_N)

# Merge non-main genera as "Other"
stacked_data <- sequence_genus_summary %>%
  mutate(
    Genus_Group = ifelse(TopTaxon.Genus %in% top_genera, TopTaxon.Genus, "Other")
  ) %>%
  group_by(process_sequence, Genus_Group) %>%
  summarise(
    total_RA = sum(total_RA),
    bin_count = sum(bin_count),
    .groups = 'drop'
  )

# 4. Calculate total abundance per sequence type (for sorting)
sequence_totals <- stacked_data %>%
  group_by(process_sequence) %>%
  summarise(total_sequence_RA = sum(total_RA)) %>%
  arrange(desc(total_sequence_RA))

# Order sequence types by total abundance
stacked_data <- stacked_data %>%
  mutate(
    process_sequence = factor(process_sequence,
                              levels = sequence_totals$process_sequence)
  )

# 5. Create stacked bar chart
library(ggplot2)
library(viridis)
library(RColorBrewer)

# Choose color scheme - assign colors to main genera
n_colors <- length(unique(stacked_data$Genus_Group))
if (n_colors <= 12) {
  # Use Set3 palette, suitable for categorical data
  color_palette <- brewer.pal(n_colors, "Set3")
} else {
  # If not enough colors, use viridis
  color_palette <- viridis(n_colors)
}

# Create stacked bar chart
stacked_plot <- ggplot(stacked_data,
                       aes(x = process_sequence, y = total_RA, fill = Genus_Group)) +
  geom_col(alpha = 0.9, color = "black", linewidth = 0.3) +
  # Add total bin count above bars
  geom_text(data = sequence_totals,
            aes(x = process_sequence, y = total_sequence_RA,
                label = paste0("n=", sapply(process_sequence, function(seq) {
                  sum(stacked_data$bin_count[stacked_data$process_sequence == seq])
                }))),
            vjust = -0.5, size = 3.0, fontface = "bold", family = "Arial",
            inherit.aes = FALSE) +  # Add inherit.aes = FALSE
  # Add total abundance value inside bars
  geom_text(data = sequence_totals,
            aes(x = process_sequence, y = total_sequence_RA,
                label = sprintf("%.4f", total_sequence_RA)),
            vjust = 1.5, size = 2.5, color = "white", fontface = "bold", family = "Arial",
            inherit.aes = FALSE) +  # Add inherit.aes = FALSE
  labs(
    title = "Community Assembly Process Sequences by Genus",
    subtitle = "N0 → N2 → N4 → N6 (Stacked by Genus)",
    x = "Assembly Process Sequence",
    y = "Cumulative Relative Abundance",
    fill = "Genus"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold", color = "black", margin = margin(t = 10)),
    axis.title.y = element_text(size = 11, face = "bold", color = "black", margin = margin(r = 10)),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "black", margin = margin(b = 10)),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(15, 15, 15, 15),
    plot.background = element_rect(color = "black", linewidth = 1, fill = "white")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = color_palette)

print(stacked_plot)

# 6. If there are too many sequence types, show only top N
top_sequences_count <- 15
if (nrow(sequence_totals) > top_sequences_count) {
  top_sequences <- head(sequence_totals$process_sequence, top_sequences_count)
  stacked_data_top <- stacked_data %>%
    filter(process_sequence %in% top_sequences)
  
  # Pre-calculate bin count per sequence type
  sequence_totals_top <- sequence_totals %>%
    filter(process_sequence %in% top_sequences) %>%
    mutate(
      bin_count = sapply(process_sequence, function(seq) {
        sum(stacked_data_top$bin_count[stacked_data_top$process_sequence == seq])
      })
    )
  
  stacked_plot_top <- ggplot(stacked_data_top,
                             aes(x = process_sequence, y = total_RA, fill = Genus_Group)) +
    geom_col(alpha = 0.9, color = "black", linewidth = 0.3) +
    geom_text(data = sequence_totals_top,
              aes(x = process_sequence, y = total_sequence_RA,
                  label = paste0("n=", bin_count)),
              vjust = -0.5, size = 3.0, fontface = "bold", family = "Arial",
              inherit.aes = FALSE) +  # Add inherit.aes = FALSE
    geom_text(data = sequence_totals_top,
              aes(x = process_sequence, y = total_sequence_RA,
                  label = sprintf("%.4f", total_sequence_RA)),
              vjust = 1.5, size = 2.5, color = "white", fontface = "bold", family = "Arial",
              inherit.aes = FALSE) +  # Add inherit.aes = FALSE
    labs(
      title = "Top Community Assembly Process Sequences by Genus",
      subtitle = paste("Top", top_sequences_count, "sequences by abundance"),
      x = "Assembly Process Sequence",
      y = "Cumulative Relative Abundance",
      fill = "Genus"
    ) +
    theme_classic() +
    theme(
      text = element_text(family = "Arial", color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.title.x = element_text(size = 11, face = "bold", color = "black", margin = margin(t = 10)),
      axis.title.y = element_text(size = 11, face = "bold", color = "black", margin = margin(r = 10)),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5, color = "black"),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "black", margin = margin(b = 10)),
      legend.position = "right",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.margin = margin(15, 15, 15, 15),
      plot.background = element_rect(color = "black", linewidth = 1, fill = "white")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    scale_fill_manual(values = color_palette)
  
  print(stacked_plot_top)
}

# 7. Save results
# Save both PNG and PDF formats
if (exists("stacked_plot_top")) {
  ggsave("AMF_Nitrogen_Stacked_Top_Assembly_Process_Sequences_by_Genus.png", stacked_plot_top,
         width = 10, height = 5, dpi = 300)
  ggsave("AMF_Nitrogen_Stacked_Top_Assembly_Process_Sequences_by_Genus.pdf", stacked_plot_top,
         width = 10, height = 5, device = cairo_pdf)
}

################# AMF Nitrogen Pie chart ######################
if (exists("sequence_totals") && exists("stacked_data")) {
  total_bins_all <- sum(stacked_data$bin_count, na.rm = TRUE)
  stable_sequences <- sequence_totals %>%
    mutate(
      processes = strsplit(process_sequence, "->"),
      is_stable = sapply(processes, function(x) length(unique(x)) == 1)
    ) %>%
    filter(is_stable)
  total_bin_stable <- 0
  stable_seq_list <- list()
  if (nrow(stable_sequences) > 0) {
    for (i in 1:nrow(stable_sequences)) {
      seq_type <- stable_sequences$process_sequence[i]
      seq_data <- stacked_data %>% filter(process_sequence == seq_type)
      bin_count_seq <- sum(seq_data$bin_count, na.rm = TRUE)
      total_bin_stable <- total_bin_stable + bin_count_seq
      
      stable_seq_list[[i]] <- data.frame(
        sequence = seq_type,
        type = "Stable",
        bin_count = bin_count_seq,
        relative_abundance = stable_sequences$total_sequence_RA[i]
      )
    }
  }
  changed_sequences <- sequence_totals %>%
    mutate(
      processes = strsplit(process_sequence, "->"),
      is_stable = sapply(processes, function(x) length(unique(x)) == 1)
    ) %>%
    filter(!is_stable)
  total_bin_changed <- 0
  changed_seq_list <- list()
  
  if (nrow(changed_sequences) > 0) {
    for (i in 1:nrow(changed_sequences)) {
      seq_type <- changed_sequences$process_sequence[i]
      seq_data <- stacked_data %>% filter(process_sequence == seq_type)
      bin_count_seq <- sum(seq_data$bin_count, na.rm = TRUE)
      total_bin_changed <- total_bin_changed + bin_count_seq
      
      changed_seq_list[[i]] <- data.frame(
        sequence = seq_type,
        type = "Changed",
        bin_count = bin_count_seq,
        relative_abundance = changed_sequences$total_sequence_RA[i]
      )
    }
  }
  
  all_sequences_info <- do.call(rbind, c(stable_seq_list, changed_seq_list))
  
  total_stable <- sum(all_sequences_info$bin_count[all_sequences_info$type == "Stable"], na.rm = TRUE)
  total_changed <- sum(all_sequences_info$bin_count[all_sequences_info$type == "Changed"], na.rm = TRUE)
  
  total_all <- total_stable + total_changed
  stable_ratio <- total_stable / total_all * 100
  changed_ratio <- total_changed / total_all * 100
  
  piechart_data <- data.frame(
    type = c("Stable", "Changed"),
    bin_count = c(total_stable, total_changed),
    percentage = c(stable_ratio, changed_ratio),
    label = c(
      paste0("Stable\n", total_stable, " bins\n", round(stable_ratio, 1), "%"),
      paste0("Changed\n", total_changed, " bins\n", round(changed_ratio, 1), "%")
    )
  )
  
  piechart_details <- list(
    summary = piechart_data,
    sequences = all_sequences_info,
    stats = data.frame(
      metric = c(
        "Total bins (stable)",
        "Total bins (changed)",
        "Total bins (all)",
        "Percentage (stable)",
        "Percentage (changed)"
      ),
      value = c(
        total_stable,
        total_changed,
        total_all,
        paste0(round(stable_ratio, 2), "%"),
        paste0(round(changed_ratio, 2), "%")
      )
    )
  )
  
} else {
  cat("Sequence type data not available\n\n")
  piechart_data <- NULL
  piechart_details <- NULL
}

if (!is.null(piechart_data)) {
  library(ggplot2)
  
  pie_chart <- ggplot(piechart_data, aes(x = "", y = bin_count, fill = type)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5),
              size = 5, color = "white", fontface = "bold") +
    scale_fill_manual(values = c("#4E79A7", "#F28E2B")) +
    labs(title = "Stable vs Changed bins",
         subtitle = paste("Based on", sum(piechart_data$bin_count), "total bins")
    ) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          legend.position = "right")
  print(pie_chart)
  # Save pie chart
  ggsave("AMF_Nitrogen_piechart.pdf", pie_chart, width = 8, height = 8, dpi = 300)
}

############### Rhizobia Stacked Bar Charts ##################################################
setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.2andFig.3/iCAMP/Draw plot/")
a <- read.csv("Rhi_1.ProcessImportance_EachBin_EachGroup.csv", header = T)
niche_nitrogen <- subset(a, Group %in% c("Rootzone", "Rhizosphere", "Root", "Nodule", "N0", "N2", "N4", "N6") & Index %in% c("DominantProcess", "DominantProcessImportance", "DominantProcessPvalue"))
niche_nitrogen$Method <- NULL
niche_nitrogen$GroupBasedOn <- NULL
niche_nitrogen$Index <- gsub("DominantProcess", "Category",
                             gsub("DominantProcessImportance", "Percentage",
                                  gsub("DominantProcessPvalue", "p",
                                       niche_nitrogen$Index)))
new_rownames <- paste(niche_nitrogen$Group, niche_nitrogen$Index, sep = "_")
rownames(niche_nitrogen) <- new_rownames
niche_nitrogen <- niche_nitrogen[, !(colnames(niche_nitrogen) %in% c("Group", "Index"))]
niche_t <- t(niche_nitrogen[grepl("Rootzone|Rhizosphere|Root|Nodule", rownames(niche_nitrogen)), ])
nitrogen_t <- t(niche_nitrogen[grepl("N0|N2|N4|N6", rownames(niche_nitrogen)), ])
niche_t <- as.data.frame(niche_t)
nitrogen_t <- as.data.frame(nitrogen_t)
niche_RZS <- niche_t[grepl("Rootzone", colnames(niche_t))]
niche_RZS$Rootzone_p <- NULL
niche_RZS$Niche <- "RZS"
niche_RZS$bin <- rownames(niche_RZS)
colnames(niche_RZS) <- c("Category", "Percentage", "Niche", "bin")

niche_RS <- niche_t[grepl("Rhizosphere", colnames(niche_t))]
niche_RS$Rhizosphere_p <- NULL
niche_RS$Niche <- "RS"
niche_RS$bin <- rownames(niche_RS)
colnames(niche_RS) <- c("Category", "Percentage", "Niche", "bin")

niche_RE <- niche_t[grepl("Root", colnames(niche_t))]
niche_RE <- niche_RE[!grepl("Rootzone", colnames(niche_RE))]
niche_RE$Root_p <- NULL
niche_RE$Niche <- "RE"
niche_RE$bin <- rownames(niche_RE)
colnames(niche_RE) <- c("Category", "Percentage", "Niche", "bin")

niche_NE <- niche_t[grepl("Nodule", colnames(niche_t))]
niche_NE$Nodule_p <- NULL
niche_NE$Niche <- "NE"
niche_NE$bin <- rownames(niche_NE)
colnames(niche_NE) <- c("Category", "Percentage", "Niche", "bin")
niche_total <- rbind(niche_RZS, niche_RS, niche_RE, niche_NE)
niche_total <- niche_total[, c("bin", "Category", "Percentage", "Niche")]
rownames(niche_total) <- NULL

nitrogen_N0 <- nitrogen_t[grepl("N0", colnames(nitrogen_t))]
nitrogen_N0$N0_p <- NULL
nitrogen_N0$nitrogen <- "N0"
nitrogen_N0$bin <- rownames(nitrogen_N0)
colnames(nitrogen_N0) <- c("Category", "Percentage", "nitrogen", "bin")

nitrogen_N2 <- nitrogen_t[grepl("N2", colnames(nitrogen_t))]
nitrogen_N2$N2_p <- NULL
nitrogen_N2$nitrogen <- "N2"
nitrogen_N2$bin <- rownames(nitrogen_N2)
colnames(nitrogen_N2) <- c("Category", "Percentage", "nitrogen", "bin")

nitrogen_N4 <- nitrogen_t[grepl("N4", colnames(nitrogen_t))]
nitrogen_N4$N4_p <- NULL
nitrogen_N4$nitrogen <- "N4"
nitrogen_N4$bin <- rownames(nitrogen_N4)
colnames(nitrogen_N4) <- c("Category", "Percentage", "nitrogen", "bin")

nitrogen_N6 <- nitrogen_t[grepl("N6", colnames(nitrogen_t))]
nitrogen_N6$N6_p <- NULL
nitrogen_N6$nitrogen <- "N6"
nitrogen_N6$bin <- rownames(nitrogen_N6)
colnames(nitrogen_N6) <- c("Category", "Percentage", "nitrogen", "bin")
nitrogen_total <- rbind(nitrogen_N0, nitrogen_N2, nitrogen_N4, nitrogen_N6)
nitrogen_total <- nitrogen_total[, c("bin", "Category", "Percentage", "nitrogen")]
rownames(nitrogen_total) <- NULL

########################## Rhizobia_Niche ########################################################
df <- niche_total
df$bin <- sub("^b", "B", df$bin)
colnames(df)[1] <- "Bin"
niche <- df[1:4]
tax <- read.csv("Rhi_1.Bin_TopTaxon.csv", header = T)
tax <- tax[c(1:2, 10)]

niche_tax <- merge(niche, tax, by = "Bin")

df <- niche_tax

# Define niche order: RZS -> RS -> RE -> NE (increasing host selection)
niche_order <- c("RZS", "RS", "RE", "NE")

# 1. Get assembly process sequences for each bin across four niches (only complete data)
bin_sequences <- df %>%
  filter(Niche %in% niche_order) %>%
  # Ensure each bin has data for all four niches
  group_by(Bin) %>%
  filter(n_distinct(Niche) == 4) %>%  # Keep only bins with data in all four niches
  arrange(Bin, factor(Niche, levels = niche_order)) %>%
  summarise(
    # Assembly process sequence, e.g., "DR->DL->DL->DR"
    process_sequence = paste(Category, collapse = "->"),
    # Get genus annotation
    TopTaxon.Genus = first(TopTaxon.Genus),
    # Get relative abundance of the bin (each bin has only one abundance value)
    BinRA = first(BinRA)
  ) %>%
  ungroup()

# 2. Statistics by process sequence type and genus
sequence_genus_summary <- bin_sequences %>%
  group_by(process_sequence, TopTaxon.Genus) %>%
  summarise(
    bin_count = n(),
    total_RA = sum(BinRA, na.rm = TRUE),  # Add abundance only once per bin
    bins = paste(Bin, collapse = ", "),
    .groups = 'drop'
  ) %>%
  # Sort by sequence type and abundance
  arrange(process_sequence, desc(total_RA))

# 3. Prepare data for stacked chart: select main genera, merge others as "Other"
# Calculate total abundance per genus
genus_total_RA <- sequence_genus_summary %>%
  group_by(TopTaxon.Genus) %>%
  summarise(total_genus_RA = sum(total_RA)) %>%
  arrange(desc(total_genus_RA))

# Select top N main genera
top_N <- 11
top_genera <- head(genus_total_RA$TopTaxon.Genus, top_N)
top_genera <- top_genera[!is.na(top_genera)]
# Merge non-main genera as "Other"
stacked_data <- sequence_genus_summary %>%
  mutate(
    Genus_Group = ifelse(TopTaxon.Genus %in% top_genera, TopTaxon.Genus, "Other")
  ) %>%
  group_by(process_sequence, Genus_Group) %>%
  summarise(
    total_RA = sum(total_RA),
    bin_count = sum(bin_count),
    .groups = 'drop'
  )

# 4. Calculate total abundance per sequence type (for sorting)
sequence_totals <- stacked_data %>%
  group_by(process_sequence) %>%
  summarise(total_sequence_RA = sum(total_RA)) %>%
  arrange(desc(total_sequence_RA))

# Order sequence types by total abundance
stacked_data <- stacked_data %>%
  mutate(
    process_sequence = factor(process_sequence,
                              levels = sequence_totals$process_sequence)
  )

# 5. Create stacked bar chart
library(ggplot2)
library(viridis)
library(RColorBrewer)

# Choose color scheme - assign colors to main genera
n_colors <- length(unique(stacked_data$Genus_Group))
if (n_colors <= 12) {
  # Use Set3 palette, suitable for categorical data
  color_palette <- brewer.pal(n_colors, "Set3")
} else {
  # If not enough colors, use viridis
  color_palette <- viridis(n_colors)
}

# Create stacked bar chart
stacked_plot <- ggplot(stacked_data,
                       aes(x = process_sequence, y = total_RA, fill = Genus_Group)) +
  geom_col(alpha = 0.9, color = "black", linewidth = 0.3) +
  # Add total bin count above bars
  geom_text(data = sequence_totals,
            aes(x = process_sequence, y = total_sequence_RA,
                label = paste0("n=", sapply(process_sequence, function(seq) {
                  sum(stacked_data$bin_count[stacked_data$process_sequence == seq])
                }))),
            vjust = -0.5, size = 3.0, fontface = "bold", family = "Arial",
            inherit.aes = FALSE) +  # Add inherit.aes = FALSE
  # Add total abundance value inside bars
  geom_text(data = sequence_totals,
            aes(x = process_sequence, y = total_sequence_RA,
                label = sprintf("%.4f", total_sequence_RA)),
            vjust = 1.5, size = 2.5, color = "white", fontface = "bold", family = "Arial",
            inherit.aes = FALSE) +  # Add inherit.aes = FALSE
  labs(
    title = "Community Assembly Process Sequences by Genus",
    subtitle = "RZS → RS → RE → NE (Stacked by Genus)",
    x = "Assembly Process Sequence",
    y = "Cumulative Relative Abundance",
    fill = "Genus"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold", color = "black", margin = margin(t = 10)),
    axis.title.y = element_text(size = 11, face = "bold", color = "black", margin = margin(r = 10)),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "black", margin = margin(b = 10)),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(15, 15, 15, 15),
    plot.background = element_rect(color = "black", linewidth = 1, fill = "white")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = color_palette)

print(stacked_plot)

# 6. If there are too many sequence types, show only top N
top_sequences_count <- 15
if (nrow(sequence_totals) > top_sequences_count) {
  top_sequences <- head(sequence_totals$process_sequence, top_sequences_count)
  stacked_data_top <- stacked_data %>%
    filter(process_sequence %in% top_sequences)
  
  # Pre-calculate bin count per sequence type
  sequence_totals_top <- sequence_totals %>%
    filter(process_sequence %in% top_sequences) %>%
    mutate(
      bin_count = sapply(process_sequence, function(seq) {
        sum(stacked_data_top$bin_count[stacked_data_top$process_sequence == seq])
      })
    )
  
  stacked_plot_top <- ggplot(stacked_data_top,
                             aes(x = process_sequence, y = total_RA, fill = Genus_Group)) +
    geom_col(alpha = 0.9, color = "black", linewidth = 0.3) +
    geom_text(data = sequence_totals_top,
              aes(x = process_sequence, y = total_sequence_RA,
                  label = paste0("n=", bin_count)),
              vjust = -0.5, size = 3.0, fontface = "bold", family = "Arial",
              inherit.aes = FALSE) +  # Add inherit.aes = FALSE
    geom_text(data = sequence_totals_top,
              aes(x = process_sequence, y = total_sequence_RA,
                  label = sprintf("%.4f", total_sequence_RA)),
              vjust = 1.5, size = 2.5, color = "white", fontface = "bold", family = "Arial",
              inherit.aes = FALSE) +  # Add inherit.aes = FALSE
    labs(
      title = "Top Community Assembly Process Sequences by Genus",
      subtitle = paste("Top", top_sequences_count, "sequences by abundance"),
      x = "Assembly Process Sequence",
      y = "Cumulative Relative Abundance",
      fill = "Genus"
    ) +
    theme_classic() +
    theme(
      text = element_text(family = "Arial", color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.title.x = element_text(size = 11, face = "bold", color = "black", margin = margin(t = 10)),
      axis.title.y = element_text(size = 11, face = "bold", color = "black", margin = margin(r = 10)),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5, color = "black"),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "black", margin = margin(b = 10)),
      legend.position = "right",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.margin = margin(15, 15, 15, 15),
      plot.background = element_rect(color = "black", linewidth = 1, fill = "white")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    scale_fill_manual(values = color_palette)
  
  print(stacked_plot_top)
}

# 7. Save results
# Save both PNG and PDF formats
if (exists("stacked_plot_top")) {
  ggsave("Rhizobia_Niche_Stacked_Top_Assembly_Process_Sequences_by_Genus.png", stacked_plot_top,
         width = 10, height = 5, dpi = 300)
  ggsave("Rhizobia_Niche_Stacked_Top_Assembly_Process_Sequences_by_Genus.pdf", stacked_plot_top,
         width = 10, height = 5, device = cairo_pdf)
}

################# Rhizobia Niche Pie chart ######################
if (exists("sequence_totals") && exists("stacked_data")) {
  total_bins_all <- sum(stacked_data$bin_count, na.rm = TRUE)
  stable_sequences <- sequence_totals %>%
    mutate(
      processes = strsplit(process_sequence, "->"),
      is_stable = sapply(processes, function(x) length(unique(x)) == 1)
    ) %>%
    filter(is_stable)
  total_bin_stable <- 0
  stable_seq_list <- list()
  if (nrow(stable_sequences) > 0) {
    for (i in 1:nrow(stable_sequences)) {
      seq_type <- stable_sequences$process_sequence[i]
      seq_data <- stacked_data %>% filter(process_sequence == seq_type)
      bin_count_seq <- sum(seq_data$bin_count, na.rm = TRUE)
      total_bin_stable <- total_bin_stable + bin_count_seq
      
      stable_seq_list[[i]] <- data.frame(
        sequence = seq_type,
        type = "Stable",
        bin_count = bin_count_seq,
        relative_abundance = stable_sequences$total_sequence_RA[i]
      )
    }
  }
  changed_sequences <- sequence_totals %>%
    mutate(
      processes = strsplit(process_sequence, "->"),
      is_stable = sapply(processes, function(x) length(unique(x)) == 1)
    ) %>%
    filter(!is_stable)
  total_bin_changed <- 0
  changed_seq_list <- list()
  
  if (nrow(changed_sequences) > 0) {
    for (i in 1:nrow(changed_sequences)) {
      seq_type <- changed_sequences$process_sequence[i]
      seq_data <- stacked_data %>% filter(process_sequence == seq_type)
      bin_count_seq <- sum(seq_data$bin_count, na.rm = TRUE)
      total_bin_changed <- total_bin_changed + bin_count_seq
      
      changed_seq_list[[i]] <- data.frame(
        sequence = seq_type,
        type = "Changed",
        bin_count = bin_count_seq,
        relative_abundance = changed_sequences$total_sequence_RA[i]
      )
    }
  }
  
  all_sequences_info <- do.call(rbind, c(stable_seq_list, changed_seq_list))
  
  total_stable <- sum(all_sequences_info$bin_count[all_sequences_info$type == "Stable"], na.rm = TRUE)
  total_changed <- sum(all_sequences_info$bin_count[all_sequences_info$type == "Changed"], na.rm = TRUE)
  
  total_all <- total_stable + total_changed
  stable_ratio <- total_stable / total_all * 100
  changed_ratio <- total_changed / total_all * 100
  
  piechart_data <- data.frame(
    type = c("Stable", "Changed"),
    bin_count = c(total_stable, total_changed),
    percentage = c(stable_ratio, changed_ratio),
    label = c(
      paste0("Stable\n", total_stable, " bins\n", round(stable_ratio, 1), "%"),
      paste0("Changed\n", total_changed, " bins\n", round(changed_ratio, 1), "%")
    )
  )
  
  piechart_details <- list(
    summary = piechart_data,
    sequences = all_sequences_info,
    stats = data.frame(
      metric = c(
        "Total bins (stable)",
        "Total bins (changed)",
        "Total bins (all)",
        "Percentage (stable)",
        "Percentage (changed)"
      ),
      value = c(
        total_stable,
        total_changed,
        total_all,
        paste0(round(stable_ratio, 2), "%"),
        paste0(round(changed_ratio, 2), "%")
      )
    )
  )
  
} else {
  cat("Sequence type data not available\n\n")
  piechart_data <- NULL
  piechart_details <- NULL
}

if (!is.null(piechart_data)) {
  library(ggplot2)
  
  pie_chart <- ggplot(piechart_data, aes(x = "", y = bin_count, fill = type)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5),
              size = 5, color = "white", fontface = "bold") +
    scale_fill_manual(values = c("#4E79A7", "#F28E2B")) +
    labs(title = "Stable vs Changed bins",
         subtitle = paste("Based on", sum(piechart_data$bin_count), "total bins")
    ) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          legend.position = "right")
  print(pie_chart)
  # Save pie chart
  ggsave("Rhizobia_Niche_piechart.pdf", pie_chart, width = 8, height = 8, dpi = 300)
}

##################### Rhizobia_Nitrogen ##########################
df <- nitrogen_total
colnames(df) <- c("Bin", "Category", "Percentage", "Nitrogen")
df$Bin <- sub("^b", "B", df$Bin)
Nitrogen <- df
tax <- read.csv("Rhi_1.Bin_TopTaxon.csv", header = T)
tax <- tax[c(1:2, 10)]
Nitrogen_tax <- merge(Nitrogen, tax, by = "Bin")

df <- Nitrogen_tax

# Define nitrogen order: N0 -> N2 -> N4 -> N6 (increasing nitrogen application)
Nitrogen_order <- c("N0", "N2", "N4", "N6")

# 1. Get assembly process sequences for each bin across four nitrogen levels (only complete data)
bin_sequences <- df %>%
  filter(Nitrogen %in% Nitrogen_order) %>%
  # Ensure each bin has data for all four nitrogen levels
  group_by(Bin) %>%
  filter(n_distinct(Nitrogen) == 4) %>%  # Keep only bins with data in all four nitrogen levels
  arrange(Bin, factor(Nitrogen, levels = Nitrogen_order)) %>%
  summarise(
    # Assembly process sequence, e.g., "DR->DL->DL->DR"
    process_sequence = paste(Category, collapse = "->"),
    # Get genus annotation
    TopTaxon.Genus = first(TopTaxon.Genus),
    # Get relative abundance of the bin (each bin has only one abundance value)
    BinRA = first(BinRA)
  ) %>%
  ungroup()

# 2. Statistics by process sequence type and genus
sequence_genus_summary <- bin_sequences %>%
  group_by(process_sequence, TopTaxon.Genus) %>%
  summarise(
    bin_count = n(),
    total_RA = sum(BinRA, na.rm = TRUE),  # Add abundance only once per bin
    bins = paste(Bin, collapse = ", "),
    .groups = 'drop'
  ) %>%
  # Sort by sequence type and abundance
  arrange(process_sequence, desc(total_RA))

# 3. Prepare data for stacked chart: select main genera, merge others as "Other"
# Calculate total abundance per genus
genus_total_RA <- sequence_genus_summary %>%
  group_by(TopTaxon.Genus) %>%
  summarise(total_genus_RA = sum(total_RA)) %>%
  arrange(desc(total_genus_RA))

# Select top N main genera (e.g., top 10)
top_N <- 11
top_genera <- head(genus_total_RA$TopTaxon.Genus, top_N)
top_genera <- top_genera[!is.na(top_genera)]
# Merge non-main genera as "Other"
stacked_data <- sequence_genus_summary %>%
  mutate(
    Genus_Group = ifelse(TopTaxon.Genus %in% top_genera, TopTaxon.Genus, "Other")
  ) %>%
  group_by(process_sequence, Genus_Group) %>%
  summarise(
    total_RA = sum(total_RA),
    bin_count = sum(bin_count),
    .groups = 'drop'
  )

# 4. Calculate total abundance per sequence type (for sorting)
sequence_totals <- stacked_data %>%
  group_by(process_sequence) %>%
  summarise(total_sequence_RA = sum(total_RA)) %>%
  arrange(desc(total_sequence_RA))

# Order sequence types by total abundance
stacked_data <- stacked_data %>%
  mutate(
    process_sequence = factor(process_sequence,
                              levels = sequence_totals$process_sequence)
  )

# 5. Create stacked bar chart
library(ggplot2)
library(viridis)
library(RColorBrewer)

# Choose color scheme - assign colors to main genera
n_colors <- length(unique(stacked_data$Genus_Group))
if (n_colors <= 12) {
  # Use Set3 palette, suitable for categorical data
  color_palette <- brewer.pal(n_colors, "Set3")
} else {
  # If not enough colors, use viridis
  color_palette <- viridis(n_colors)
}

# Create stacked bar chart
stacked_plot <- ggplot(stacked_data,
                       aes(x = process_sequence, y = total_RA, fill = Genus_Group)) +
  geom_col(alpha = 0.9, color = "black", linewidth = 0.3) +
  # Add total bin count above bars
  geom_text(data = sequence_totals,
            aes(x = process_sequence, y = total_sequence_RA,
                label = paste0("n=", sapply(process_sequence, function(seq) {
                  sum(stacked_data$bin_count[stacked_data$process_sequence == seq])
                }))),
            vjust = -0.5, size = 3.0, fontface = "bold", family = "Arial",
            inherit.aes = FALSE) +  # Add inherit.aes = FALSE
  # Add total abundance value inside bars
  geom_text(data = sequence_totals,
            aes(x = process_sequence, y = total_sequence_RA,
                label = sprintf("%.4f", total_sequence_RA)),
            vjust = 1.5, size = 2.5, color = "white", fontface = "bold", family = "Arial",
            inherit.aes = FALSE) +  # Add inherit.aes = FALSE
  labs(
    title = "Community Assembly Process Sequences by Genus",
    subtitle = "N0 → N2 → N4 → N6 (Stacked by Genus)",
    x = "Assembly Process Sequence",
    y = "Cumulative Relative Abundance",
    fill = "Genus"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.x = element_text(size = 11, face = "bold", color = "black", margin = margin(t = 10)),
    axis.title.y = element_text(size = 11, face = "bold", color = "black", margin = margin(r = 10)),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5, color = "black"),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "black", margin = margin(b = 10)),
    legend.position = "right",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(15, 15, 15, 15),
    plot.background = element_rect(color = "black", linewidth = 1, fill = "white")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = color_palette)

print(stacked_plot)

# 6. If there are too many sequence types, show only top N
top_sequences_count <- 15
if (nrow(sequence_totals) > top_sequences_count) {
  top_sequences <- head(sequence_totals$process_sequence, top_sequences_count)
  stacked_data_top <- stacked_data %>%
    filter(process_sequence %in% top_sequences)
  
  # Pre-calculate bin count per sequence type
  sequence_totals_top <- sequence_totals %>%
    filter(process_sequence %in% top_sequences) %>%
    mutate(
      bin_count = sapply(process_sequence, function(seq) {
        sum(stacked_data_top$bin_count[stacked_data_top$process_sequence == seq])
      })
    )
  
  stacked_plot_top <- ggplot(stacked_data_top,
                             aes(x = process_sequence, y = total_RA, fill = Genus_Group)) +
    geom_col(alpha = 0.9, color = "black", linewidth = 0.3) +
    geom_text(data = sequence_totals_top,
              aes(x = process_sequence, y = total_sequence_RA,
                  label = paste0("n=", bin_count)),
              vjust = -0.5, size = 3.0, fontface = "bold", family = "Arial",
              inherit.aes = FALSE) +  # Add inherit.aes = FALSE
    geom_text(data = sequence_totals_top,
              aes(x = process_sequence, y = total_sequence_RA,
                  label = sprintf("%.4f", total_sequence_RA)),
              vjust = 1.5, size = 2.5, color = "white", fontface = "bold", family = "Arial",
              inherit.aes = FALSE) +  # Add inherit.aes = FALSE
    labs(
      title = "Top Community Assembly Process Sequences by Genus",
      subtitle = paste("Top", top_sequences_count, "sequences by abundance"),
      x = "Assembly Process Sequence",
      y = "Cumulative Relative Abundance",
      fill = "Genus"
    ) +
    theme_classic() +
    theme(
      text = element_text(family = "Arial", color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.title.x = element_text(size = 11, face = "bold", color = "black", margin = margin(t = 10)),
      axis.title.y = element_text(size = 11, face = "bold", color = "black", margin = margin(r = 10)),
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5, color = "black"),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "black", margin = margin(b = 10)),
      legend.position = "right",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.margin = margin(15, 15, 15, 15),
      plot.background = element_rect(color = "black", linewidth = 1, fill = "white")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    scale_fill_manual(values = color_palette)
  
  print(stacked_plot_top)
}

# 7. Save results
# Save both PNG and PDF formats
if (exists("stacked_plot_top")) {
  ggsave("Rhizobia_Nitrogen_Stacked_Top_Assembly_Process_Sequences_by_Genus.png", stacked_plot_top,
         width = 10, height = 5, dpi = 300)
  ggsave("Rhizobia_Nitrogen_Stacked_Top_Assembly_Process_Sequences_by_Genus.pdf", stacked_plot_top,
         width = 10, height = 5, device = cairo_pdf)
}

################# Rhizobia Nitrogen Pie chart ######################
if (exists("sequence_totals") && exists("stacked_data")) {
  total_bins_all <- sum(stacked_data$bin_count, na.rm = TRUE)
  stable_sequences <- sequence_totals %>%
    mutate(
      processes = strsplit(process_sequence, "->"),
      is_stable = sapply(processes, function(x) length(unique(x)) == 1)
    ) %>%
    filter(is_stable)
  total_bin_stable <- 0
  stable_seq_list <- list()
  if (nrow(stable_sequences) > 0) {
    for (i in 1:nrow(stable_sequences)) {
      seq_type <- stable_sequences$process_sequence[i]
      seq_data <- stacked_data %>% filter(process_sequence == seq_type)
      bin_count_seq <- sum(seq_data$bin_count, na.rm = TRUE)
      total_bin_stable <- total_bin_stable + bin_count_seq
      
      stable_seq_list[[i]] <- data.frame(
        sequence = seq_type,
        type = "Stable",
        bin_count = bin_count_seq,
        relative_abundance = stable_sequences$total_sequence_RA[i]
      )
    }
  }
  changed_sequences <- sequence_totals %>%
    mutate(
      processes = strsplit(process_sequence, "->"),
      is_stable = sapply(processes, function(x) length(unique(x)) == 1)
    ) %>%
    filter(!is_stable)
  total_bin_changed <- 0
  changed_seq_list <- list()
  
  if (nrow(changed_sequences) > 0) {
    for (i in 1:nrow(changed_sequences)) {
      seq_type <- changed_sequences$process_sequence[i]
      seq_data <- stacked_data %>% filter(process_sequence == seq_type)
      bin_count_seq <- sum(seq_data$bin_count, na.rm = TRUE)
      total_bin_changed <- total_bin_changed + bin_count_seq
      
      changed_seq_list[[i]] <- data.frame(
        sequence = seq_type,
        type = "Changed",
        bin_count = bin_count_seq,
        relative_abundance = changed_sequences$total_sequence_RA[i]
      )
    }
  }
  
  all_sequences_info <- do.call(rbind, c(stable_seq_list, changed_seq_list))
  
  total_stable <- sum(all_sequences_info$bin_count[all_sequences_info$type == "Stable"], na.rm = TRUE)
  total_changed <- sum(all_sequences_info$bin_count[all_sequences_info$type == "Changed"], na.rm = TRUE)
  
  total_all <- total_stable + total_changed
  stable_ratio <- total_stable / total_all * 100
  changed_ratio <- total_changed / total_all * 100
  
  piechart_data <- data.frame(
    type = c("Stable", "Changed"),
    bin_count = c(total_stable, total_changed),
    percentage = c(stable_ratio, changed_ratio),
    label = c(
      paste0("Stable\n", total_stable, " bins\n", round(stable_ratio, 1), "%"),
      paste0("Changed\n", total_changed, " bins\n", round(changed_ratio, 1), "%")
    )
  )
  
  piechart_details <- list(
    summary = piechart_data,
    sequences = all_sequences_info,
    stats = data.frame(
      metric = c(
        "Total bins (stable)",
        "Total bins (changed)",
        "Total bins (all)",
        "Percentage (stable)",
        "Percentage (changed)"
      ),
      value = c(
        total_stable,
        total_changed,
        total_all,
        paste0(round(stable_ratio, 2), "%"),
        paste0(round(changed_ratio, 2), "%")
      )
    )
  )
  
} else {
  cat("Sequence type data not available\n\n")
  piechart_data <- NULL
  piechart_details <- NULL
}

if (!is.null(piechart_data)) {
  library(ggplot2)
  
  pie_chart <- ggplot(piechart_data, aes(x = "", y = bin_count, fill = type)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5),
              size = 5, color = "white", fontface = "bold") +
    scale_fill_manual(values = c("#4E79A7", "#F28E2B")) +
    labs(title = "Stable vs Changed bins",
         subtitle = paste("Based on", sum(piechart_data$bin_count), "total bins")
    ) +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          legend.position = "right")
  print(pie_chart)
  # Save pie chart
  ggsave("Rhizobia_Nitrogen_piechart.pdf", pie_chart, width = 8, height = 8, dpi = 300)
}
