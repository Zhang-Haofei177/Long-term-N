rm(list = ls())
setwd("C:/Users/ADMIN/Desktop/Longterm N code/Fig.S5/Rhizobia/")

####### Ori compartment ##########
niche <- "Ori"
metadata <- read.table(paste0("mappingRhi_", niche, ".txt"), sep = '\t', h = TRUE, row.names = 1, check.names = FALSE, comment.char = '')

otus <- read.table('Rhizobia_ASVtable.txt', sep = '\t', header = TRUE, row.names = 1, check.names = FALSE, skip = 1, comment.char = '')
otus <- round(otus)
otus <- t(as.matrix(otus))
otus[is.na(otus)] <- 0

common.sample.ids <- intersect(rownames(metadata), rownames(otus))
otus <- otus[common.sample.ids, ]
metadata <- metadata[common.sample.ids, ]

# Error message if no common samples between mapping file and input file
if (length(common.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}

train.ix <- which(metadata$SourceSink == 'source')
test.ix <- which(metadata$SourceSink == 'sink')
envs <- metadata$Env

if (is.element('Description', colnames(metadata))) 
  desc <- metadata$Description

source('SourceTracker2.r')
alpha1 <- alpha2 <- 0.001

st <- sourcetracker(otus[train.ix, ], envs[train.ix])
st

results <- predict(st, otus[test.ix, ], alpha1 = alpha1, alpha2 = alpha2)
write.table(results$proportions, paste0("sourceTracker_", niche, ".txt"), sep = "\t", quote = FALSE)

p <- plot(results, type = 'pie', include.legend = TRUE, env.colors = c('#47697E', '#5B7444', '#CC6666', '#79BEDB', '#885588'))
p

####### RZ compartment ##########
niche <- "RZ"
metadata <- read.table(paste0("mappingRhi_", niche, ".txt"), sep = '\t', h = TRUE, row.names = 1, check.names = FALSE, comment.char = '')

otus <- read.table('Rhizobia_ASVtable.txt', sep = '\t', header = TRUE, row.names = 1, check.names = FALSE, skip = 1, comment.char = '')
otus <- round(otus)
otus <- t(as.matrix(otus))
otus[is.na(otus)] <- 0

common.sample.ids <- intersect(rownames(metadata), rownames(otus))
otus <- otus[common.sample.ids, ]
metadata <- metadata[common.sample.ids, ]

# Error message if no common samples between mapping file and input file
if (length(common.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}

train.ix <- which(metadata$SourceSink == 'source')
test.ix <- which(metadata$SourceSink == 'sink')
envs <- metadata$Env

if (is.element('Description', colnames(metadata))) 
  desc <- metadata$Description

source('SourceTracker2.r')
alpha1 <- alpha2 <- 0.001

st <- sourcetracker(otus[train.ix, ], envs[train.ix])
st

results <- predict(st, otus[test.ix, ], alpha1 = alpha1, alpha2 = alpha2)
write.table(results$proportions, paste0("sourceTracker_", niche, ".txt"), sep = "\t", quote = FALSE)

p <- plot(results, type = 'pie', include.legend = TRUE, env.colors = c('#47697E', '#5B7444', '#CC6666', '#79BEDB', '#885588'))
p

####### RH compartment ##########
niche <- "RH"
metadata <- read.table(paste0("mappingRhi_", niche, ".txt"), sep = '\t', h = TRUE, row.names = 1, check.names = FALSE, comment.char = '')

otus <- read.table('Rhizobia_ASVtable.txt', sep = '\t', header = TRUE, row.names = 1, check.names = FALSE, skip = 1, comment.char = '')
otus <- round(otus)
otus <- t(as.matrix(otus))
otus[is.na(otus)] <- 0

common.sample.ids <- intersect(rownames(metadata), rownames(otus))
otus <- otus[common.sample.ids, ]
metadata <- metadata[common.sample.ids, ]

# Error message if no common samples between mapping file and input file
if (length(common.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}

train.ix <- which(metadata$SourceSink == 'source')
test.ix <- which(metadata$SourceSink == 'sink')
envs <- metadata$Env

if (is.element('Description', colnames(metadata))) 
  desc <- metadata$Description

source('SourceTracker2.r')
alpha1 <- alpha2 <- 0.001

st <- sourcetracker(otus[train.ix, ], envs[train.ix])
st

results <- predict(st, otus[test.ix, ], alpha1 = alpha1, alpha2 = alpha2)
write.table(results$proportions, paste0("sourceTracker_", niche, ".txt"), sep = "\t", quote = FALSE)

p <- plot(results, type = 'pie', include.legend = TRUE, env.colors = c('#47697E', '#5B7444', '#CC6666', '#79BEDB', '#885588'))
p

####### RT compartment ##########
niche <- "RT"
metadata <- read.table(paste0("mappingRhi", niche, ".txt"), sep = '\t', h = TRUE, row.names = 1, check.names = FALSE, comment.char = '')

otus <- read.table('Rhizobia_ASVtable.txt', sep = '\t', header = TRUE, row.names = 1, check.names = FALSE, skip = 1, comment.char = '')
otus <- round(otus)
otus <- t(as.matrix(otus))
otus[is.na(otus)] <- 0

common.sample.ids <- intersect(rownames(metadata), rownames(otus))
otus <- otus[common.sample.ids, ]
metadata <- metadata[common.sample.ids, ]

# Error message if no common samples between mapping file and input file
if (length(common.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}

train.ix <- which(metadata$SourceSink == 'source')
test.ix <- which(metadata$SourceSink == 'sink')
envs <- metadata$Env

if (is.element('Description', colnames(metadata))) 
  desc <- metadata$Description

source('SourceTracker2.r')
alpha1 <- alpha2 <- 0.001

st <- sourcetracker(otus[train.ix, ], envs[train.ix])
st

results <- predict(st, otus[test.ix, ], alpha1 = alpha1, alpha2 = alpha2)
write.table(results$proportions, paste0("sourceTracker_", niche, ".txt"), sep = "\t", quote = FALSE)

p <- plot(results, type = 'pie', include.legend = TRUE, env.colors = c('#47697E', '#5B7444', '#CC6666', '#79BEDB', '#885588'))
p

####### Nod compartment ##########
niche <- "Nod"
metadata <- read.table(paste0("mappingRhi", niche, ".txt"), sep = '\t', h = TRUE, row.names = 1, check.names = FALSE, comment.char = '')

otus <- read.table('Rhizobia_ASVtable.txt', sep = '\t', header = TRUE, row.names = 1, check.names = FALSE, skip = 1, comment.char = '')
otus <- round(otus)
otus <- t(as.matrix(otus))
otus[is.na(otus)] <- 0

common.sample.ids <- intersect(rownames(metadata), rownames(otus))
otus <- otus[common.sample.ids, ]
metadata <- metadata[common.sample.ids, ]

# Error message if no common samples between mapping file and input file
if (length(common.sample.ids) <= 1) {
  message <- paste(sprintf('Error: there are %d sample ids in common '),
                   'between the metadata file and data table')
  stop(message)
}

train.ix <- which(metadata$SourceSink == 'source')
test.ix <- which(metadata$SourceSink == 'sink')
envs <- metadata$Env

if (is.element('Description', colnames(metadata))) 
  desc <- metadata$Description

source('SourceTracker2.r')
alpha1 <- alpha2 <- 0.001

st <- sourcetracker(otus[train.ix, ], envs[train.ix])
st

results <- predict(st, otus[test.ix, ], alpha1 = alpha1, alpha2 = alpha2)
write.table(results$proportions, paste0("sourceTracker_", niche, ".txt"), sep = "\t", quote = FALSE)

p <- plot(results, type = 'pie', include.legend = TRUE, env.colors = c('#47697E', '#5B7444', '#CC6666', '#79BEDB', '#885588'))
p

# Process and combine results for Rhizobia analysis
rz <- read.table("SourceTracker_RZ.txt", header = TRUE, sep = "\t")
colnames(rz) <- c("sink", colnames(rz)[1:5])
means <- colMeans(rz[, 2:6], na.rm = TRUE)
rz_result <- data.frame(
  sink = "rz",
  as.list(means))

rh <- read.table("SourceTracker_RH.txt", header = TRUE, sep = "\t")
colnames(rh) <- c("sink", colnames(rh)[1:5])
means <- colMeans(rh[, 2:6], na.rm = TRUE)
rh_result <- data.frame(
  sink = "rh",
  as.list(means))

rt <- read.table("SourceTracker_RT.txt", header = TRUE, sep = "\t")
colnames(rt) <- c("sink", colnames(rt)[1:5])
means <- colMeans(rt[, 2:6], na.rm = TRUE)
rt_result <- data.frame(
  sink = "rt",
  as.list(means))

nod <- read.table("SourceTracker_Nod.txt", header = TRUE, sep = "\t")
colnames(nod) <- c("sink", colnames(nod)[1:5])
means <- colMeans(nod[, 2:6], na.rm = TRUE)
nod_result <- data.frame(
  sink = "nod",
  as.list(means))

library(dplyr)
combined_result <- bind_rows(rz_result, rh_result, rt_result, nod_result)
combined_result[is.na(combined_result)] <- 0
print(combined_result)
write.csv(combined_result, "Rhizobia source track result.csv")