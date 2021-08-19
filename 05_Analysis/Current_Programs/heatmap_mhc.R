library(pheatmap)
library(RColorBrewer)
library(viridisLite)

# Heatmap from csv

df <- read.csv("~/PycharmProjects/TCRDOCK/Fran_Breakdown/Docking/1OGA/Heatmap/MHC/Heatmap.csv", row.names = 1, stringsAsFactors = FALSE)
df[df >= 0] <- 0.0  # Removes positive values
df <- df[rowSums(df[])<0,]  # Removes rows with only zeros for interactions between MHC and TCR chains


#bk = unique(c(seq(-1,min(df) - 1, length=5),0.1,-0.0001,rev(seq(1,max(df) + 1, length=5))))
bk = unique(c(-0.00001, -0.0001,-0.01, seq(-0.01, min(df) - 1, length=7)))


# Generate the alpha and beta labeling based on alpha until backward count
tcr_aa <- colnames(df)  # All labels of TCR AA
numbers <- tcr_aa  # Numbers correlating to TCR AA

# Removes everything but AA num
for (i in 1:length(tcr_aa)){
  split <- strsplit(tcr_aa[i], "\\.")
  numbers[i] <- split[[1]][2]
}

# Pull numbers of AA from MHC AA
mhc_aa <- rownames(df)
mhc_num <- mhc_aa
for (i in 1:length(mhc_aa)) {
  split <- strsplit(mhc_aa[i], " ")
  mhc_num[i] <- split[[1]][2]
}



# Lables alpha beta based on alpha first until backward count
CDR <- rep("alpha cdr1", length(numbers))
count <- numbers[1]
cdrCount <- 1
tempLabel <- "alpha cdr"
flag <- TRUE
for (i in 2:length(numbers)){
  if (as.numeric(numbers[i]) < as.numeric(numbers[i - 1])){
    tempLabel <- "beta cdr"
  }
  CDR[i] <- tempLabel
}
for (i in 2:(length(numbers))){
  if ((CDR[i] == "beta cdr") & (flag)){
    cdrCount = 1
    flag <- FALSE
  }
  if ((as.numeric(numbers[i]) - as.numeric(numbers[i - 1])) > 1) {
    cdrCount <- cdrCount + 1
  }
  CDR[i] <- paste0(CDR[i], cdrCount)
}

# Adjust for two heatmaps, one AH1, another AH2
flag1 <- FALSE
cut <- 0
for (i in 2:length(mhc_num)){
  if ((as.numeric(mhc_num[i]) - as.numeric(mhc_num[i - 1])) > 40) {
    cut <- as.numeric(i) - 1
  }
}

trDF <- t(df)
df_ah1 <- trDF[,c(1:cut)]
df_ah2 <- trDF[,c((cut + 1):ncol(trDF))]

colabLabels <- as.data.frame(CDR)
row.names(colabLabels) <- colnames(df)

# Setting colors for CDR regions
ann_color <- list("CDR" = c("alpha cdr1" = "pink", "alpha cdr2" = "red", "alpha cdr3" = "darkred", "beta cdr1" = "cyan", "beta cdr2" = "blue", "beta cdr3" = "darkblue"))

pheatmap(t(df), main = "MHC Alpha Helix 1", color = colorRampPalette(c("darkgreen", "white"))(length(bk)), breaks = rev(bk), cluster_cols = F, cluster_rows = F, annotation_row = colabLabels, annotation_colors = ann_color)

#AH1
pheatmap(df_ah1, main = "MHC Alpha Helix 1", color = colorRampPalette(c("darkgreen", "white"))(length(bk)), breaks = rev(bk), cluster_cols = F, cluster_rows = F, annotation_row = colabLabels, annotation_colors = ann_color)
#AH2
pheatmap(df_ah2, main = "MHC Alpha Helix 2", color = colorRampPalette(c("darkgreen", "white"))(length(bk)), breaks = rev(bk), cluster_cols = F, cluster_rows = F, annotation_row = colabLabels, annotation_colors = ann_color)

