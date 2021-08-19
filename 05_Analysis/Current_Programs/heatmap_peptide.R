library(pheatmap)
library(RColorBrewer)
library(viridisLite)

# Heatmap from csv

df <- read.csv("~/PycharmProjects/TCRDOCK/Fran_Breakdown/Docking/1OGA/Heatmap/Peptide/heatmap.csv", row.names = 1, stringsAsFactors = FALSE)
df[df >= 0] <- 0.0


# Adjusting breaks
bk = unique(c(-0.00001, -0.0001,-0.01, seq(-0.01, min(df) - 1, length=7)))


# Generate the alpha and beta labeling based on alpha until backward count
tcr_aa <- colnames(df)  # All labels of TCR AA
numbers <- tcr_aa  # Numbers correlating to TCR AA

# Removes everything but AA num
for (i in 1:length(tcr_aa)){
  split <- strsplit(tcr_aa[i], "\\.")
  numbers[i] <- split[[1]][2]
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

colCDRs <- as.data.frame(CDR)
row.names(colCDRs) <- colnames(df)


# Set CDR Labels
ann_color <- list("CDR" = c("alpha cdr1" = "pink", "alpha cdr2" = "red", "alpha cdr3" = "darkred", "beta cdr1" = "cyan", "beta cdr2" = "blue", "beta cdr3" = "darkblue"))

#pheatmap(df, main = "Energy Breakdown", color = c(inferno(length(bk - 1)), "white"), breaks = rev(bk), cluster_cols = F, cluster_rows = F, annotation_col = colLabels)
# Orignal
#pheatmap(matrixDf, main = "Energy Breakdown", color = colorRampPalette(c(rev(brewer.pal(n=length(bk - 1), name="YlOrRd")), "white"))(length(bk)), breaks = rev(bk), cluster_cols = F, cluster_rows = F, annotation_col = colCDRs)
# Updated
pheatmap(t(df), main = "Peptide", color = colorRampPalette(c("darkgreen", "white"))(length(bk)), breaks = rev(bk), cluster_cols = F, cluster_rows = F, annotation_row = colCDRs, annotation_colors = ann_color)

#color = colorRampPalette(rev(c("navy","white","red")))(50)
#color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(100)
#pheatmap(df, main = "Energy_Breakdown", breaks = bk, color = rev(mycol), cluster_cols = F, cluster_rows = F, annotation_col = colLabels)

