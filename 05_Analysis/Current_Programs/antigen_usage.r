library(ggplot2)

score_df <- read.table("/Users/austinseamann/PycharmProjects/TCRDOCK/AntigenList/results.tsv", header=T)
df_tcr <- c()
df_antigen <- c()

for (i in 1:nrow(score_df)){
  temp_var <- strsplit(as.character(score_df[i,1]), "_")
  df_tcr <- c(df_tcr, temp_var[[1]][1])
  df_antigen <- c(df_antigen, temp_var[[1]][2])
}
score_df$tcr <- df_tcr
score_df$antigen <- df_antigen

ggplot(data = score_df, aes(x=tcr, y=I_sc)) + geom_boxplot(aes(fill=tcr)) + labs(title = "Antigen Usage")


######DENSITY PLOT FUNCTION OR SOMETHING##########3
plot(density(score_df[score_df$tcr==unique(score_df$tcr)[1],]$I_sc))

for (i in unique(score_df$tcr)[2:length(unique(score_df$tcr))]){
  lines(density(score_df[score_df$tcr==i,]$I_sc))
}
legend("topright", unique(score_df$tcr), fill=1:length(unique(score_df$tcr)))


TCR_df <- score_df[score_df$TCR == "4prh",]

# I have restricted this value to 0, so it can never be lower :)
minAlphaBeta <- 0
maxAlphaBeta <- max(TCR_df$alpha, TCR_df$beta)
difAlphaBeta <- maxAlphaBeta - minAlphaBeta

minIntScore <- min(TCR_df$I_sc)
maxIntScore <- max(TCR_df$I_sc)
difIntScore <- abs(maxIntScore - minIntScore)

colors <- c("Interface Score" = "red", "Alpha Score" = "orange", "Beta Score" = "purple")
shapes <- c("Interface Score" = "triangle", "Alpha Score" = "circle", "Beta Score" = "circle")

ggplot(TCR_df, aes(x=pMHC, y = I_sc, group=I_sc, colour=I_sc)) + 
  
  geom_point(aes(y=I_sc, color = "Interface Score", shape = "Interface Score"), size = 5) + 
  geom_point(aes(y=alpha * difIntScore/difAlphaBeta + minIntScore, color = "Alpha Score", shape = "Alpha Score"), size = 3) +
  geom_point(aes(y=beta * difIntScore/difAlphaBeta + minIntScore, color = "Beta Score", shape = "Beta Score"), size = 3) +
  
  scale_y_continuous(name = "Interface Score", 
                     sec.axis = sec_axis(trans = ~.*difAlphaBeta/difIntScore + 
                                           -1*difAlphaBeta/difIntScore*minIntScore, 
                                         name = "Alpha/Beta Scores")) + 
  labs(color = "Legend") + scale_color_manual(name = "Legend", values = colors) + scale_shape_manual(name = "Legend", values = shapes)

