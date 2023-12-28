#Analysis of nCounter Data
#Prior to analysis please pass RCC data through Rosalind to perform normalization and differential analysis
library(EnhancedVolcano)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(openxlsx)
library(org.Hs.eg.db)
library(ggrepel)

library("png")
rm(list=ls())

files<- c("statistic_1vs2.xlsx", "statistic_1vs3.xlsx", "statistic_1vs4.xlsx",
          "statistic_2vs3.xlsx", "statistic_2vs4.xlsx", "statistic_3vs4.xlsx")
titles <- c("First vs Second", "First vs Third", "First vs Fourth",
            "Second vs Third", "Second vs Fourth",
            "Third vs Fourth")
namespng <- c("Differettil expressed genes First vs Second.png", 
              "Differettil expressed genes First vs Third.png", 
              "Differettil expressed genes First vs Fourth.png",
              "Differettil expressed genes Second vs Third.png", 
              "Differettil expressed genes Second vs Fourth.png",
              "Differettil expressed genes Third vs Fourth.png")

i=6
j=-0.48
k=.48
for (i in 1:length(files)) {
  statistic <- read.xlsx(files[i])
  statistic$DEG <- "NO"
  statistic$DEG[statistic$log2fc < j & statistic$p <.05] <- "DOWN"
  statistic$DEG[statistic$log2fc > k & statistic$p <.06] <- "UP"

  statistic$delabel <- NA
  statistic$delabel[statistic$DEG != "NO"] <- statistic$.y.[statistic$DEG != "NO"]
  
  img <- ggplot(data=statistic, aes(x=log2fc, y=-log10(p), col= DEG, label=delabel)) +
    geom_point() + 
    theme_minimal() +
    geom_text_repel() +
    geom_vline(xintercept=c(-0.48, 0.48), col="black", linetype=2) +
    geom_hline(yintercept=-log10(0.005), col="black", linetype=2)+
    xlab("Log2 Fold Change")+
    ggtitle(titles[i])
  
  namepng<- namespng[i]
  png(filename=namepng)
  plot(img)
  dev.off()
  i=i+1
}
#################### naming the important genes only ##########
for (i in 1:length(files)) {
  statistic <- read.xlsx(files[i])
  statistic$DEG <- "NO"
  statistic$DEG[statistic$log2fc < j & statistic$p <.05] <- "DOWN"
  statistic$DEG[statistic$log2fc > k & statistic$p <.06] <- "UP"
  
  statistic$delabel <- NA
  statistic$delabel[statistic$DEG != "NO"] <- statistic$.y.[statistic$DEG != "NO"]
  statistic$genelabels <- ifelse(statistic$delabel == "CDK7"
                                    | statistic$delabel == "MGMT", statistic$delabel,  "")
  
  img <- ggplot(data=statistic, aes(x=log2fc, y=-log10(p), col= DEG, label=genelabels)) +
    geom_point() + 
    theme_minimal() +
    geom_text_repel() +
    geom_vline(xintercept=c(-0.48, 0.48), col="black", linetype=2) +
    geom_hline(yintercept=-log10(0.05), col="black", linetype=2)+
    xlab("Log2 Fold Change")+
    ggtitle(titles[i])
  
  namepng<- namespng[i]
  png(filename=namepng)
  plot(img)
  dev.off()
  i=i+1
}

