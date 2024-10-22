#####################################
## Code to generate Supplementary Figure 1
##      in the manuscript
##      "Is inflammaging a universal aging process across human populations?"
##      by Franck et al.   October 2024
##
##
## Note: The authors do not have permission to share the data used in the manuscript
##       Simulated data is provided so that the code can be run but the results are
##       not expected to match those in the paper
#####################################


##########################
## Load libraries

library(FactoMineR) ## for PCA
library(corrplot)   ## for correlation plots
library(psych)      ## for corr.test function
library(ggplot2)    ## for biplots
library(ggrepel)    ## for biplots
library(cowplot)    ## for biplots
library(matrixcalc) ## for inverting a matrix
library(splines)    ## for cubic spline
library(broom)      ## for working with GLM
library(vioplot)
library(RColorBrewer)

#########################
## Source helper functions

source("graph_helperFxns.R")
source("Figures_2-3-4.R")    ##we will reuse some of this data

##########################
## Loaded data files in "Figures-2_3_4.R"

## We have two simulated datasets of 8 biomarkers
## dat1 will be the core dataset (analogous to InCHIANTI in the manuscript)
#dat1 <- readRDS("sampleData1.RDS")  
#dat2 <- readRDS("sampleData2.RDS")  
## 

################################
##
## Figure S1
## PCA biplots
## (Figure S3 uses the same code but different data)
## 
################################

## Factor analysis on dataset 1
##  first on all 8 biomarkers
pca1_all <- PCA(dat1[,1:8])
##  then on a subset of 6 biomarkers
pca1_subset <- PCA(dat1[,1:6])

## Factor analysis on dataset 2
pca2_all <- PCA(dat2[,1:8])
## Find maximum values for PC1 and PC2 contributions in dat1 using all biomarkers
max_PC1 <- max((pca1_all$var$coord[,"Dim.1"]), na.rm = TRUE)
max_PC2 <- max((pca1_all$var$coord[,"Dim.2"]), na.rm = TRUE)
min_PC1 <- min((pca1_all$var$coord[,"Dim.1"]), na.rm = TRUE)
min_PC2 <- min((pca1_all$var$coord[,"Dim.2"]), na.rm = TRUE)

# Dark to pale red for PC1
color_PC1 <- scales::rescale(pca1_all$var$coord[,"Dim.1"], to = c(0.2, 1), from = c(min_PC1, max_PC1))
# Dark to light blue for PC2
color_PC2 <- scales::rescale(pca1_all$var$coord[,"Dim.2"], to = c(0.2, 1), from = c(min_PC1, max_PC2))

## Mix the colors of PC1 (red) and PC2 (blue) to obtain the final color for each point
color_final <- mapply(mix_colors, color_PC1, color_PC2)
colorFinal <- data.frame(color=color_final, bm = rownames(pca1_all$var$coord))
cGray <- "#A6A6A6"

pca_biplot <- list()
for(i in 1:3){
  if(i == 1){
    fa_loadings <- pca1_all$var$coord[,1:2]
    colnames(fa_loadings) <- c("Factor1", "Factor2")
    factor1_variance <- round(pca1_all$eig[1,2],1)
    factor2_variance <- round(pca1_all$eig[2,2],1)
  } else if(i == 2){
    fa_loadings <- pca1_subset$var$coord[,1:2]
    colnames(fa_loadings) <- c("Factor1", "Factor2")
    factor1_variance <- round(pca1_subset$eig[1,2],1)
    factor2_variance <- round(pca1_subset$eig[2,2],1)
  } else if(i == 3){
    fa_loadings <- pca2_all$var$coord[,1:2]
    colnames(fa_loadings) <- c("Factor1", "Factor2")
    factor1_variance <- round(pca2_all$eig[1,2],1)
    factor2_variance <- round(pca2_all$eig[2,2],1)
  }
  loadings_df <- as.data.frame(fa_loadings)
  loadings_df$Variable <- rownames(loadings_df)
  loadings_df$color_final <- apply(as.array(loadings_df$Variable), MARGIN=1, FUN=colorLookup, colorFinal)
  
  ## GENERATE GRAPH WITH FINAL COLORS
  arrow_size <- 0.5
  pca_biplot[[i]] <- ggplot(loadings_df, aes(x = 0, y = 0, xend = Factor1, yend = Factor2, label = Variable)) +
    geom_segment(aes(color = I(color_final)), arrow = arrow(type = "closed", length = unit(0.1,"inches")), size = arrow_size) +  
    geom_text_repel(aes(color = I(color_final), x = Factor1, y = Factor2), size = 7, box.padding = unit(0.35, "lines")) +
    scale_color_identity() +
    labs(x = paste0("Factor 1 (", round(factor1_variance, 1), "%)"), 
         y = paste0("Factor 2 (", round(factor2_variance, 1), "%)")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_blank(),  
      axis.text.x = element_blank(), 
      axis.text.y = element_blank(), 
      axis.ticks = element_blank(),  
      axis.text = element_text(size = 12),
      axis.title.x = element_text(color = "black", size = 20, face = "bold"),
      axis.title.y = element_text(color = "black", size = 20, face = "bold"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      plot.background = element_rect(fill = "white", colour = NA), 
      plot.margin = margin(1, 1, 1, 1, "cm")
    ) +
    annotate("text", x = Inf, y = Inf, label = "", color = "#800020", 
             hjust = 1, vjust = 1, size = 6, fontface = "bold") 
}
## Combine the plots
length(pca_biplot)
plot_grid(pca_biplot[[1]], pca_biplot[[2]], pca_biplot[[3]],
          ncol=3,
          labels=c("Data1-All", "Data1-Subset", "Data2-All"), label_size=17, 
          label_colour="black", hjust=c(-2, -1, -2))



###############################
##
## Figure S2
## Factor score violin plot
##
###############################

c1 <- brewer.pal(n=8, name="Set2")
c2 <- brewer.pal(n=8, name="Set1")

par(mfrow=c(1,2))
vioplot(dat1$f1_score, dat2$f1_score, 
        main = "Inflamm-aging factor scores",  col=c2[1:2], 
        names=c("Dat1", "Dat2") , cex.axis=1.5,
        rectCol=c("indianred", "steelblue2"))

###############################
##
## Figure S4
## Biomarker violin plots
##
###############################


c1 <- brewer.pal(n=8, name="Set2")
c2 <- brewer.pal(n=8, name="Set1")

par(mfrow=c(2,4))
for(i in 1:8){##loop over biomarkers
  
  vioplot(dat1[,1], dat2[,1], 
        main = colnames(dat1)[i],  col=c2[1:2], 
        names=c("Dat1", "Dat2") , cex.axis=1.5,
        rectCol=c("indianred", "steelblue2"))
  
}

