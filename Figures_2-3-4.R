#####################################
## Code to generate Figures 2, 3 and 4
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

library(corrplot)   ## for correlation plots
library(psych)      ## for corr.test function
library(ggplot2)    ## for biplots
library(ggrepel)    ## for biplots
library(cowplot)    ## for biplots
library(matrixcalc) ## for inverting a matrix
library(splines)    ## for cubic spline
library(broom)      ## for working with GLM

#########################
## Source helper functions

source("graph_helperFxns.R")

##########################
## Load data files

## We have two simulated datasets of 8 biomarkers
## dat1 will be the core dataset (analogous to InCHIANTI in the manuscript)
dat1 <- readRDS("sampleData1.RDS")  
dat2 <- readRDS("sampleData2.RDS")  


###############################
##
## Figure 2
## Biomarker correlation plot
##
###############################

m1.dat <- cor(dat1, use="pairwise.complete.obs", method="spearman")
m1.dat.p <- corr.test(dat1, method="spearman", adjust="none") 

m2.dat <- cor(dat2, use="pairwise.complete.obs", method="spearman")
m2.dat.p <- corr.test(dat2, method="spearman", adjust="none") 

par(mfrow=c(1,2))
corrplot(m1.dat, tl.col="black", tl.cex=1.5, cl.cex=1.3, diag=FALSE,  cl.pos="n",
         mar=c(1.5, 0, 0, 0), sig.level=0.05, p.mat=m1.dat.p$p, insig="pch", pch=4, pch.col="gray", pch.cex=2.75)

corrplot(m2.dat, tl.col="black", tl.cex=1.5, cl.cex=1.3, diag=FALSE,  cl.pos="n",
         mar=c(1.5, 0, 0, 0), sig.level=0.05, p.mat=m2.dat.p$p, insig="pch", pch=4, pch.col="gray", pch.cex=2.75)


################################
##
## Figure 3 (A- J)
## Factor analysis biplots
## 
################################

## Factor analysis on dataset 1
##  first on all 8 biomarkers
fa1_all <- factanal(dat1[,1:8], factors=2, rotation="varimax", scores="regression")
##  then on a subset of 6 biomarkers
fa1_subset <- factanal(dat1[,1:6], factors=2, rotation="varimax", scores="regression")

## Factor analysis on dataset 2
fa2_all <- factanal(dat2[,1:8], factors=2, rotation="varimax", scores="regression")

## Find maximum values for PC1 and PC2 contributions in dat1 using all biomarkers
max_PC1 <- max((fa1_all$loadings[,"Factor1"]), na.rm = TRUE)
max_PC2 <- max((fa1_all$loadings[,"Factor2"]), na.rm = TRUE)
min_PC1 <- min((fa1_all$loadings[,"Factor1"]), na.rm = TRUE)
min_PC2 <- min((fa1_all$loadings[,"Factor2"]), na.rm = TRUE)

# Dark to pale red for PC1
color_PC1 <- scales::rescale(fa1_all$loadings[,"Factor1"], to = c(0.2, 1), from = c(min_PC1, max_PC1))
# Dark to light blue for PC2
color_PC2 <- scales::rescale(fa1_all$loadings[,"Factor2"], to = c(0.2, 1), from = c(min_PC1, max_PC2))

## Mix the colors of PC1 (red) and PC2 (blue) to obtain the final color for each point
color_final <- mapply(mix_colors, color_PC1, color_PC2)
colorFinal <- data.frame(color=color_final, bm = rownames(fa1_all$loadings))
cGray <- "#A6A6A6"

pca_biplot <- list()
for(i in 1:3){
  if(i == 1){
    fa_loadings <- fa1_all$loadings[,1:2]
    total_variance <- sum(fa_loadings^2) + sum(fa1_all$uniquenesses)
    factor1_variance <- sum(fa_loadings[, "Factor1"]^2) / total_variance * 100
    factor2_variance <- sum(fa_loadings[, "Factor2"]^2) / total_variance * 100
  } else if(i == 2){
    fa_loadings <- fa1_subset$loadings[,1:2]
    total_variance <- sum(fa_loadings^2) + sum(fa1_subset$uniquenesses)
    factor1_variance <- sum(fa_loadings[, "Factor1"]^2) / total_variance * 100
    factor2_variance <- sum(fa_loadings[, "Factor2"]^2) / total_variance * 100
  } else if(i == 3){
    fa_loadings <- fa2_all$loadings[,1:2]
    total_variance <- sum(fa_loadings^2) + sum(fa2_all$uniquenesses)
    factor1_variance <- sum(fa_loadings[, "Factor1"]^2) / total_variance * 100
    factor2_variance <- sum(fa_loadings[, "Factor2"]^2) / total_variance * 100
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


################################
##
## Figure 3 (K)
## Axis stability plots
## 
################################

## The idea is to perform the factor analysis on subgroups: female, male, <65, > 65

## First for dat1
## Female
F1 <- subset(dat1, Female==1)
dat1.F <- factanal(F1[,1:8],2, rotation="varimax", scores = "regression") 
## Male
M1 <- subset(dat1, Female==0)
dat1.M <- factanal(M1[,1:8],2, rotation="varimax", scores = "regression") 
## Less than 65 years of age
L1 <- subset(dat1, Age < 65)
dat1.L <- factanal(L1[,1:8],2, rotation="varimax", scores = "regression") 
## 65 years of age or more
S1 <- subset(dat1, Age >= 65)
dat1.S <- factanal(S1[,1:8],2, rotation="varimax", scores = "regression") 

## Then for dat2
## Female
F2 <- subset(dat2, Female==1)
dat2.F <- factanal(F2[,1:8],2, rotation="varimax", scores = "regression") 
## Male
M2 <- subset(dat2, Female==0)
dat2.M <- factanal(M2[,1:8],2, rotation="varimax", scores = "regression") 
## Less than 65 years of age
L2 <- subset(dat2, Age < 65)
dat2.L <- factanal(L2[,1:8],2, rotation="varimax", scores = "regression") 
## 65 years of age or more
S2 <- subset(dat2, Age >= 65)
dat2.S <- factanal(S2[,1:8],2, rotation="varimax", scores = "regression") 

par(mfrow=c(1,2))
## Create correlation matrix for Factor 1 loadings
## First for dat1
stab1 <- data.frame(biomarkers = rownames(fa1_all$loadings), "1-All" = fa1_all$loadings[,"Factor1"])
stab1$"1-Female" <-dat1.F$loadings[,"Factor1"]
stab1$"1-Male" <- dat1.M$loadings[,"Factor1"]
stab1$"1-Age<65" <- dat1.L$loadings[,"Factor1"]
stab1$"1-Age 65+" <- dat1.S$loadings[,"Factor1"]
  
stab1.corr <- cor(stab1[2:6])
corrplot(stab1.corr, tl.col="black", tl.cex=2, cl.cex=1.3,
         sig.level=0.05, insig="pch", pch=4, pch.col="gray", pch.cex=1)


## First for dat2
stab2 <- data.frame(biomarkers = rownames(fa2_all$loadings), "2-All" = fa2_all$loadings[,"Factor1"])
stab2$"2-Female" <-dat2.F$loadings[,"Factor1"]
stab2$"2-Male" <- dat2.M$loadings[,"Factor1"]
stab2$"2-Age<65" <- dat2.L$loadings[,"Factor1"]
stab2$"2-Age 65+" <- dat2.S$loadings[,"Factor1"]

stab2.corr <- cor(stab2[2:6])
corrplot(stab2.corr, tl.col="black", tl.cex=2, cl.cex=1.3,
         sig.level=0.05, insig="pch", pch=4, pch.col="gray", pch.cex=1)


################################
##
## Figure 4 (A)
## Age versus score
## 
################################

## Get age and factor 1 score for each sample for each dataset

dat1.scores <- data.frame(Age=dat1$Age, f1_score=fa1_all$scores[,"Factor1"], f2_score=fa1_all$scores[,"Factor2"])
dat2.scores <- data.frame(Age=dat2$Age, f1_score=fa2_all$scores[,"Factor1"], f2_score=fa2_all$scores[,"Factor2"])

## Get best fit lines we will need for the plot
lm1 <- lm(dat1.scores$f1_score~dat1.scores$Age)
out1 <- as.data.frame(predict(lm1, newdata=dat1.scores, interval="confidence"))
out1$Age <- dat1.scores$Age
out1 <- out1[order(out1$Age),]

lm2 <- lm(dat2.scores$f1_score~dat2.scores$Age)
out2 <- as.data.frame(predict(lm2, newdata=dat2.scores, interval="confidence"))
out2$Age <- dat2.scores$Age
out2 <- out2[order(out2$Age),]

##Get some colors
library(RColorBrewer)
c1 <- brewer.pal(n=8, name="Set2")
c2 <- brewer.pal(n=8, name="Set1")
par(mfrow=c(1,2))
## Dat1
plot(x=out1$Age, y=out1$fit, col=c2[1], 
     type="l", ylim=c(-1, 0.5), xlim=c(20, 100), lwd=2.8, cex.lab=1.9, cex.axis=1.6,
     ylab="", xlab="")
title(ylab="Inflamm-aging factor score", line=2.5, cex.lab=1.9, font=2)
title(xlab="Age (years)", line=2.5, cex.lab=1.9, font=2)
polygon(x = c(x=out1$Age, rev(x=out1$Age)),
        y = c(out1$upr, rev(out1$lwr)),
        col =  adjustcolor(c2[1], alpha.f = 0.10), border = NA)

## Dat2
lines(x=out2$Age, y=out2$fit, type="l", col=c2[2], lwd=2.5)
polygon(x=c(out2$Age, rev(out2$Age)),
        y=c(out2$upr, rev(out2$lwr)),
        col=adjustcolor(c2[2], alpha.f=0.08), border=NA)
# Legend
t1 <- paste0("Dat1: slope=0.004")  ##Information from summary(lm1)
t2 <- paste0("Dat2: slope=-0.001") ##Information from summary(lm2)
legend("bottomright", legend=c(t1, t2), 
       col=c2[1:2], lwd=2, cex=1.5)


################################
##
## Figure 4 (B-D)
## Association with health outcomes
## 
################################

## Assemble the data

## We use dat1 factor1 scores as is because dat1 is our main dataset
## We recalculate dat2 scores using dat1 loadings
## (This is analogous to applying the InCHIANTI loadings to another
##  dataset in the manuscript)

rownames(fa1_all$loadings) # dat1
colnames(dat2[,1:8]) # They match
dat2_matrix <- as.matrix(dat2[,1:8])
load1_matrix <- as.matrix(fa1_all$loadings[,1:2])
sigma_inv <- matrix.inverse(fa1_all$correlation)
W <- t(load1_matrix) %*% sigma_inv
projected_data <- W %*% t(dat2_matrix) 
projected_data <- as.data.frame(t(projected_data))

## Now, simulate some health outcomes 
## Let H1 will be associated with Factor 1 in both datasets
dat1[,c("f1_score", "f2_score")] <- scale(dat1.scores[,c("f1_score", "f2_score")])
dat1$H1 <- rbinom(n=nrow(dat1), size=1, prob=ifelse(dat1$f1_score>0, 0.75, 0.25))
dat2[,c("f1_score_proj", "f2_score_proj")] <- projected_data[,c("Factor1", "Factor2")]
dat2$H1 <- rbinom(n=nrow(dat2), size=1, prob=ifelse(dat2$f1_score_proj>0, 0.75, 0.25))

## Let H2 be associated with Factor 1 in dat1 only
dat1$H2 <- rbinom(n=nrow(dat1), size=1, prob=ifelse(dat1$f1_score>0, 0.85, 0.25))
dat2$H2 <- rbinom(n=nrow(dat2), size=1, prob=0.5)

## Let H3 be mildly associated with Factor 1 in dat1,2
dat1$H3 <- rbinom(n=nrow(dat1), size=1, prob=ifelse(dat1$f1_score>0, 0.41, 0.59))
dat2$H3 <- rbinom(n=nrow(dat2), size=1, prob=ifelse(dat2$f1_score_proj>0, 0.45, 0.55))

## Now run logistic regression models to uncover the association

# Prepare table for results
results <- data.frame(Outcome = character(), Controls = character(), Measure = character(),
                      Effect = numeric(), p_value = numeric(), LCI = numeric(), UCI = numeric(),
                      Sign_level = integer(), stringsAsFactors = FALSE)


for (comorbidity in c("H1", "H2", "H3")) {
  
  model1 <- glm(dat1[,comorbidity] ~ f1_score + f2_score + bs(Age) + Female, family = binomial, data = dat1)
  model2 <- glm(dat2[,comorbidity] ~ f1_score_proj + f2_score_proj + bs(Age) + Female, family = binomial, data = dat2)
  tidy1 <- tidy(model1, conf.int = TRUE)
  tidy2 <- tidy(model2, conf.int = TRUE)
  for (i in 2:3) {  ## model1
    p_value_1 <- tidy1$p.value[i]
    results <- rbind(results, data.frame(Outcome = comorbidity, 
                                         Controls = "Model 1",
                                         Measure = names(model1$coefficients)[i],
                                         Effect = exp(tidy1$estimate[i]),
                                         p_value = p_value_1,
                                         LCI = exp(tidy1$conf.low[i]),
                                         UCI = exp(tidy1$conf.high[i]),
                                         Sign_level = sign_level(p_value_1)))
  }
  for (i in 2:3) { ## model2
    p_value_2 <- tidy2$p.value[i]
    results <- rbind(results, data.frame(Outcome = comorbidity, 
                                         Controls = "Model 2",
                                         Measure = names(model2$coefficients)[i],
                                         Effect = exp(tidy2$estimate[i]),
                                         p_value = p_value_2,
                                         LCI = exp(tidy2$conf.low[i]),
                                         UCI = exp(tidy2$conf.high[i]),
                                         Sign_level = sign_level(p_value_2)))
  }
}

## Now create the forest plot

# Definition of colors for different diseases
outcome_colors <- c("H1" = "#377EB8",
                    "H2" = "#984EA3",
                    "H3" = "#FF7F00")

# Line types for different cohorts
line_types <- c("f1_score" = 2,
                "f1_score_proj" = 1)

point_types <- c("f1_score" = 16,
                 "f1_score_proj" = 17)

# Filter data to include selected measures and outcomes
selected_measures <- c("f1_score","f1_score_proj")
ordered_outcomes <- c("H3", "H2", "H1")
# Filter by measurements and order of diseases
data_filtered1 <- results[results$Measure %in% selected_measures & results$Outcome %in% ordered_outcomes,]


common_xlim <- c(min(data_filtered1$LCI, na.rm = TRUE), 6.5)

plot_frame(data_filtered1, "", add_legend = TRUE, x_axis_title = "Odds Ratio", 
           lw=2, figLetter="", c_names=c("Data2", "Data1"))




