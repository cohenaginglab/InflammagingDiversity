################################
## Create a set of simulated data to be used
## with sample code for the manuscript
##      "Is inflammaging a universal aging process across human populations?"
##      by Franck et al.   October 2024
##
##
## We will create two sets of 8 simulated biomarkers
##    these biomarkers all have a Normal(0,1) distribution as if
##    they had been transformed and scaled in preparation for
##    principal components analysis
##    50% female
##    Age uniformly distributed 
##
#################################


############
## Load libraries

library(MASS)  ##for the mvrnorm function


############
## Define simulation parameters

N_samples <- 1000
mu <- rep(0, 8)
sigma1 <- matrix( 
  #  1     2   3    4      5    6    7    8   
  c(1.0, 0.4, 0.5, 0.6,   0.5, 0.4, 0.3, 0.2, 
    0.4, 1.0, 0.3, 0.2,   0.2, 0.1, 0.1, 0.3, 
    0.5, 0.3, 1.0, 0.4,   0.5, 0.6, 0.5, 0.2, 
    0.6, 0.2, 0.4, 1.0,   0.5, 0.0, 0.0, 0.0, 
    0.5, 0.2, 0.5, 0.5,   1.0, 0.3, 0.4, 0.1, 
    0.4, 0.1, 0.6, 0.0,   0.3, 1.0, 0.7, 0.4, 
    0.3, 0.1, 0.5, 0.0,   0.4, 0.7, 1.0, 0.4, 
    0.2, 0.3, 0.2, 0.0,   0.1, 0.3, 0.4, 1.0),
    nrow=8, ncol=8) 
    
set.seed(19)
df1<-as.data.frame(mvrnorm(n=N_samples, mu=mu, Sigma=sigma1))
df1$Female <- rbinom(n=N_samples, size=1, prob=0.5)
df1$Age <- runif(n=N_samples, min=50, max=85)

## check the data
head(df1)
cor(df1)

saveRDS(df1, file="sampleData1.RDS")


N_samples <- 1000
mu <- rep(0, 8)
sigma1 <- matrix( 
  #  1     2   3    4      5    6    7    8   
  c(1.0, 0.4, 0.3, 0.2,   0.5, 0.4, 0.3, 0.2, 
    0.4, 1.0, 0.3, 0.2,   0.2, 0.1, 0.1, 0.3, 
    0.3, 0.3, 1.0, 0.4,   0.2, 0.4, 0.5, 0.2, 
    0.2, 0.2, 0.4, 1.0,   0.5, 0.0, 0.0, 0.0, 
    0.5, 0.2, 0.2, 0.5,   1.0, 0.3, 0.4, -0.1, 
    0.4, 0.1, 0.4, 0.0,   0.3, 1.0, 0.2, 0.3, 
    0.3, 0.1, 0.5, 0.0,   0.4, 0.2, 1.0, 0.2, 
    0.2, 0.3, 0.2, 0.0,   -0.1, 0.3, 0.2, 1.0),
  nrow=8, ncol=8) 

df2<-as.data.frame(mvrnorm(n=N_samples, mu=mu, Sigma=sigma1))
df2$Female <- rbinom(n=N_samples, size=1, prob=0.5)
df2$Age <- runif(n=N_samples, min=50, max=85)

## check the data
head(df2)
cor(df2)

saveRDS(df2, file="sampleData2.RDS")



