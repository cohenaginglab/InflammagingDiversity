#####################################
## Code for multiple impuation of biomarker data that
##      is below the limit of detection (BLD)
##
## See Hargarten and Wheeler (2020). "Accounting for the uncertainty due to
##      chemicals below the detection limit in mixture analysis". Envrinomental Research 186.
##      https://doi.org/10.1016/j.envres.2020.109466
##
##
##    Method used in the manuscript:
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

library(miWQS)   ##for the imputation

#########################
## Load simulated data and add BLD data

## We have two simulated datasets of 8 biomarkers
## dat1 will be the core dataset (analogous to InCHIANTI in the manuscript)
## assume this is the one with the BLD data

dat1 <- readRDS("sampleData1.RDS")  

## We simulated this data without any BLD values so we add one more data series
dat1$newSeries <- rnorm(n=nrow(dat1), mean=10, sd=2)
## replace lowest 20% with detection limit to indicate it is BLD
detectionLimit <- quantile(dat1$newSeries, probs=0.2)
dat1$newSeries <- ifelse(dat1$newSeries < detectionLimit, detectionLimit, dat1$newSeries )
hist(dat1$newSeries, main="new series with BLD data")

#########################
## Multiple imputation of BLD data


# Step 1: Replace BLD values with NA
dat1_na <- dat1
dat1_na$newSeries <- ifelse(dat1_na$newSeries== detectionLimit, NA, dat1_na$newSeries)


# Step 2: Construct separate data frames of covariates and biomarkers
# We'll use age and sex as covariates in the imputation model
# (Note: in this example, they are not predictive of the BLD data)

Zcov <- dat1_na[,c("Female", "Age")]

# This is the matrix of biomarkers we will be imputing
Xbm <- dat1_na$newSeries 

# Step 3: How many datasets to impute?
NUM_MI <- 200

# Step 4: Impute
# set a seed so we can replicate the results
set.seed(401)
miDat <- impute.boot(X = Xbm, DL = detectionLimit, Z = Zcov,
                     K = NUM_MI, verbose = TRUE)

# miDat is a list of 3 -- the data is a 3-dimensional matrix in [[1]]
# check the data
(miDat$X.imputed[,,1])[1:12] # This is the first imputed dataset
(dat1_na$newSeries)[1:12]          # Only BLD values should be replaced with imputed values
(miDat$X.imputed[,,2])[1:12] # This is the second imputed dataset

# Step 5: 
# Put the data back together into NUM_MI complete imputed datasets and get NUM_MI covariance matrices
# The variables we did not impute need to be added back together:
mAdd <- as.data.frame(dat1[,1:8])
##initialise
cvm <- matrix(data=0, nrow=9, ncol=9)

for(i in 1:NUM_MI){  
  V9 <- miDat$X.imputed[,,i]
  m <- cbind(mAdd, V9) 
  ##standardize the imputed data 
  m$V9 <- scale(m$V9, center=T, scale=T)
  cvm <- cvm + cov(m)
}

# Now take the average of the covariance matrices
# This can be used in principal components or PCA analyses
cvm <- cvm/NUM_MI







