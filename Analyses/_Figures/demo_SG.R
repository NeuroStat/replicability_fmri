####################
#### TITLE:     OLS group analysis of N = 20 subjects, STOP - GO task
#### Contents: 	
#### 
#### Source Files: /replicability_fmri/
#### First Modified: 20/06/2018
#### Notes: 
#################


# Take random 20 subjects from the stop and signal task.
# Pool the subjects using OLS.
# Do we find activation?

##
###############
### Preparation
###############
##

# Function to get p-value from t.test
getP <- function(x){
  pval <- try(t.test(x)$p.val, silent = TRUE)
  if(class(pval) == 'try-error'){
    pval <- NA
  }
  return(pval)
}

# Load in libraries
library(oro.nifti)
library(lattice)


##
###############
### Read in data
###############
##

# Empty map
con_map <- array(NA, dim = c(53,63,46,20))

# Location of data 
rawDat <- '/check_SG_data/'
for(i in 1:20){
  con_map[,,,i] <- readNIfTI(paste(getwd(), rawDat, i, '/EPI_stop_signal/con_0005_stop_success_-_stop_failure.nii.gz', sep = ''))[,,]
}

# OLS approach
grMap <- apply(X = con_map, MARGIN = c(1,2,3), FUN = getP)
hist(grMap)

# Uncorrected map
uncThresh <- array(NA, dim = c(53,63,46))
uncThresh[grMap <= 0.001] <- 1
uncThresh[grMap > 0.001] <- 0
pdf('GroupAnalysis_unc.pdf')
levelplot(uncThresh,
          main='Thresholded image (0/1): uncorrected p-values < 0.001',
          pretty=TRUE, col.regions = terrain.colors, xlab = 'X', ylab = 'Y', 
          xlim=c(0,53),ylim=c(0,63))
dev.off()

# Corrected
corMap <- array(p.adjust(p = c(grMap), method = 'fdr'), dim = c(53,63,46))
hist(corMap)

# Plot
threshMap <- array(NA, dim = c(53,63,46))
threshMap[corMap <= 0.05] <- 1
threshMap[corMap > 0.05] <- 0

pdf('GroupAnalysis_cor.pdf')
levelplot(threshMap,
          main='Thresholded image (0/1): FDR corrected p-values < 0.05',
          pretty=TRUE, col.regions = terrain.colors, xlab = 'X', ylab = 'Y', 
          xlim=c(0,53),ylim=c(0,63))
dev.off()

# Some checks
levelplot(con_map[,,,20],
          main='Subject con map',
          pretty=TRUE, col.regions = terrain.colors, xlab = 'X', ylab = 'Y', 
          xlim=c(0,53),ylim=c(0,63))

levelplot(grMap,
          main='Group P-value map',
          pretty=TRUE, col.regions = terrain.colors, xlab = 'X', ylab = 'Y', 
          xlim=c(0,53),ylim=c(0,63))
