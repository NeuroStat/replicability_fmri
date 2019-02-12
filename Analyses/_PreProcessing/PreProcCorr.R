####################
#### TITLE:     Unconditional mutual information between two images
#### Contents:
####
#### Source Files: /Volumes/2_TB_WD_Elements_10B8_Han/PhD/I****DATA/
####                Data/FreddieFreeloader/Script/Analyses/Faces/Unconditional
#### First Modified: 24/11/2015
#### Notes:
#################



##
###############
### Notes
###############
##

# Data comes from SplitSubjects
# We will calculate the correlation between t-maps with increasing sample size.


##
###############
### Preparation
###############
##


# Source paths
source('blind_PreProcessing.R')

# Possible contrasts: default = MATH > LANGUAGE
contrast <- c('ML', 'Faces', 'Incentive', 'StopGo')
contr <- contrast[4]

# Stability is plotted on separate figures for each contrast
# Therefor, we can use different objects
if(contr == 'ML'){
  dat <- RawDatSplit
}
if(contr == 'Faces'){
  dat <- RawDatSplitF
}
if(contr == 'Incentive'){
  dat <- RawDatSplitI
}
if(contr == 'StopGo'){
  dat <- RawDatSplitS
}


# Load in libraries
library(lattice)
library(gridExtra)
library(oro.nifti)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(Hmisc)
library(NeuRRoStat)

# Dimension of the brains
DIM <- c(53,63,46)
# Number of runs
NRUNS <- 50
# Number of steps
NSTEP <- 70

##
###############
### For loop: start with sample size S, calculate correlation between all images, then move to sample size S+1
###############
##

# Matrix were data will get into
MatrixCorrelation <- array(NA,dim=c(NSTEP,NRUNS))

# Print statement
PriST <- (c(1:NRUNS)/NRUNS)[seq(1,NRUNS,length.out=10)][-10]

# For loop over all runs
for(i in 1:NRUNS){
  IDprint <- c(i/NRUNS)==PriST
  if(any(IDprint)){
    print(paste(round(PriST,2)[IDprint]*100, "% done"))
  }
  # For loop over all steps
  for(j in 1:NSTEP){
    # We have two images (G1 and G2) from group 1 and group 2: both are the t-maps
    #-------------------------------------------------------------------------------------------------------------#
    imageG1 <- try(readNIfTI(paste(dat,'/Run_',i,'/Step_',j,'/Group1/stats','/tstat1.nii.gz',sep=''))[,,],silent=TRUE)
    if(!class(imageG1)=='array') next
    if(length(table(imageG1)) == 1) next
    maskG1 <- try(readNIfTI(paste(dat,'/Run_',i,'/Step_',j,'/Group1','/masked.nii.gz',sep=''))[,,],silent=TRUE)
    if(!class(maskG1)=='array') next
    
    #-------------------------------------------------------------------------------------------------------------#
    # Group 2
    imageG2 <- try(readNIfTI(paste(dat,'/Run_',i,'/Step_',j,'/Group2/stats','/tstat1.nii.gz',sep=''))[,,],silent=TRUE)
    if(!class(imageG2)=='array') next
    if(length(table(imageG1)) == 2) next
    maskG2 <- try(readNIfTI(paste(dat,'/Run_',i,'/Step_',j,'/Group2','/masked.nii.gz',sep=''))[,,],silent=TRUE)
    if(!class(maskG2)=='array') next
    
    #-------------------------------------------------------------------------------------------------------------#
    # Mask and put the masked values to NA (otherwise we have inflation of 0)
    idMaskG1 <- maskG1==0
    imageG1[idMaskG1] <- NA
    idMaskG2 <- maskG2==0
    imageG2[idMaskG2] <- NA
    
    #-------------------------------------------------------------------------------------------------------------#
    # Calculate correlation
    imageG1 <- array(imageG1,dim=prod(DIM))
    imageG2 <- array(imageG2,dim=prod(DIM))
    PearsCorr <- cor(imageG1,imageG2,method='pearson',use='complete.obs')
    
    #-------------------------------------------------------------------------------------------------------------#
    # Put correlation in matrix
    MatrixCorrelation[j,i] <- round(PearsCorr,6)
    
    #-------------------------------------------------------------------------------------------------------------#
    # Remove objects
    rm(imageG1,imageG2,PearsCorr, maskG1,maskG2)
  }
  if(i==NRUNS) print("100% done")
}

# Transform NaN values to 0
MatrixCorrelation[is.nan(MatrixCorrelation)] <- 0

# Check for missing values (TRUE = complete)
complete.cases(t(MatrixCorrelation))

##
###############
### Save intermediate results
###############
##

# Save the correlations
saveRDS(MatrixCorrelation, 
    file = paste(SaveLoc, '/', contr, '/MatrixCorrelation.rda',sep=''))


(i - 1)*70*2 + ((j - 1) * 2) + 1

fileCon <- paste(getwd(),"/runJobs.txt",sep="")
# Check data
for(i in 1:NRUNS){
  # For loop over all steps
  for(j in 1:NSTEP){
    ttestG1 <- try(readNIfTI(paste(dat,'Run_',i,'/Step_',j,'/Group1/thresh_zstat1.nii',sep=''))[,,],silent=TRUE)
    if(class(ttestG1) == 'try-error'){
      cat(paste((i - 1)*70*2 + ((j - 1) * 2) + 1, ' \n', sep = ''), file=fileCon, append = TRUE)
    }
    ttestG2 <- try(readNIfTI(paste(dat,'Run_',i,'/Step_',j,'/Group2/thresh_zstat1.nii',sep=''))[,,],silent=TRUE)
    if(class(ttestG2) == 'try-error'){
      cat(paste((i - 1)*70*2 + ((j - 1) * 2) + 2, ' \n', sep = ''), file=fileCon, append = TRUE)
    }
  }
}


##
###############
### Create plots
###############
##

quartz.options(width=18,height=12)

## Prepare matrix
Correlation.tmp <- matrix(MatrixCorrelation,ncol=1)
sampleSize <- rep(seq(10,700,by=10),NRUNS)
Correlation <- data.frame('PearsonCorr' = Correlation.tmp, 'size' = sampleSize)

# Variables for plotting
subjBreak <- c(seq(0,100,by=20),seq(100,700,by=50))

# Mean r at N=60
aggregate(PearsonCorr~size, data=Correlation,mean,na.rm=TRUE)


#################
## ggplot: points
#################
ggplot(Correlation, aes(x=factor(size),y=PearsonCorr)) +
  geom_point(colour='#006d2c',size=1.3) +
  scale_x_discrete(breaks=subjBreak, name="Sample size")+
  scale_y_continuous(name='Pearson Correlation')+
  theme(plot.title = element_text(lineheight=.2, face="bold")) +
  ggtitle('Correlation between t-maps of increasing sample size (N=42 per sample size).')



##
###############
### Plot for JSM 2016
###############
##

# Take average per sample size.
Avg <- aggregate(PearsonCorr ~ size, data = Correlation, FUN = mean)
Correlation$Average <- rep(Avg$PearsonCorr, NRUNS)

# Add horizontal line on 0.8290392 and vertical line on size = 60
ggplot(Correlation, aes(x = size,y=PearsonCorr)) +
  geom_point(colour='#377eb8',size=1.2) +
  #geom_segment(aes(x=60,xend=60,y=0,yend=0.8290392)) +
  #geom_segment(aes(x=0,xend=60,y=0.8290392,yend=0.8290392)) +
  geom_vline(xintercept=60, size = 1.2,linetype = 'longdash') +
  geom_hline(yintercept=0.8290392, size = 1.2,linetype = 'longdash') +
  geom_line(aes(x = size, y = Average), colour = '#d95f0e', size = 2) +
  scale_x_discrete(breaks=subjBreak, name="Sample size")+
  scale_y_continuous(breaks = seq(0,1,by = 0.2),name='Pearson correlation')+
  theme(axis.text=element_text(size=20, face="bold"),
        axis.title=element_text(size=24))








#################
## ggplot: boxplots
#################
ggplot(Correlation, aes(x=factor(sampleSize), y=PearsonCorr)) +
  geom_boxplot(fill='#006d2c') +
  scale_x_discrete(breaks=subjBreak, name="Sample size")+
  scale_y_continuous(name='Pearson Correlation')+
  theme(plot.title = element_text(lineheight=.2, face="bold")) +
  ggtitle('Correlation between t-maps of increasing sample size (N=42 per sample size).')





##
###############
### Calculate unconditional mutual information in amygdala
###############
##

## Created map for the Amygdala
map <- readNIfTI(mapAmygdala)[,,]

# Matrix were data will get into
MatCorrAmyg <- array(NA,dim=c(NSTEP,NRUNS))

# Print statement
PriST <- (c(1:NRUNS)/NRUNS)[seq(1,NRUNS,length.out=10)][-10]

# For loop over all runs
for(i in 1:NRUNS){
  IDprint <- c(i/NRUNS)==PriST
  if(any(IDprint)){
    print(paste(round(PriST,2)[IDprint]*100, "% done"))
  }
  # For loop over all steps
  for(j in 1:NSTEP){
    # We have two images (G1 and G2) from group 1 and group 2: both are the t-maps
    #-------------------------------------------------------------------------------------------------------------#
    imageG1 <- try(readNIfTI(paste(datawd,'/Run_',i,'/Step_',j,'/Group1/stats','/tstat1.nii.gz',sep=''))[,,],silent=TRUE)
    if(!class(imageG1)=='array') next
    maskG1 <- map
    
    #-------------------------------------------------------------------------------------------------------------#
    # Group 2
    imageG2 <- try(readNIfTI(paste(datawd,'/Run_',i,'/Step_',j,'/Group2/stats','/tstat1.nii.gz',sep=''))[,,],silent=TRUE)
    if(!class(imageG2)=='array') next
    maskG2 <- map
    
    #-------------------------------------------------------------------------------------------------------------#
    # Mask and put the masked values to NA (otherwise we have inflation of 0)
    idMaskG1 <- maskG1==0
    imageG1[idMaskG1] <- NA
    idMaskG2 <- maskG2==0
    imageG2[idMaskG2] <- NA
    
    #-------------------------------------------------------------------------------------------------------------#
    # Calculate correlation
    imageG1 <- array(imageG1,dim=prod(DIM))
    imageG2 <- array(imageG2,dim=prod(DIM))
    PearsCorr <- cor(imageG1,imageG2,method='pearson',use='complete.obs')
    
    #-------------------------------------------------------------------------------------------------------------#
    # Put correlation in matrix
    MatCorrAmyg[j,i] <- round(PearsCorr,6)
    
    #-------------------------------------------------------------------------------------------------------------#
    # Remove objects
    rm(imageG1,imageG2,PearsCorr,maskG1,maskG2)
  }
  if(i==NRUNS) print("100% done")
}


# Save objects: saveobjWD comes from blind_PreProcessing.R
save(MatCorrAmyg,file=paste(saveobjWD,'/MatrixCorrelationFaceAmygdala',sep=''))

# Or load them in
load(paste(saveobjWD,'/MatrixCorrelationFaceAmygdala',sep=''))


##
###############
### Create plots
###############
##



## Prepare matrix
CorrelationAm.tmp <- matrix(MatCorrAmyg,ncol=1)
sampleSize <- rep(seq(10,700,by=10),NRUNS)
CorrelationAm <- data.frame('PearsonCorr' = CorrelationAm.tmp, 'size' = sampleSize)

# Variables for plotting
subjBreak <- c(seq(0,100,by=20),seq(100,700,by=50))



#################
## ggplot: points
#################
ggplot(CorrelationAm, aes(x=factor(sampleSize),y=PearsonCorr)) +
  geom_point(colour='#006d2c',size=1.3) +
  scale_x_discrete(breaks=subjBreak, name="Sample size")+
  scale_y_continuous(name='Pearson Correlation')+
  theme(plot.title = element_text(lineheight=.2, face="bold")) +
  ggtitle('Pearson Correlation between t-maps of increasing sample size (N=42 per sample size) in Amygdala .')


#################
## ggplot: boxplots
#################
ggplot(CorrelationAm, aes(x=factor(sampleSize), y=PearsonCorr)) +
  geom_boxplot(fill='#006d2c') +
  scale_x_discrete(breaks=subjBreak, name="Sample size")+
  scale_y_continuous(name='Pearson Correlation')+
  theme(plot.title = element_text(lineheight=.2, face="bold")) +
  ggtitle('Correlation between t-maps of increasing sample size (N=42 per sample size) in Amygdala.')



# Now comparing amygdala with complete brain
bothCorr.tmp <- rbind(Correlation,CorrelationAm)
bothCorr <- data.frame(bothCorr.tmp, 'Structure' = rep(c(1,2),each=dim(Correlation)[1]))
bothCorr$Structure <- as.factor(bothCorr$Structure)

ggplot(bothCorr, aes(x=factor(size), y=PearsonCorr, fill=Structure)) +
  geom_boxplot() +
  scale_x_discrete(breaks=subjBreak, name="Sample size") +
  scale_fill_manual(values = c('#2b8cbe', '#016c59'), name='Brain Structure', labels = c('Whole brain', 'Amygdala')) +
  scale_y_continuous(name='Pearson Correlation') +
  theme(plot.title = element_text(lineheight=.2, face="bold")) +
  ggtitle('Correlation between t-maps of increasing sample size (N=42 per sample size) for whole brain and amygdala.')
