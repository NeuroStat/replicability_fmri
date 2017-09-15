

IMGtmp <- readNIfTI('/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Data/SplitSubjects/ScenarioA/Run_5/Step_35/Group1/stats/zstat1.nii.gz')[,,]
mas <- readNIfTI('/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Data/SplitSubjects/ScenarioA/Run_5/Step_35/Group1/masked.nii.gz')[,,]


	IMGtmp <- round(IMGtmp,6)
	IMGtmp <- matrix(IMGtmp,ncol=1)
	mas <- matrix(mas,ncol=1)
	Masked <- matrix(IMGtmp[mas==0],ncol=1)
IMGtmp[mas==0] <- apply(Masked,1,function(x){x <- x+runif(n=1,min=-17,max=17)})
tval <- (sign(IMGtmp)*(-1))*qt(pnorm(-abs(IMGtmp)),df=df)
tval <- array(tval,dim=DIM)	

tval[mas==0] <- NA
img[mas==0] <- NA
tval <- round(tval,4)
tval[mas==0] <- 255
apply(matrix(tval[mas==0],ncol=1),1,runif,n=1,min=-0.5,max=0.5)

convert.datatype(datatype.code='FLOAT32')
tets <- nifti(tval,dim=DIM,datatype=16)
datatype(tets)
writeNIfTI(tets,filename='/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Data/SplitSubjects/ScenarioA/Run_5/Step_35/Group1/stats/tstat1',gzipped=FALSE)
writeNIfTI(img.nifti,filename='/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Data/SplitSubjects/ScenarioA/Run_5/Step_35/Group1/stats/tstat1',gzipped=FALSE)

norm <- dnorm(seq(-5, 5, length=32), sd=2)
norm <- (norm-min(norm)) / max(norm-min(norm))
img <- outer(outer(norm, norm), norm)
img <- round(255 * img)
img[17:32,,] <- 255 - img[17:32,,]
img.nifti <- nifti(img)
save(tval)

b <- readNIfTI('/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Data/SplitSubjects/ScenarioA/Run_5/Step_35/Group1/stats/tstat1.nii')[,,]
IMGtmp
a <- array(rnorm(prod(DIM),0,5),dim=DIM)


img <- nifti(a,dim=DIM)
writeNIfTI(a,filename='/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Data/SplitSubjects/ScenarioA/Run_5/Step_35/Group1/stats/tstat1',gzipped=FALSE)






##
###############
### For loop: start with group A, sample size S within a RUN, calculate percentage, adapt if necessary, 
### move to group B, calculate overlap between the two images, then move to sample size S+1
###############
##

# Vector with degrees of freedom
DF <- c(seq(9,699,by=10))

# Matrix were overlap results will get into
MatrixOverlapAdap <- array(NA,dim=c(NSTEP,NRUNS))

# Vector where we will store the final percentage of activated masked voxels
	# As we have two images per step and run, we will take the average of the percentage. 
PercActAdap <- array(NA,dim=c(NSTEP,NRUNS))

# Vector to store average of percentage activated masked voxels in each step and run
SignLevels <- array(NA,dim=c(NSTEP,NRUNS))

# Vector for activation map of first run only
SPMRun1G1 <- array(NA,dim=c(prod(DIM),NSTEP))
SPMRun1G2 <- array(NA,dim=c(prod(DIM),NSTEP))


# For loop over all runs
for(i in 1:NRUNS){
	print(paste('--------------- Run ', i,'--------------',sep=''))
	# For loop over all steps
	for(j in 1:NSTEP){
		print(paste('Step ',j,sep=''))
		###################################################################
		# We have two images (G1 and G2) from group 1 and group 2: 0 is non significant and 1 is significant
			# We also need to check whether both of the thresholded maps are present (if not, skip step)
		imageG1 <- try(readNIfTI(paste(dat,'/Run_',i,'/Step_',j,'/Group1','/stats/tstat1.nii',sep=''))[,,],silent=TRUE)
			if(!class(imageG1)=='array') next
			maskG1 <- try(readNIfTI(paste(dat,'/Run_',i,'/Step_',j,'/Group1','/masked.nii.gz',sep=''))[,,],silent=TRUE)
				if(!class(maskG1)=='array') next
			idMaskG1 <- maskG1==0
			imageG1[idMaskG1] <- NA
			# Check second image (there is no point in calculating percentages and stuff if not both images are available)
		imageG2 <- try(readNIfTI(paste(dat,'/Run_',i,'/Step_',j,'/Group2','/stats/tstat1.nii',sep=''))[,,],silent=TRUE)
			if(!class(imageG2)=='array') next
			maskG2 <- try(readNIfTI(paste(dat,'/Run_',i,'/Step_',j,'/Group2','/masked.nii.gz',sep=''))[,,],silent=TRUE)
				if(!class(maskG2)=='array') next
			idMaskG2 <- maskG2==0
			imageG2[idMaskG2] <- NA

		###################################################################
		# Now that we have t-statistic images with its masks, calculate BF-adjusted p-values using baseline significance level
			# Degrees of freedom
			df <- DF[j]
		pvalsG1 <- 1-pt(imageG1,df)
		pvalsG2 <- 1-pt(imageG2,df)
	
		# Function AdapThresholding: calculates adjusted p-values, then finds correct significance level according to percentage.
		adapG1 <- AdapThresholding(pvalsG1,DIM,signLevel,idMaskG1)
			signLevel <- adapG1$signLevel
			perc <- adapG1$percentage

		# Image G2
		adapG2 <- AdapThresholding(pvalsG2,DIM,signLevel,idMaskG2)
			signLevel <- round(mean(c(signLevel,adapG2$signLevel)),6)
			perc <- round(mean(c(perc,adapG2$percentage)),4)


		###################################################################
		# Calculate overlap: summing image K and K-1 to know the voxels in both maps (1+1 = 2)
		sumMap <- adapG1$SPM+adapG2$SPM
		# Minus map: know the voxels different in both images
		minusMap <- adapG1$SPM-adapG2$SPM
			# Mask this image
		idMask <- maskG1*maskG2
			minusMap[idMask] <- NA
	
		# Union of activated voxels, voxels in image G1 and voxels in image G2
		Vjt <- length(sumMap[which(sumMap==2)])
		Vt <- length(adapG1$SPM[which(adapG1$SPM==1)])
		Vj <- length(adapG2$SPM[which(adapG2$SPM==1)])
		# Voxels different from both images
		VjtS <- length(minusMap[which(minusMap==0)])


		# Put overlap, percentage and signLevel in matrix
		MatrixOverlapAdap[j,i] <- round((Vjt)/(Vj + Vt - Vjt),6)
		PercActAdap[j,i] <- perc
		SignLevels[j,i] <- signLevel



		###################################################################
		# From run 1: take the thresholded maps for levelplotting
		if(i==1){
			assign(paste('ThrMap1_',j,sep=''),
			levelplot(adapG1$SPM,main=paste('Group 1 analysis of step ' ,j,sep=''),
			pretty=TRUE, col.regions = terrain.colors, xlab = 'X', ylab = 'Y', 
			xlim=c(0,DIM[1]),ylim=c(0,DIM[2]))
			)
			assign(paste('ThrMap2_',j,sep=''),
			levelplot(adapG2$SPM,main=paste('Group 2 analysis of step ' ,j,sep=''),
			pretty=TRUE, col.regions = terrain.colors, xlab = 'X', ylab = 'Y', 
			xlim=c(0,DIM[1]),ylim=c(0,DIM[2]))
			)
			# Activation maps
			SPMRun1G1[,j] <- adapG1$SPM
			SPMRun1G2[,j] <- adapG2$SPM
		}


		# Remove objects
		rm(Vjt,Vt,Vj,imageG1,imageG2,sumMap,maskG1,idMaskG1,maskG2,idMaskG2,adapG1,adapG2)
		gc()
	}
if(i==NRUNS) print("100% done")
}

