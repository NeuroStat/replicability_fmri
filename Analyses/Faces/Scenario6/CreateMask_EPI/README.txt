#### Creating the amygdala ROI in MNI space
Go to FSLVIEW, tools, atlases. Then select the regions you want.
For each selected region, go to save as, and save that region.
Then we add the regions that we saved and binarize the map according to a chosen probability (here larger then 10%). 
This is done through command line:
fslmaths amygdala1L -add amygdala1R -thr 10 -bin amygdala_thresh.nii.gz

	btw, I checked and this command line results in the same map as 
	fslmaths amygdala1L -add amygdala1R -thrp 10 -bin amygdala_thresh.nii.gz


#### Creation of amygdala ROI with IMAGEN data.
We flirt the MNI_brain T1_ 1 mm brain to the zstat, saving the transformation matrix. Then we use this matrix to flirt the amygdala_thresh.nii.gz image to the same zstat.
This Z-map comes from Freddie, analysis based on 700 subjects.


Command line code, using FSL:
flirt -in MNI152_T1_1mm_brain.nii.gz -ref ImagenEPI200_3mm.nii.gz -out flirtedMNI_EPI.nii -omat flirted_MNI_EPI_matrix

flirt -in amygdala_thresh.nii.gz -ref ImagenEPI200_3mm.nii.gz -init flirted_MNI_EPI_matrix -out flirted_amygdala_EPI.nii -omat flirted_MNI_EPI_matrix -applyxfm