# Sensitive information - not on Github

# Provide working directory
wd <- "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/SampleSizes/Incentive/Coherence/"

# Load the file with information about dataset
  # This file contains subject number without and with leading zero (1,2), 
  # then center ID (3) and the number I have of the subject in my folder.
# In addition, subjects that we don't need are removed from the R object. 
load('/Volumes/Elements/IMAGEN_DATA/INCENTIVE/BL/processed/OurIDsInt')
# Some columns are obsolote (refer to research questions earlier)
OurIDs <- FiltOurIDsInt[,-c(5:11)]

