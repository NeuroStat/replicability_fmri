# README

For some reason, the HPC had insufficient memory to read in the large masks (4D slices when N is large).
This caused the maskSlicing.R file to crash. This has no effect on the results as we only 
used this file to keep the first slice from the mask and convert it to 3D.
We still have the 4D masks, which are unfortunately bigger in size. 

The remove_cores.sh file is used on the HPC to remove the core.* files which were created
when the slicing failed. These are binary files created by HPC.