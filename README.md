# sobermanclusters
soberman lab superresolution tools

UNBIASED CLUSTER ANALYSIS FOR SINGLE MOLECULE LOCALIZATION MICROSCOPY

This software was used for the following publications: AB Schmider, HL Elliott, MD Godin, H Sunwoo, JT Lee, RJ Soberman. "Higher order organization of leukotriene synthetic proteins on the nuclear envelope" eLIFE (2018)

REQUIREMENTS

MATLAB 2013b or later
Distributed computing, Image Processing and statistical analyses toolboxes

QUICK START

Clone or download all files into the desired folder via link or through git clone
https://github.com/Sobermanclusters.git

Navigate to local cloned repository in MATLAB file path
Begin by editing .txt files so that only columns A-R exist 
STEP 1 - Define ROI
Execute by calling 'saveStormROI' at command and select folder with .txt file
Select ROI and save (now will have a .sroi file)
For two channels, a previous ROI from one channel and be applied to another channel.  Execute in command window saveStormROI('UseExistingROI',true)
In GUI window, select the .txt file to which you want to apply the ROI, then in subsequent GUI window select the .sroi file that has the ROI you want to use, save
STEP 2 - Cluster Analysis
Execute by calling 'sobermanBatchClusterStormROIs'
Select folder of files to be run, in subsequent listSelectGUI window, highlight files and move to the right side window, hit ok (now will have clustering.mat file)







