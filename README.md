# sobermanclusters
soberman lab superresolution analyses tools

**UNBIASED CLUSTER AND LOCALIZATION DENSITY ANALYSIS FOR SINGLE MOLECULE LOCALIZATION MICROSCOPY**
---------------------------------------------------------------------

This software was used for the following publications: AB Schmider, HL Elliott, MD Godin, H Sunwoo, JT Lee, RJ Soberman. "Higher order organization of leukotriene synthetic proteins on the nuclear envelope" eLIFE (2018)

**REQUIREMENTS**

MATLAB 2013b or later
Toolboxes: Parallel Computing, Image Processing, Bioinformatics, Optimization, Signal Processing, Statistics and Machine Learning

**QUICK START**

Clone or download all files into the desired folder via link or through git clone. 
https://github.com/bairangie/sobermanclusters.git

Navigate to local cloned repository in MATLAB file path.  Add code files to path.
Begin by editing .txt files so that only columns A-R exist
For unbaised cluster analysis:
* STEP 1 - Define ROI
Execute by calling 'saveStormROI' at command and select folder with .txt file
Select ROI and save (now will have a .sroi file)
For two channels, a previous ROI from one channel and be applied to another channel.  Execute in command window saveStormROI('UseExistingROI',true)
In GUI window, select the .txt file to which you want to apply the ROI, then in subsequent GUI window select the .sroi file that has the ROI you want to use, save
* STEP 2 - Cluster Analysis
Execute by calling 'sobermanBatchClusterStormROIs'
Select folder of files to be run, in subsequent listSelectGUI window, highlight files and move to the right side window, hit ok (now will have clustering.mat file)  This step takes some time, depending on how many files you run and how many workers you have available.
* STEP 3 - Per-ROI Cluster Statistics
First, go to the completed files from STEP 2 and shorten the file names, they were duplicated and become too long.
Execute by calling 'batchStormClusterStats'
Select folder with clustering.mat files and run (now will have clustering_stats files in separate folders with 7 data figures for each file)
* STEP 4 - Per-Condition Combined 
make new folder for each condition and copy corresponding clustering.mat files into each folder.
Execute by calling 'combineStormClusterStats'
Select cluster_stats folder from Step 3, select same folder to save. (now will have combined_stats.mat file, 6 new data figures)
* To make normalized-point weighted histograms:
make subfolders for each condition using combined_stats.mat files (one file/folder) then put these folders(s) into one main folder
Execute by calling 'stormCombinedStatsAndConditionComparisonFigures('','condStr',{'type here what','ever','words','are','used','for','naming','folders'},'DataType','Point-weighted')

For localization density:
Execute by opening the LocalizationDensity.m file in the command window. (The readbinfile.m file is included in path)
Select .bin file to be analyzed.
Save
You can change the heat map percentage within the code for varied pseduocoloring results.












