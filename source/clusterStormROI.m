function clusterStormROI(roiFile,varargin)
%CLUSTERSTORMROI runs clustering on the ROI specified in an SROI file created by saveStormROI
%
% clusterStormROI
% clusterStormROI(roiFile)
% clusterStormROI(roiFile,'ClusteringOptionName1',clusteringOptionValue1,...)
%
% Hunter Elliott
% 8/2014


if nargin < 1
    roiFile = '';
end

[roiFilePath,roiFile] = optionalFileInput(roiFile,'*.sroi','Select a storm ROI file to cluster:');

%Load the ROI info
roiDat = load([roiFilePath roiFile],'-mat');

%Get the storm file and path
stormFilePath = roiDat.filePath;
stormFile = roiDat.fileName;

fileFound = false;%Flag if we've found the (possibly relocated) file
fileRelocated = false;
if exist([stormFilePath stormFile],'file')
    fileFound = true;
    fileRelocated = false;
end

if ~fileFound
    %Check if the ROI file was moved, and see if the storm file was
    %similarly relocated
    stormFilePath = relocatePath(stormFilePath,roiDat.outFilePath,roiFilePath);
    stormFilePath = [stormFilePath filesep];        
    if exist([stormFilePath stormFile],'file')
        fileFound = true;
        fileRelocated = true;
        oldPath = roiDat.outFilePath;        
    end
end

if ~fileFound
    %Compare old STORM file path and new ROI file path to see if this
    %relocation explains it.
    stormFilePath = relocatePath(stormFilePath,roiDat.filePath,roiFilePath);
    stormFilePath = [stormFilePath filesep];        
    if exist([stormFilePath stormFile],'file')
        fileFound = true;
        fileRelocated = true;
        oldPath = roiDat.filePath;        
    end
end
 
if ~fileFound
    %grrrrr.... just check in the current directory of the ROI
    stormFilePath = roiFilePath;
    stormFilePath = [stormFilePath filesep];        
    if exist([stormFilePath stormFile],'file')
        fileFound = true;
        fileRelocated = true;
        oldPath = '??';        
    end
end

if ~fileFound
    %god damnit why are people so retarded? as a last resort, check the ROI
    %directory for similarly named .txt files, hoping they at least used
    %the default ROI naming convention..
    possStormFile = dir([roiFilePath filesep roiFile(1:min(strfind(roiFile,'_ROI')-1)) '*.txt']);
    if numel(possStormFile) == 1
        %Only if an unambiguous match was found..
        fileFound = true;
        fileRelocated = true;
        oldPath = '??';
        stormFilePath = roiFilePath;
        stormFile = possStormFile.name;
    end
end        

if ~fileFound
    error('Cannot find storm file for this ROI!!!!')
elseif fileRelocated
    warning('stormROIClust:relocation',['Detected relocation and or renaming!: ROI file and storm file were moved from ' oldPath ' to ' stormFilePath])
end
    
%Call the clustering function, passing additioinal inputs
clusterStormData([stormFilePath stormFile],'ROI',roiDat.cropPoly,'OutputFile',[stormFilePath filesep roiDat.fileName(1:end-4) '_ROI_' roiFile(1:end-5) '_clustering.mat'],varargin{:})
