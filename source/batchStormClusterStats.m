function batchStormClusterStats(dirOrFileList,varargin)
%

%% ------------- Input ----------- %%

ip = inputParser;

ip.addParamValue('Verbose',true,@islogical);%If true, text showing progress will be displayed
ip.addParamValue('CloseFigs',true,@islogical);%If true, figures are closed after each dataset

ip.parse(varargin{:});

p = ip.Results;



if nargin < 1 || isempty(dirOrFileList)    
    dirOrFileList = uigetdir(pwd,'Select a directory containig clustering files:');    
end

if ~iscell(dirOrFileList)
    
    if exist(dirOrFileList,'dir')
        fileList = dir([dirOrFileList filesep '*clustering.mat']);
        fileList = {fileList.name};
        fileList = cellfun(@(x)([dirOrFileList filesep x]),fileList,'Unif',0);
    else
        error('First input must be a valid directory or a cell array of file names!')
    end    
    
else
    fileList = dirOrFileList;
end


%% ----------- Init ----------- %%

nFiles = numel(fileList);

if nFiles == 0
    error('No valid storm clustering files specified!')
end


%% ---------- Processing ---------- %%



for iFile = 1:nFiles
    
    if p.Verbose;tic;disp(['Processing file ' num2str(iFile) ' of ' num2str(nFiles) '...']);end
    
    
    
    currOutDir = [fileList{iFile}(1:max(strfind(fileList{iFile},'.'))-1) '_stats'];
    if ~exist(currOutDir,'dir')
        mkdir(currOutDir);
    end
    
    try
        %Load parameters used for this file
        filePar = load(fileList{iFile},'p');filePar = filePar.p;
        %Check for ROI so the area can be passed if present
        if filePar.CropData
            
            if p.Verbose;disp(['Found ROI, applying.']);end
            
            %Load the ROI polygon
            cropPoly = load(fileList{iFile},'cropPoly');cropPoly = cropPoly.cropPoly;
            wholeIm = false;
        else
            if p.Verbose;disp(['No ROI found, calculating whole-image stats.']);end
            cropPoly = [];
            wholeIm = true;            
        end    
        
        %Run the stats on this file.
        stormROIClusterStats(fileList{iFile},'OutputDirectory',currOutDir,'ROIPoly',cropPoly,'WholeImage',wholeIm);
    catch err
        disp(['!!! Error calculating stats for file ' num2str(iFile) ' : ' err.message]);
    end
    
    if p.CloseFigs
        close all
    end
    
    
    if p.Verbose;disp(['Finished, took ' num2str(toc) ' seconds.']);end    
    
end


