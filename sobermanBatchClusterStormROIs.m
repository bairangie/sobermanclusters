function sobermanBatchClusterStormROIs(parentDir)


if nargin < 1 || isempty(parentDir)    
    parentDir = uigetdir(pwd,'Select a directory with the sroi files to analyze:');
end

allROIs = searchFiles('.sroi','',parentDir,1,'all',1);


selFiles = listSelectGUI(allROIs,'','move');

allROIs = allROIs(selFiles);
nFiles = numel(allROIs);


for iFile = 1:nFiles
        
    
    try
        
        
        %clusterStormROI(allROIs{iFile},'NumParallel',4,'FramesUse',1:16e3,'flagUseKDTree',1);
        clusterStormROI(allROIs{iFile},'NumParallel',6,'flagUseKDTree',1);
        
        
    catch err
        
        disp(['!!!Error on file ' num2str(iFile) ' :' err.message])
        save([allROIs{iFile}(1:end-5) '_errorReport.mat'],'err','allROIs','parentDir')
        
    end
    
end


