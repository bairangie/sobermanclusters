function clusterStormData(fileName,varargin)
%CLUSTERSTORMDATA bandwidth estimation, variable bandwidth mean shift clustering on 2D STORM file.
%
% clusterStormData
% clusterStormData(fileName)
% clusterStormData(fileName,'OptionName',optionValue,...)  (Options are described in input section of this function)
%   
%   This function performs per-point optimal bandwidth estimation and then
%   variable-bandwidth clustering on the selected STORM dataset. The
%   results are saved to a file and can then be explored using
%   stormROIClusterStats.m
% 
%Hunter Elliott
%5/2013

%% ------------- Input ---------- %%

ip = inputParser;
ip.addParamValue('CropData',false,@islogical);%If true, ask the user to select an ROI for points to be clustered.
ip.addParamValue('ROI',[],@(x)(isempty(x) || (size(x,2) == 2 && size(x,1) >= 3)));%Can directly input an ROI polygon.
ip.addParamValue('NumParallel',[],@isposint);%Number of matlab workers to use for parallelization (number of cores to use)
ip.addParamValue('HRange',[20 250],@(x)(numel(x) == 2 && diff(x) > 0));%Default is based on localization limit with primary+secondary antibody and mean inter-cluster distance
ip.addParamValue('nH',100,@isposint);%Number of bandwidths to try within specified range (will be log spaced)
ip.addParamValue('FramesUse',[],@(x)(isempty(x) || all(isposint(x))));%Frames to use data from. Empty means use all frames
ip.addParamValue('w',3,@isposint);%Number of bandwidths over which to calculate the Jensen-shannon divergence (measure of cluster sigma stability)
ip.addParamValue('OptimizedEstimation',true,@islogical);%If true, use the performance-optimized bandwidth estimation.
ip.addParamValue('OutputFile','',@ischar);%Name of file to write results to. Default is same as STORM file with _clustering.mat at end.
ip.KeepUnmatched = true;
ip.parse(varargin{:});

p = ip.Results;
op = ip.Unmatched;%For passing extra params.

outVars = {};%List of variables to write to file

%Get the filename here so we can store it in output
if nargin < 1 || isempty(fileName)
    [fileName, filepath] = uigetfile('*.txt','Select the storm localization file to open:');
    if filepath == 0        
        return
    end
    fileName = [filepath fileName];
elseif ~exist(fileName,'file')
    error(['"' fileName '" is not a valid file! Please specify a valid file to open!'])
end

if ~isempty(p.ROI)
    %Force cropping if an ROI polygon was input
    p.CropData = true;
end

outVars = [outVars 'fileName','p'];

%% ------------ Init ------------ %%



disp('Reading STORM data file...')
stormData = readSTORMTxtFile(fileName);
if isempty(stormData)
    %If the user clicked "cancel"
    return
end
nP = numel(stormData.X);
disp(['Loaded dataset with ' num2str(nP) ' points.'])

%Use only selected frames
if ~isempty(p.FramesUse)    
    stFields = fieldnames(stormData);    
    goodFrame = ismember(stormData.Frame,p.FramesUse);
    nFields = numel(stFields);
    for j = 1:nFields
        if numel(stormData.(stFields{j})) == nP            
            stormData.(stFields{j}) = stormData.(stFields{j})(goodFrame);
        end
    end
    nP = numel(stormData.X);
    disp(['Frame restriction limited dataset to ' num2str(nP) ' points.'])    
end

if p.CropData
   [stormData,cropPoly] = cropStormData(stormData,p.ROI);     %#ok<NASGU>
   outVars = [ outVars 'cropPoly'];%Save the ROI   
   nP = numel(stormData.X);
   disp(['Cropped dataset to ' num2str(nP) ' points.'])
end

outVars = [outVars 'stormData'];


X = [stormData.X stormData.Y];%Extract positions for analysis

if isempty(p.OutputFile)    
    %Default is to name after the storm data file
    p.OutputFile = [fileName(1:end-4) '_clustering.mat'];
end

outVars = [outVars 'p'];

%% ------------ Optimal Bandwidth Estimation ---------- %%

disp('Starting bandwidth estimation...')
tID = tic;
if p.OptimizedEstimation
    [Hest,ppSig,ppMu,ljsAll,Htry,HtryAll,isOutlier] = optimizedOptimalBandwidthEstimation(X,p.HRange,'nH',p.nH,'NumParallel',p.NumParallel,'w',p.w,'MinClustSize',3,op); %#ok<ASGLU>
    outVars = [outVars 'HtryAll'];
else
    [Hest,ppSig,ppMu,ljsAll,Htry] = estimateOptimalBandwidth(X,p.HRange,'nH',p.nH,'NumParallel',p.NumParallel,'w',p.w,op); %#ok<NASGU,ASGLU>
    isOutlier = false(size(X,1),1);
end    
totEstTime = toc(tID);
disp(['Finished bandwidth estimation, took ' num2str(totEstTime) ' seconds.'])

%For points where a sigma could not be estimated, we use the median.
noHEstimate = isnan(Hest);
if nnz(noHEstimate) > 0
    for j= 1:2
        for k = 1:2        
            Hest(j,k,isnan(Hest(j,k,:))) = nanmedian(Hest(j,k,:));
        end
    end
end
   
outVars = [outVars 'Hest','ppSig','ppMu','ljsAll' 'Htry','noHEstimate','isOutlier','totEstTime'];

%% ------ Variable Bandwidth clustering using estimated Bandwidths ---- %%

tID = tic;
disp('Starting variable bandwidth clustering...')
[clusterInfo, pointIndInlier, pointTrajInlier] = VariableBandwidthMeanShiftClustering(X(~isOutlier,:),Hest(:,:,~isOutlier),'method','standard','flagUseKDTree',1);
totClustTime = toc(tID);

%Pad the clustering output to include the outlier points
pointInd = zeros(nP,1);
pointInd(~isOutlier) = pointIndInlier; %#ok<NASGU>
pointTraj = cell(nP,1);
pointTraj(~isOutlier) = pointTrajInlier; %#ok<NASGU>
%And fix the indices in the cluster info to match the original data
iOrigInd = find(~isOutlier);
for j = 1:numel(clusterInfo)
    clusterInfo(j).ptIdData = iOrigInd(clusterInfo(j).ptIdData);
    
end

outVars = [outVars 'clusterInfo','pointInd','pointTraj','totClustTime'];

disp(['Finished variable bandwidth clustering, took ' num2str(totClustTime) ' seconds.'])

%% ----- Save to Disk ----- %%


save(p.OutputFile,outVars{:});



