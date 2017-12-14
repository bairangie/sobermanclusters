function stormROIClusterStats(fileName,varargin)
%STORMROICLUSTERSTATS produces figures and stats for an ROI from an already-clustered STORM dataset
%
% stormROIClusterStats
% stormROIClusterStats(fileName);
% stormROIClusterStats(fileName,'OutputDirectory','output directory')
%
%   This function requires that the STORM dataset has already been
%   clustered using the clusterStormData function. The user then selects a
%   region of interest (ROI) and several figures are produced visualizing
%   the clustering and various cluster statistics (size,
%   density etc.). The figures are saved to the output directory, and the
%   cluster statistics are also written to a .csv spreadsheet in the same
%   directory.
%
%
%Hunter Elliott
%5/2013

%% ---- Parameters ----- %%

%Cluster stat fields to write to spreadsheet
statFields = {    
    'NumPoints'    
    'Area'
    'Density'};

nFields = numel(statFields);

outFile = 'roi cluster stats.mat';

%% ------------- Input -------------- %%

ip = inputParser;
ip.addParamValue('OutputDirectory','',@ischar);
ip.addParamValue('WholeImage',false,@islogical);%If true, no ROI is selected and stats are generated for the entire image.
ip.addParamValue('ROIPoly',[],@(x)(isempty(x) || (size(x,1) > 2 && size(x,2) == 2)));%Can directly input an ROI polygon. If not input and WholeImage is false, user will be asked to select an ROI interactively.
ip.parse(varargin{:});

p = ip.Results;



%Get the filename here so we can store it in output
if nargin < 1
    fileName = '';
end

[filePath,fileName] = optionalFileInput(fileName,'*clustering.mat','Select the storm clustering file to open:');


if isempty(p.OutputDirectory)
    
     p.OutputDirectory = uigetdir(filePath,'Select a directory to save results to:');
    if p.OutputDirectory == 0
        return
    end    
    
end

%% -------------- Init ---------- %%


if ~exist(p.OutputDirectory,'dir')
    mkdir(p.OutputDirectory);
end

clustData = load([filePath fileName]);
nCTot = numel(clustData.clusterInfo);

if ~p.WholeImage
    %Let user select region to analyze
    [stormData,roiPoly,isInROI] = cropStormData(clustData.stormData,p.ROIPoly);    
else
    stormData = clustData.stormData;
    %Simplify processing below by creating an artificial ROI that includes
    %all points
    roiPoly = convhull(stormData.X,stormData.Y);
    roiPoly = [stormData.X(roiPoly) stormData.Y(roiPoly)];
    isInROI = true(numel(stormData.X),1);
end
%Get area of ROI
roiArea = polyarea(roiPoly(:,1),roiPoly(:,2));
    
%Make and save ROI figure
currFig = fsFigure(.5);
if any(~isInROI)    
    plot(clustData.stormData.X(~isInROI),clustData.stormData.Y(~isInROI),'.')
    legStr = {'Unselected Localizations','Selected Localizations','ROI'};
else
    legStr = {'Selected Localizations','ROI'};
end
hold on
xlabel('X, nm')
ylabel('Y, nm')
xlim auto,ylim auto;%In case the ROI includes all points, this fixes the limits
plot(clustData.stormData.X(isInROI),clustData.stormData.Y(isInROI),'r.')
patch(roiPoly(:,1),roiPoly(:,2),'r','FaceAlpha',.3)
axis image
legend(legStr)
title(['Area of selected ROI: ' num2str(roiArea) ' nm^2'])
figName = [p.OutputDirectory filesep 'Selected ROI area'];
mfFigureExport(currFig,figName);


%Fix the cluster and indices given the crop
clustIsIn = false(nCTot,1);
iInROI = find(isInROI);
n = numel(iInROI);
ptMapVec = zeros(n,1);
ptMapVec(iInROI) = 1:n;%Vector for converting indices
for j = 1:nCTot
    
    %Check if any of the points in this cluster were selected
    if any(isInROI(clustData.clusterInfo(j).ptIdData))
       clustIsIn(j) = true;
       %Remove the points which were cropped
       clustData.clusterInfo(j).ptIdData(~isInROI(clustData.clusterInfo(j).ptIdData)) = [];
       clustData.clusterInfo(j).ptIdData = ptMapVec(clustData.clusterInfo(j).ptIdData);
       clustData.clusterInfo(j).numPoints = numel(clustData.clusterInfo(j).ptIdData);    
    end
end

clusterInfo = clustData.clusterInfo(clustIsIn);

nC = numel(clusterInfo);

iClustIn = find(clustIsIn);
clMapVec = zeros(max(iClustIn),1);
clMapVec(iClustIn) = 1:nC;%Vector for converting indices
pointInd = zeros(n,1);
inInd = clustData.pointInd(isInROI);
pointInd(inInd~=0) = clMapVec(inInd(inInd~=0));

Hest = clustData.Hest(:,:,isInROI);

%% ---------- Get stats for this ROI ------ %%

clustStats = stormClusterStats(stormData,clusterInfo,pointInd);

stormClusterStatFigures(stormData,clustStats,pointInd,Hest,'SaveFigures',true,'OutputDirectory',p.OutputDirectory);


%% --------- Write to Disk ---- %%

for j = 1:nFields
    
   currFileName = [p.OutputDirectory filesep statFields{j} '.csv'];   
   csvwrite(currFileName,clustStats.(statFields{j}))    
    
end

save([p.OutputDirectory filesep outFile],'roiPoly','stormData','iInROI','clustStats','fileName','filePath','iClustIn','pointInd','Hest','clusterInfo','nC','roiArea');

