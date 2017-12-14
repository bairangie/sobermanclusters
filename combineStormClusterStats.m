function combStat = combineStormClusterStats(fileListOrDir,varargin)

%NOTE!!!! The MinArea and MinNumPoints filters will only affect the
%AllClust and CombinedClust fields and will NOT affect the AllDataSet and
%CombinedDataSet fields!!

%% -------- Input ------- %%

ip = inputParser;

ip.addParamValue('Verbose',true,@islogical);%If true, text showing progress will be displayed
ip.addParamValue('MinNumPoints',3,@(x)(isposint(x)));%Minimum number of points in a cluster for it to be included
ip.addParamValue('MinArea',eps,@isscalar);%Minimum area of clusters to include
ip.addParamValue('NumBinsClust',50,@(x)(isposint(x)));%Number of bins to use in per-cluster combined histograms
ip.addParamValue('NumBinsDataset',10,@(x)(isposint(x)));%Number of bins to use in per-cluster combined histograms
ip.addParamValue('OutputDirectory','',@ischar);%Directory to save combined data and figures to 

ip.parse(varargin{:});

p = ip.Results;


%% -------- Init ------- %%

if nargin < 1 || isempty(fileListOrDir)
    fileListOrDir = uigetdir(pwd,'Select a parent directory with ROI stats:');
end

if ~iscell(fileListOrDir);
    fileList = findFilesInSubDirs(fileListOrDir,'roi cluster stats.mat');    
else
    fileList = fileListOrDir;
end


nFiles = numel(fileList);

combFun = {@nanmean,@nanstd,@nanmedian};
combFunNames = {'Mean','STD','Median'};
nCombFun = numel(combFun);


if isempty(p.OutputDirectory)    
    p.OutputDirectory =uigetdir(pwd,'Select a directory for saving combined results:');    
elseif ~exist(p.OutputDirectory,'dir')
    error('Invalid output directory!')
end

%% ------- Loading ------ %%
if p.Verbose;tic;disp(['Loading ' num2str(nFiles) ' stat files...']);end

allInfo = cell(nFiles,1);

for j = 1:nFiles
    
    allInfo{j} = load(fileList{j});    
    
    if ~isfield(allInfo{j},'clustStats')
        error(['Invalid file: ' fileList{j}]);
    end
                
end

allInfo = vertcat(allInfo{:});
allStat = [allInfo.clustStats];

if p.Verbose;disp(['Finished, took ' num2str(toc) ' seconds.']);end    
%% ----- Processing ------- %%

if p.Verbose;tic;disp('Combining all per-cluster statistics...');end

%TEMP - auto-find fields to combine based on size and numclust!!
perClustFields = {'NumPoints','Area','Density'};
nPCFields = numel(perClustFields);


for j = 1:nPCFields   
    combStat.(['AllClust' perClustFields{j}]) = vertcat(allStat.(perClustFields{j}));    
end

%Filter the clusters using the input criteria
usePoints = combStat.AllClustNumPoints >= p.MinNumPoints & combStat.AllClustArea >= p.MinArea;%% TEMP - make this generic?

%We loop twice so we can use stats to filter points (e.g. numPoints)
for j = 1:nPCFields      
    %Remove filtered points
    combStat.(['AllClust' perClustFields{j}]) = combStat.(['AllClust' perClustFields{j}])(usePoints);
    for k = 1:nCombFun
        combStat.(['CombinedClust' combFunNames{k} perClustFields{j}]) = combFun{k}(combStat.(['AllClust' perClustFields{j}]));   
    end        
end

if p.Verbose;disp(['Finished, took ' num2str(toc) ' seconds.']);end    

if p.Verbose;tic;disp('Combining all per-datset statistics...');end

perSetFields = {'MeanDistToCent','MeanNumPoints','nClust'};
nPSFields = numel(perSetFields);

for j = 1:nPSFields
    
    combStat.(['AllDataSet' perSetFields{j}]) = vertcat(allStat.(perSetFields{j}));
    
    for k = 1:nCombFun        
        combStat.(['CombinedDataSet' combFunNames{k} perSetFields{j}]) = combFun{k}(combStat.(['AllDataSet' perSetFields{j}]));
    end    
end

%Store these for posterity
combStat.NumFiles = nFiles;
combStat.PointsUsed = usePoints;

if p.Verbose;disp(['Finished, took ' num2str(toc) ' seconds.']);end    


%% ----- Per-Cluster Stat Figures ---- %%
if p.Verbose;tic;disp('Making per-cluster stat figures...');end

for j = 1:nPCFields
        
    currFig = figure;
    hist(combStat.(['AllClust' perClustFields{j}]),p.NumBinsClust);
    xlabel(['Per-Cluster ' perClustFields{j} ', All Datasets'])
    ylabel('#  Clusters')
    mfFigureExport(currFig,[p.OutputDirectory filesep 'Per_Cluster_' perClustFields{j}]);                

end

if p.Verbose;disp(['Finished, took ' num2str(toc) ' seconds.']);end  

%% ----- Per-Dataset Stat Figures ---- %%
if p.Verbose;tic;disp('Making per-dataset stat figures...');end

for j = 1:nPSFields
        
    currFig = figure;
    hist(combStat.(['AllDataSet' perSetFields{j}]),p.NumBinsDataset);
    xlabel(['Per-Dataset ' perSetFields{j} ', All Datasets'])
    ylabel('#  Datasets')
    mfFigureExport(currFig,[p.OutputDirectory filesep 'Per_Dataset_' perSetFields{j}]);                

end

if p.Verbose;disp(['Finished, took ' num2str(toc) ' seconds.']);end  

%% ---------- Output --------- %%


outFile = [p.OutputDirectory filesep 'combined_stats.mat'];

save(outFile,'combStat','p','fileList');

