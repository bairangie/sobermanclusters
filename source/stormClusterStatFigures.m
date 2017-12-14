function figHans = stormClusterStatFigures(stormData,clustStats,pointInd,Hest,varargin)




%% ------------ Input ---------- %%

ip = inputParser;
ip.addParamValue('SaveFigures',false,@islogical);
ip.addParamValue('OutputDirectory','',@ischar);
ip.addParamValue('PtSize',15,@isposint);
ip.addParamValue('SatPct',1,@isposint);
ip.parse(varargin{:});
p = ip.Results;

if nargin < 4
    Hest = [];
end


%% ---------- Init --------- %%

%Get the cooridnates of inliers only to ease indexing below.
XY = [stormData.X stormData.Y];
isOutlier = pointInd == 0;
XYout = XY(isOutlier,:);%But keep outliers so whe know what was ecxcluded
XY = XY(~isOutlier,:);
pointInd = pointInd(~isOutlier);
nP = numel(pointInd);
if ~isempty(Hest)
    Hest = Hest(:,:,~isOutlier);
end

figHans = [];
figNames = {};

if p.SaveFigures
    
    if isempty(p.OutputDirectory)
        p.OutputDirectory = uigetdir(pwd,'Select a directory to save figures to:');
    elseif ~exist(p.OutputDirectory,'dir')
        mkdir(p.OutputDirectory)
    end
    
end

satPct = [p.SatPct 100-p.SatPct];

%We use the same per-cluster coloring throughout
clustCols = randomColormap(clustStats.nClust,42);

%% ----- Make FIgures ---- %%

% ------ Cluster Assignment Figure ---- %

fH = fsFigure(.5);
hold on
scatter(XY(:,1),XY(:,2),p.PtSize,pointInd);
plot(XYout(:,1),XYout(:,2),'bx');
colormap(clustCols);
%arrayfun(@(x)(plot(clustStats.Center(x,1),clustStats.Center(x,2),'o','color',clustCols(x,:),'MarkerSize',15)),1:clustStats.nClust);
xlabel('X, nm')
ylabel('Y, nm')
title({['n = ' num2str(clustStats.nClust) ' clusters total'],'Color indicates cluster assignment, x indicates outlier (unclustered)'})
axis image
figHans = [figHans fH];
figNames = [ figNames 'cluster assignment plot'];

% ------ # Points Per-Cluster Figures ----- %

fH = fsFigure(.5);
hold on
sizeCols = jet(max(clustStats.NumPoints));    
sizePerPoint = clustStats.NumPoints(pointInd);%Size of associated clsuter for each point
scatter(XY(:,1),XY(:,2),p.PtSize,sizePerPoint);
colormap(sizeCols)
cRange = prctile(clustStats.NumPoints,satPct);
caxis(cRange);
colorbar
xlabel('X, nm')
ylabel('Y, nm')
title({['Mean # localizations per cluster = ' num2str(clustStats.MeanNumPoints)], ...
        'Color indicates # localizations for each cluster'})    
axis image
figHans = [figHans fH];
figNames = [ figNames 'number of points per cluster plot'];


fH = fsFigure(.5);
nlBins = 1:5:max(clustStats.NumPoints);
nlHist = histc(clustStats.NumPoints,nlBins);
bar(nlBins,nlHist)
xlabel('# Localizations')
ylabel('# Clusters')
xlim([min(nlBins) max(nlBins)])
title('Number of localizations per cluster histogram')
figHans = [figHans fH];   
figNames = [ figNames 'number of points per cluster histogram'];

% ------- Cluster Density Figure ------ %

fH = fsFigure(.5);
hold on
statCols = jet(512);    
statPerPoint = clustStats.Density(pointInd);%Size of associated clsuter for each point
scatter(XY(:,1),XY(:,2),p.PtSize,statPerPoint);
colormap(statCols)
cRange = prctile(clustStats.Density(clustStats.NumPoints>4),satPct);%The triangular and quadrilateral clusters have "Artificially" high density, so we ignore for visualization
caxis(cRange);
colorbar
xlabel('X, nm')
ylabel('Y, nm')
title({['Density per cluster (localizations / nm^2). Mean = ' num2str(nanmean(clustStats.Density))], ...
        'Color indicates Density for each cluster'})    
axis image
figHans = [figHans fH];
figNames = [ figNames 'density per cluster plot'];

% ------ Cluster Boundaries Figure ----- %

fH = fsFigure(.5);
hold on
hasBound = find(~cellfun(@isempty,clustStats.Boundary));
ptNoBound = ~any(bsxfun(@eq,hasBound',pointInd),2);
scatter(XY(ptNoBound,1),XY(ptNoBound,2),p.PtSize,pointInd(ptNoBound));
colormap(clustCols);
for j = hasBound(:)'
        
    patch(clustStats.Boundary{j}(:,1),clustStats.Boundary{j}(:,2),clustCols(j,:),'FaceAlpha',.2)    
    
end
xlabel('X, nm')
ylabel('Y, nm')
axis image
title({'CLuster convex hulls. Color indicates cluster number', 'Clusters with < 3 points shown as points'})
figHans = [figHans fH];
figNames = [ figNames 'cluster convex hulls plot'];


if ~isempty(Hest)
    % ------ Estimated Sigma Figure ------ %
    
    
    statPerPoint = nan(nP,1);
    for j = 1:nP
        statPerPoint(j) = mean(diag(Hest(:,:,j)));%Since we don't currently estimate off-diagonals, use mean of diag
    end

    fH=fsFigure(.5);
    hold on
    statCols = jet(512);        
    scatter(XY(:,1),XY(:,2),p.PtSize,statPerPoint);
    colormap(statCols)
    cRange = prctile(statPerPoint,satPct);%The triangular and quadrilateral clusters have "Artificially" high density, so we ignore for visualization
    caxis(cRange);
    colorbar
    xlabel('X, nm')
    ylabel('Y, nm')
    title({['Mean of Estimated XY Sigma per Point. All point mean = ' num2str(nanmean(statPerPoint))], ...
            'Color indicates sigma for each point'})    
    axis image
    figHans = [figHans fH];
    figNames = [ figNames 'estimated bandwidth per point plot'];
    
end

%% ------- Save to Disk --------- %%

if p.SaveFigures
    
    nFig = numel(figHans);
    
    for j = 1:nFig
        mfFigureExport(figHans(j),[p.OutputDirectory filesep figNames{j}]);
    end    
    
end
    
    






    
