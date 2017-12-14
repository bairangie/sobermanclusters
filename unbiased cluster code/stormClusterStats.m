function clustStats = stormClusterStats(stDat,clustInf,pointInd,varargin)


%% ---------------- Input ------------ %%

% ip = inputParser;
% 
% ip.parse(varargin{:});
% p = ip.Results;


%% ------------- Init -------------- %%

nP = numel(stDat.X);
XY = [stDat.X stDat.Y];


%% --------------- Clustering ----------- %%

clustStats.nClust = numel(clustInf);
clustStats.Center = vertcat(clustInf.ptClusterCenter);
clustStats.NumPoints = vertcat(clustInf.numPoints);



%% ------- Cluster Stats --------------- %%


%Get distance from each point to cluster center
clustCentMat = nan(nP,2);
clustCentMat(pointInd~=0,:) = clustStats.Center(pointInd(pointInd~=0),:);
dToCent = sqrt(sum((clustCentMat - XY) .^2,2));
%And the mean distance to the cluster center for each cluster
clustStats.MeanDistToCent = arrayfun(@(x)(mean(dToCent(pointInd==x))),1:clustStats.nClust)';

clustStats.MeanDistToCent = mean(clustStats.MeanDistToCent);
clustStats.MeanNumPoints = mean(clustStats.NumPoints);
clustStats.STDNumPoints = std(clustStats.NumPoints);

%Get per-cluster density
clustStats.Boundary = cell(clustStats.nClust,1);
clustStats.Area = nan(clustStats.nClust,1);
clustStats.Density = nan(clustStats.nClust,1);
for j = 1:clustStats.nClust    
    if clustStats.NumPoints(j) >= 3    
        chInd = convhull(XY(clustInf(j).ptIdData,1),XY(clustInf(j).ptIdData,2));
        clustStats.Boundary{j} = XY(clustInf(j).ptIdData(chInd),:);
        clustStats.Area(j) = polyarea(clustStats.Boundary{j}(:,1),clustStats.Boundary{j}(:,2));
        clustStats.Density(j) = clustStats.NumPoints(j)   / clustStats.Area(j);
    end
end


