function [H,ppSig,ppMu,ljsAll,Htry,HtryAll,isOutlier] = optimizedOptimalBandwidthEstimation(X,HRange,varargin)

%% ---------- Input ----------- %%

ip = inputParser;
ip.addRequired('X');
ip.addRequired('HRange',@(x)(numel(x) == 2 && diff(x) > 0 && all(x) > 0));%Minimum and maximum bandwidths to test for optimality
ip.addParamValue('DBSCANRadius',[],@(x)(isempty(x) || (x>0 && isscalar(x))));%Neighborhood radius to use in DBSCAN pre-processing. If not input, will be estimated from data.
ip.addParamValue('DBSCANNNPercentile',99,@(x)(x>0 && x <=100));%If DBSCAN radius is not input, it will be set at this percentile of the nearest-neighbor distances in the data. Yes, I know this is the most annoying parameter name ever created.
ip.addParamValue('nH',20,@isposint);%Total number of bandwidths to spread across the specified range.
ip.addParamValue('Verbose',true,@islogical);%If true, text progress/output displayed
ip.addParamValue('MinClustSize',2,@(x)(isequal(x,round(x)) &&isscalar(x)));%Minimum cluster size to use in DBSCAN pre-processing.
ip.addParamValue('ConstrainH',true,@islogical);%If true, H (bandwidth) ranges will be constrained based on the size of data sub-sets.
ip.addParamValue('RangePercentile',95,@(x)(x <= 100 && x > 0));%Bandwidth range for each sub-set is determined using the range of this percentile of the points in the sub-set. (If ConstrainH is enabled)

ip.KeepUnmatched = true;%For passing parameters to sub-functions.

ip.parse(X,HRange,varargin{:});
p = ip.Results;
op = ip.Unmatched;%Parameters for passing to estimation function

%% ------------ Init ------------ %%

nH = p.nH;%decreases annoyingness of code.
[n,d] = size(X);%Number of points and dimensionality.

%Setup the vector of possible bandwidths to try.
HtryAll = logspace(log10(HRange(1)),log10(HRange(2)),nH);

pTiles = [ (100-p.RangePercentile) / 2 100-(100-p.RangePercentile)/2];%Split specified percentile to high and low values
%Convert percentile to equivalent sigma so we can convert range to
%bandwidth
pTileSigma = abs(icdf('Normal',(100-p.RangePercentile)/200,0,1));

if isempty(p.DBSCANRadius)    
    
    %Estimate radius based on percentile of NN distances with specified
    %cluster size    
    if p.Verbose;tic;disp('Determinding DBSCAN Radius from NN distances...');end
    
    kd = KDTreeSearcher(X);
    [~,nnDist] = knnsearch(kd,X,'K',p.MinClustSize+1);%Self will be considered neighbor, so add 1 to min clust size
    
    nnDist = nnDist(:,2:end);
    p.DBSCANRadius = prctile(nnDist(:),p.DBSCANNNPercentile);
    
    if p.Verbose;disp(['Finished, took ' num2str(toc) ' seconds, selected DBSCAN radius of ' num2str(p.DBSCANRadius) '.']);end
end

%% ------- Pre-Processing ------- %%
%Sub-divides the dataset and removes isolated points using the DBSCAN
%algorithm.

if p.Verbose;tic;disp('Running DBSCAN Pre-processing...');end

cID = dbscan(X,p.MinClustSize,p.DBSCANRadius);

nSubSets = numel(unique(cID(cID~=0)));%Number of clusters excluding "outliers"/isolated points
isOutlier = cID == 0;
nOutliers = nnz(isOutlier);

if p.Verbose;disp(['Finished, took ' num2str(toc) ' seconds, partitioned data into ' num2str(nSubSets) ' sub-sets and excluded ' num2str(nOutliers) ' points as outliers.']);end

%% ----- Sub-divided Bandwidth Estimation ---- %%

HmaxAll = nan(nSubSets,1);
Htry = cell(nSubSets,1);
H = nan(2,2,n);
ppSig = nan(n,nH,2);
ppMu = nan(n,nH,2);
ljsAll = nan(n,nH);

if p.Verbose;tic;disp('Starting bandwidth estimation on sub-sets...');end

for j = 1:nSubSets
    
    
    ind = cID == j;%Get points in current sub-set to keep it simple
    Xcurr = X(ind,:);
    
    if p.ConstrainH
    
        %Get bandwidth range for this sub-set based on the range in the data
        %across all dimensions.
        allRange = nan(d,1);
        for k = 1:d
            allRange(k) = diff(prctile(Xcurr(:,k),pTiles));
        end
        %Use the smallest dimension to constrain the bandwidth since
        %clustering is done with scalar bandwidth.
        HmaxAll(j) = min(allRange) / (2*pTileSigma);%Divide by 2*sigma because we are using range
        %Use this maximum bandwidth to select from the specified test
        %bandwidths
        if HmaxAll(j) > HRange(2)
            Htry{j} = HtryAll;
        else
            Htry{j} = HtryAll(1:find(HtryAll>HmaxAll(j),1,'first'));       
        end
    else
        Htry{j} = HtryAll;
    end
    nHcurr = numel(Htry{j});    
        
    %Run the bandwidth estimation on this sub-set
    tic;
    op.HTry = Htry{j};%Add current bandwidths to structure for passing
    [H(:,:,ind),ppSig(ind,1:nHcurr,:),ppMu(ind,1:nHcurr,:),ljsAll(ind,1:nHcurr)] = estimateOptimalBandwidth(Xcurr,[],op);    
    
    if p.Verbose;disp(['Finished subset ' num2str(j) ', took ' num2str(toc) ' seconds']);end
           
end



if nnz(isnan(H(1,1,:))) > 0
    if p.Verbose;disp(['Unable to estimate bandwidth for ' num2str(nnz(isnan(H(1,1,:)))) ' points!']);end            
end




