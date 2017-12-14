function stormCombinedStatsAndConditionComparisonFigures(allCondDir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------ PARAMS------ %%
%%%%%%%%%%%%%%%%%%%%%%%%%

statFieldNames = {'AllClustNumPoints','AllClustArea','AllClustDensity'};
statNames = {'Cluster No. Localizations','Cluster Area, nm^2','Cluster Density, Localizations/nm^2'};

%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------ INPUT ------ %%
%%%%%%%%%%%%%%%%%%%%%%%%%


ip = inputParser;
ip.addParameter('condStr',{}); %Optional. Can specify a string to search the combined_statistics file name for determining the condition and creating figure legends.
ip.addParameter('condNames',{});%Optional - if strings for file matching are different than those to use in figures, specify here (cell array of same size as condStr)
ip.addParameter('OutputFile',[]);%If not input, user will be asked - folder to store output in.
ip.addParameter('DataType','Raw',@(x)(ismember(x,{'Raw','Point-Weighted','Per-ROI-Mean','Per-ROI-Mean-Point-Weighted','Per-ROI-NumClusters','Per-ROI-NumLocalizations'})));%Type of data to use in figure. Median is only compatible with bar-graphs
ip.addParameter('StatToCalc','',@(x)(ismember(x,{'','Median'})));%Optional statistic of data to calculate and compare
ip.addParameter('HistNorm',true,@islogical);%If true, normalized histograms (frequency) will be displayed, if false raw counts.
ip.addParameter('FigureType','FilledHistogram',@(x)(ismember(x,{'Bar','Box','BarHistogram','FilledHistogram','Scatter'})));%Type of figure / comparison to display.
ip.addParameter('PlotArgs',{});%Optional - additional arguments to pass to plotting function.
ip.addParameter('iStatFile',[],@(x)(all(isposint(x))));%If not input, user will be asked - index of file(s) to produce figures/comparisons from.
ip.addParameter('DataFilterFun',[]);%Optionally input a function to use to filter data before plotting stats, e.g. to select only large clusters etc. The function should produce a boolean vector with length corresponding to the number of relevant data elements (e.g. clusters), where the true elements will be retained in the figure/stats.
ip.addParameter('iStatField',[],@(x)(all(isposint(x))));%If not input, user will be asked - index of cluster statistic to use in figures/comparisons, 1=# points, 2=Area, 3=density. Specify one number for all plot types except scatter
ip.addParameter('nBins',25,@isposint);%Number of bins to use in histograms
ip.addParameter('FilledHistType','Stairs',@(x)(ismember(x,{'Line','Stairs'})));%Type of curves to plot for filled histogram figures.
ip.addParameter('PctLimHist',[0 95]);%Percentiles to set limits of histograms at.
ip.addParameter('PctLimPlot',[0.5 99.5]);%Percentiles to set limits of plots at.

ip.parse(varargin{:});
p = ip.Results;

if nargin < 1
    allCondDir = '';
end

if isempty(allCondDir)
    allCondDir = uigetdir(pwd,'Select the directory containing the combined statistics files:');
    if allCondDir == 0,return;end        
end
   
if isempty(p.iStatField)
    [iStatField,OK] = listdlg('Name','Select a statistic:','ListString',statNames,'ListSize',[300 100]);    
    if ~OK,return;end    
else
    iStatField = p.iStatField;
end

if ~isempty(p.DataFilterFun) && ~ismember(p.DataType,{'Per-ROI-Mean','Per-ROI-NumClusters','Per-ROI-NumLocalizations'})
    error('This data type does not support data filtering yet!')
end
    

%%%%%%%%%%%%%%%%%%%%%%%%
%% ------ INIT ------ %%
%%%%%%%%%%%%%%%%%%%%%%%%

nStat = numel(iStatField);
nBins = p.nBins;

if p.HistNorm
    histArgs = {'Normalization','probability'};
end

if ~isempty(strfind(p.FigureType,'Histogram'))
    if ~p.HistNorm
        histLab1 = 'No. ';
        histLab3 = '';
    else
        histLab1 = '';
        histLab3 = 'Frequency';
    end
    if strcmp(p.DataType,'Raw')
        histLab2 = 'Cluster ';
    elseif strcmp(p.DataType,'Point-Weighted')
        histLab2 = 'Localization ';
    end
    histYLabel = [histLab1 histLab2 histLab3];
    
end

%% ---- find stat files and determine their conditions ---- %%

statFiles = searchFiles('combined_stats.mat','',allCondDir,true,'new',1);
nStatFiles = numel(statFiles);

if isempty(p.condStr)
    condStr = statFiles;%If the user didn't specify, just name the conditions after the file
else
    condStr = p.condStr;
end

condStr = condStr(:)';%Make sure we have a row cell array

if isempty(p.condNames)
    p.condNames = condStr;
end
condNames = p.condNames;
nCond = numel(condStr);

%Look for each condition matching string in each stat file name            
condMat = cellfun(@(y)(~cellfun(@isempty,cellfun(@(x)(strfind(y,x)),condStr,'Unif',0))),statFiles,'Unif',0); %It's one line but it's one ugly line of code.....
condMat = vertcat(condMat{:});

%Generate condition names for each file
fileCondNames = cell(nStatFiles,1);
for j = 1:nStatFiles    
    fileCondNames{j} = strjoin(condNames(condMat(j,:)),' ');
end

%Load all the stats
stats = cellfun(@load,statFiles);


if isempty(p.iStatFile)
    [iStatFile,OK] = listdlg('Name','File selection','ListString',statFiles,'ListSize',[400 300]);
    if ~OK,return;end
else    
    iStatFile = p.iStatFile;
end

nFiles = numel(iStatFile);%Number of files to compare


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------ DATA PREP ------ %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get N for each condition for reporting
nLocTot = nan(nFiles,1);
nClustTot = nan(nFiles,1);
nROITot = nan(nFiles,1);
condTitleStr = cell(nFiles,1);
for i = 1:nFiles
    nLocTot(i) = sum(stats(iStatFile(i)).combStat.AllClustNumPoints);
    nClustTot(i) = numel(stats(iStatFile(i)).combStat.AllClustNumPoints);
    nROITot(i) = numel(stats(iStatFile(i)).combStat.AllDataSetnClust);
    condTitleStr{i} = [fileCondNames{iStatFile(i)} ': ' num2str(nLocTot(i)) ' localizations, ' num2str(nClustTot(i)) ' clusters, ' num2str(nROITot(i)) ' ROIs'];
end

switch p.DataType        
    
    
    case 'Raw'
        

        % Raw per-cluster data 
        %Biological/statistical unit = cluster
        condData = cell(nFiles,1);

        for i = 1:nFiles
            nPts = numel(stats(iStatFile(i)).combStat.(statFieldNames{1}));
            condData{i} = zeros(nPts,nStat);
            for j = 1:nStat
                condData{i}(:,j) = stats(iStatFile(i)).combStat.(statFieldNames{iStatField(j)});
            end
        end


    case 'Point-Weighted'
        
        %For statistics that are weighted based on the number of localizations
        %Biological/statistical unit = localization

        condData = cell(nFiles,1);

        for i = 1:nFiles

            rawDat = stats(iStatFile(i)).combStat.(statFieldNames{iStatField});
            condData{i} = cell(numel(rawDat),1);
            for j = 1:numel(rawDat)
                condData{i}{j} = repmat(rawDat(j),[stats(iStatFile(i)).combStat.AllClustNumPoints(j) 1]);
            end
            condData{i} = vertcat(condData{i}{:});    

        end
        
    case 'Per-ROI-Mean'
    
        %For cluster statistics that are averaged over each cell
        %Biological/statistical unit = ROI (with sub-unit as cluster)
        
        condData = cell(nFiles,1);        
        for i = 1:nFiles

            rawDat = stats(iStatFile(i)).combStat.(statFieldNames{iStatField});
            
            %Create index relating clusters to ROIs
            roiInd = arrayfun(@(x)(repmat(x,stats(iStatFile(i)).combStat.AllDataSetnClust(x),1)),1:nROITot(i),'Unif',0);
            roiInd = vertcat(roiInd{:});
            roiInd = roiInd(stats(iStatFile(i)).combStat.PointsUsed);%Account for excluded clusters 
            
            %If requested, apply filtering to these clusters
            if ~isempty(p.DataFilterFun)
                keepData = p.DataFilterFun(stats(iStatFile(i)));
                roiInd = roiInd(keepData);
                rawDat = rawDat(keepData);
            end
            
            condData{i} = nan(nROITot(i),1);            
            for j = 1:nROITot(i)
                condData{i}(j) = mean(rawDat(roiInd == j));                
            end
            

        end                
        
    case 'Per-ROI-Mean-Point-Weighted'
    
        %For statistics that are weighted based on the number of
        %localizations and averaged over each cell
        %Biological/statistical unit = ROI (with sub-unit as localization)
        
        condData = cell(nFiles,1);        
        for i = 1:nFiles

            rawDat = stats(iStatFile(i)).combStat.(statFieldNames{iStatField});
            rawNumLoc = stats(iStatFile(i)).combStat.AllClustNumPoints(roiInd == j);
            
            %Create index relating clusters to ROIs
            roiInd = arrayfun(@(x)(repmat(x,stats(iStatFile(i)).combStat.AllDataSetnClust(x),1)),1:nROITot(i),'Unif',0);
            roiInd = vertcat(roiInd{:});
            roiInd = roiInd(stats(iStatFile(i)).combStat.PointsUsed);%Account for excluded clusters 
            
            
            %If requested, apply filtering to these clusters
            if ~isempty(p.DataFilterFun)
                keepData = p.DataFilterFun(stats(iStatFile(i)));
                roiInd = roiInd(keepData);
                rawDat = rawDat(keepData);
                rawNumLoc = rawNumLoc(keepData);
            end
            
            condData{i} = nan(nROITot(i),1);            
            for j = 1:nROITot(i)
                roiRawData = rawDat(roiInd == j);
                roiNumLoc = rawNumLoc(roiInd == j);
                roiDat = cell(numel(roiRawData),1);
                for k = 1:numel(roiRawData)
                    roiDat{k} = repmat(roiRawData(k),[roiNumLoc(k) 1]);                    
                end
                roiDat = vertcat(roiDat{:});
                condData{i}(j) = mean(roiDat);
            end
            

        end                                                        
        
    case 'Per-ROI-NumClusters'
        
        %For cluster counts in each ROI. Can be combined with filtering to
        %count clusters with specific properties
        %Biological/statistical unit = ROI (with sub-unit as cluster)
        
        condData = cell(nFiles,1);        
        for i = 1:nFiles            
            
            %Create index relating clusters to ROIs
            roiInd = arrayfun(@(x)(repmat(x,stats(iStatFile(i)).combStat.AllDataSetnClust(x),1)),1:nROITot(i),'Unif',0);
            roiInd = vertcat(roiInd{:});
            roiInd = roiInd(stats(iStatFile(i)).combStat.PointsUsed);%Account for excluded clusters 
            
            %If requested, apply filtering to these clusters
            if ~isempty(p.DataFilterFun)
                keepData = p.DataFilterFun(stats(iStatFile(i)));
                roiInd = roiInd(keepData);                
            end
            
            condData{i} = nan(nROITot(i),1);            
            for j = 1:nROITot(i)
                condData{i}(j) = nnz(roiInd == j);                
            end
            

        end    
        
    case 'Per-ROI-NumLocalizations'
    
        %For statistics on numbers of localizations averaged over each cell
        %Biological/statistical unit = ROI (with sub-unit as localization)
        
        condData = cell(nFiles,1);        
        for i = 1:nFiles

            %rawDat = stats(iStatFile(i)).combStat.(statFieldNames{iStatField});
            rawNumLoc = stats(iStatFile(i)).combStat.AllClustNumPoints;
            
            %Create index relating clusters to ROIs
            roiInd = arrayfun(@(x)(repmat(x,stats(iStatFile(i)).combStat.AllDataSetnClust(x),1)),1:nROITot(i),'Unif',0);
            roiInd = vertcat(roiInd{:});
            roiInd = roiInd(stats(iStatFile(i)).combStat.PointsUsed);%Account for excluded clusters 
            
            
            %If requested, apply filtering to these clusters
            if ~isempty(p.DataFilterFun)
                keepData = p.DataFilterFun(stats(iStatFile(i)));
                roiInd = roiInd(keepData);
                %rawDat = rawDat(keepData);
                rawNumLoc = rawNumLoc(keepData);
            end
            
            condData{i} = nan(nROITot(i),1);            
            for j = 1:nROITot(i)                                
                condData{i}(j) = sum(rawNumLoc(roiInd == j));                
            end
            

        end            

        
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------ STATISTICS ------ %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch p.StatToCalc
    
    case 'Median'
        
        % ------- Median of data ------ %


        isNorm = false(nFiles,1);
        condDataGroup = cell(nFiles,1);


        condScalStat = nan(nFiles,1); %Scalar statistic for each condition
        condScalErr = nan(nFiles,1); %Scalar uncertainty for each condition
        statName = 'Median';
        statFun = @median;
        statFunErr = @(x)(std(x) / sqrt(numel(x)));


        for i = 1:nFiles

            condDataGroup{i} = i * ones(numel(condData{i}),1);
            condScalStat(i) = statFun(condData{i});
            condScalErr(i) = statFunErr(condData{i});
            isNorm(i) = ~adtest(condData{i});

        end

        isNorm

        [condP,condStatTab,condCompStats]= kruskalwallis(vertcat(condData{:}),vertcat(condDataGroup{:}))

        [condCompMat,condCompStat] = multcompare(condCompStats)
        
        statTable = table(fileCondNames(iStatFile(condCompMat(:,1))),fileCondNames(iStatFile(condCompMat(:,2))),condCompMat(:,6),...
                    'VariableNames',{'CondtionA','ConditionB','PValue'});
        
        
    case ''
        

    otherwise
            
            error('Unrecognized data type!')
            
end

            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------ Figures / Stats ------ %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch p.FigureType

    
    case 'Bar'
        % ---------- Bar graph condition comparison -------- %

        cf = figure;
        hold on
        plotCols = cf.CurrentAxes.ColorOrder;
        for j = 1:nFiles

            bHan = bar(j,condScalStat(j),'FaceColor',plotCols(j,:),p.PlotArgs{:});
            errorbar(j,condScalStat(j),condScalErr(j),'k')

        end

        cf.CurrentAxes.XTick = [];           
        cf.CurrentAxes.YLabel.String  = statNames{iStatField};

        ylabel([statName ' ' statNames{iStatField}])
        
    case 'Box'
        % ---------- Box plot condition comparison -------- %

        tmpDat = vertcat(condData{:});
        plotLim = prctile(tmpDat,p.PctLimPlot);

        cf = figure;
        hold on
        plotCols = cf.CurrentAxes.ColorOrder;

        boxplot(vertcat(condData{:}),vertcat(condDataGroup{:}),'colors',plotCols(1:nFiles,:),'symbol','k.',p.PlotArgs{:})        
        
        allObj = cf.CurrentAxes.Children;
        set(allObj.Children,'LineWidth',2);

        cf.CurrentAxes.YLim = [0 plotLim(2)];
        cf.CurrentAxes.XTick = [];
        cf.CurrentAxes.YLabel.String  = statNames{iStatField};

    case 'BarHistogram'


        % ------------ Bar Histogram Condition Comparision -------- %

        tmpDat = vertcat(condData{:});
        binLim = prctile(tmpDat,p.PctLimHist);
        binEdges = linspace(binLim(1),binLim(2),nBins+1);

        cf = figure;
        hold on
        clear histHan
        for j = 1:nFiles
            
            histHan(j) = histogram(condData{j},binEdges,p.PlotArgs{:},histArgs{:});
            histHan(j).FaceAlpha = .6;
        end

        xlabel(statNames{iStatField})
        ylabel(histYLabel)
        legend(fileCondNames(iStatFile))
        title([p.DataType ' histograms' condTitleStr'])

    case 'FilledHistogram'
        
        % ------------ Filled Line Histogram Condition Comparision -------- %        

        tmpDat = vertcat(condData{:});
        binLim = prctile(tmpDat,p.PctLimHist);
        binEdges = linspace(binLim(1),binLim(2),nBins+1);

        cf = figure;
        hold on

        plotCols = cf.CurrentAxes.ColorOrder;

        histVals = nan(nFiles,nBins);
        clear histHan
        for j = 1:nFiles
            histVals(j,:) = histcounts(condData{j},binEdges,histArgs{:});    


            switch p.FilledHistType

                case 'Line'            
                    xVals = [binEdges(1:end-1) binEdges(end-1) binEdges(1)];
                    yVals = [histVals(j,:) 0 0];
                case 'Stairs'
                    [xVals,yVals] = stairs(binEdges(1:end-1),histVals(j,:));
                    xVals = [xVals' xVals(end) xVals(1)];
                    yVals = [yVals' 0 0];
            end

            histHan(j) = patch(xVals,yVals,plotCols(j,:),'FaceAlpha',.6,p.PlotArgs{:}); %Black edges, colored faces
            %histHan(j) = patch(xVals,yVals,plotCols(j,:),'EdgeColor',plotCols(j,:),'FaceColor','none','LineWidth',2); %Colored edges no faces


        end

        xlim([binEdges(1) binEdges(end-1)])
        xlabel(statNames{iStatField})
        ylabel(histYLabel)
        legend(fileCondNames(iStatFile))
        title([p.DataType ' histograms' condTitleStr'])
        
    case 'Scatter'


        % -------- 3D scatter plot Condition Comparison --------- %

        cf = figure;
        hold on

        tmpDat = vertcat(condData{:});
        axLim = zeros(nStat,2);
        for j = 1:nStat
            axLim(j,:) = prctile(tmpDat(:,j),p.PctLimPlot);
        end   

        for j = 1:nFiles
            plot3(condData{j}(:,1),condData{j}(:,2),condData{j}(:,3),'.',p.PlotArgs{:})
        end
        legend(fileCondNames(iStatFile))

        view([132.5 32])
        xlim(axLim(1,:))
        xlabel(statNames{iStatField(1)})
        ylim(axLim(2,:))
        ylabel(statNames{iStatField(2)})
        zlim(axLim(3,:))
        zlabel(statNames{iStatField(3)})

     
    otherwise
        error('Unrecognized figure type!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ------ Figures / Stats ------ %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(p.OutputFile)
    [outFile,outPath] = uiputfile('','Save figure to:');
    if outFile == 0
        return
    else
        p.OutputFile = [outPath outFile];
    end
end
    
mfFigureExport(cf,p.OutputFile);

if exist('statTable','var')
    writetable(statTable,[p.OutputFile '.csv']);
end


