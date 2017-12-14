
%allCondDir =
%'/home/he19/files/CellBiology/IDAC/Hunter/Soberman/Angie/Angie_Figure_Data
%for Hunter_6_21_2016/'; %Directory stat files for first two figures
allCondDir = '/home/he19/files/CellBiology/IDAC/Hunter/Soberman/Angie/Angie_Combined_Paper_Data'; %Directory containing all the cluster stat files
condStr   = {'5LO','FLAP','Control','GM','C5a','0''','5''','2 color','cPLA2 Inh','5LO knockout'};  %Strings to search for in the file names
statFileNames = {'5LO','FLAP','Control','+ GM-CSF','+ C5a','0''','5''','2 color','+ cPLA2 Inh','5LO knockout'}; %Strings for figure / output naming


outDir = '/home/he19/files/CellBiology/IDAC/Hunter/Soberman/Angie/Figures/Figs_10_25_2016';


mkdir(outDir);

%%%%%%%%%%%%%%%%%%%%%%%%
%% ===== 5 - LO  ==== %%

%% ---- Figure 3c-e ----- %%
%5-LO activation
%Histogram overlays and bar graphs


iStatFile = [19 17 18];
outFileNameBase = [outDir filesep 'Figure 3'];
iStatFields = 1:3;
statFieldNames = {'number','area','density'};
panelLetters = {'C','D','E'};

for iStat = 1:numel(iStatFields)


    %Make histogram overlays
    outFileName = [outFileNameBase panelLetters{iStat} ' 5-LO activation histogram ' statFieldNames{iStat}];
        
    stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',statFileNames,'OutputFile',outFileName,...
                    'iStatFile',iStatFile,'iStatField',iStatFields(iStat),'DataType','Point-Weighted','FigureType','FilledHistogram')

    %Make per-cell bar graphs for inset
    outFileName = [outFileNameBase panelLetters{iStat} ' 5-LO activation per cell point weighted boxplot ' statFieldNames{iStat}];                
                    
    
    stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',statFileNames,'OutputFile',outFileName,...
                    'iStatFile',iStatFile,'iStatField',iStatFields(iStat),'DataType','Per-ROI-Mean','FigureType','Box','StatToCalc','Median',...
                    'PctLimPlot',[0 95],'PlotArgs',{'notch','on'});

end



%% ---- Figure 3f ----- %%
%5-LO activation & cpla2 inhibitor
%Histogram overlay and bar graphs

iStatFile = [7 5 6];
outFileNameBase = [outDir filesep 'Figure 3'];
iStatFields = 1;
statFieldNames = {'number'};
panelLetters = {'F'};
iStat = 1;

%Make histogram overlays
outFileName = [outFileNameBase panelLetters{iStat} ' 5-LO cpla2inh histogram ' statFieldNames{iStat}];

stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',statFileNames,'OutputFile',outFileName,...
                'iStatFile',iStatFile,'iStatField',iStatFields(iStat),'DataType','Point-Weighted','FigureType','FilledHistogram')

%Make bar graphs for inset
outFileName = [outFileNameBase panelLetters{iStat} ' 5-LO cpla2inh boxplot ' statFieldNames{iStat}];                


stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',statFileNames,'OutputFile',outFileName,...
                'iStatFile',iStatFile,'iStatField',iStatFields(iStat),'DataType','Point-Weighted','FigureType','Box','StatToCalc','Median',...
                'PctLimPlot',[0 95]);

%% ---- Figure 3?? ----- %%
%Pairwise comparison of 5-LO activation & cpla2 inhibitor
%Histogram overlay and bar graphs




outFileNameBase = [outDir filesep '5-LO cpla2inh vs '];
iStatFields = 1:3;
statFieldNames = {'number','area','density'};


statFilesInd = {[18 6],...
                 [17 5],...
                 [19 7]};
statFileIndNames = {'control','C5A','C5a+CSF'};
            
             

for iStat = 1:numel(iStatFields)
    
    for iCond = 1:numel(statFileIndNames)
    
        iStatFile = statFilesInd{iCond};        

        %Make histogram overlays
        outFileName = [outFileNameBase statFileIndNames{iCond} ' histogram pairwise ' statFieldNames{iStat}];

        stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',statFileNames,'OutputFile',outFileName,...
                        'iStatFile',iStatFile,'iStatField',iStatFields(iStat),'DataType','Point-Weighted','FigureType','FilledHistogram')

        %Make bar graphs for inset
        outFileName = [outFileNameBase statFileIndNames{iCond} ' boxplot pairwise ' statFieldNames{iStat}];                


        stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',statFileNames,'OutputFile',outFileName,...
                        'iStatFile',iStatFile,'iStatField',iStatFields(iStat),'DataType','Point-Weighted','FigureType','Box','StatToCalc','Median',...
                        'PctLimPlot',[0 95]);            

    end


end
%% ---- 5LO Activation #Big Clusters ----- %%

minClustLoc = 437; %This is set equal to the 95th percentile of localization-weighted cluster sizes in the control 5-LO data.
%minClustLoc = 140; %This is set equal to the 95th percentile of cluster sizes in the control 5-LO data.
clustFiltFun = @(x)(x.combStat.AllClustNumPoints >= minClustLoc);


iStatFile = [19 17 18];
outFileNameBase = [outDir filesep 'Figure X'];
iStatFields = 1;



%     %Make histogram overlays
%     outFileName = [outFileNameBase panelLetters{iStat} ' 5-LO activation clusters over ' num2str(minClustLoc) ' localizations histogram ' statFieldNames{iStat}];
%         
%     stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',condNames,'OutputFile',outFileName,...
%                     'iStatFile',iStatFile,'iStatField',iStatFields(iStat),'DataType','Point-Weighted','FigureType','FilledHistogram',...
%                     'DataFilterFun',clustFiltFun)

    %Make per-cell bar graphs for inset
outFileName = [outFileNameBase ' 5-LO activation number of clusters over ' num2str(minClustLoc) ' localizations per cell boxplot'];                
                    
    
stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',statFileNames,'OutputFile',outFileName,...
                'iStatFile',iStatFile,'iStatField',iStatFields,'DataType','Per-ROI-NumClusters','FigureType','Box','StatToCalc','Median',...
                'PctLimPlot',[0 95],'PlotArgs',{'notch','on'},'DataFilterFun',clustFiltFun);

%% ---- 5LO Activation #Localizations in Big Clusters ----- %%

%minClustLoc = 437; %This is set equal to the 95th percentile of localization-weighted cluster sizes in the control 5-LO data.
minClustLoc = 140; %This is set equal to the 95th percentile of cluster sizes in the control 5-LO data.
clustFiltFun = @(x)(x.combStat.AllClustNumPoints >= minClustLoc);


iStatFile = [19 17 18];
outFileNameBase = [outDir filesep 'Figure X'];
iStatFields = 1;



%     %Make histogram overlays
%     outFileName = [outFileNameBase panelLetters{iStat} ' 5-LO activation clusters over ' num2str(minClustLoc) ' localizations histogram ' statFieldNames{iStat}];
%         
%     stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',condNames,'OutputFile',outFileName,...
%                     'iStatFile',iStatFile,'iStatField',iStatFields(iStat),'DataType','Point-Weighted','FigureType','FilledHistogram',...
%                     'DataFilterFun',clustFiltFun)

    %Make per-cell bar graphs for inset
outFileName = [outFileNameBase ' 5-LO activation number of localizations in clusters over ' num2str(minClustLoc) ' localizations per cell boxplot'];                
                    
    
stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',statFileNames,'OutputFile',outFileName,...
                'iStatFile',iStatFile,'iStatField',iStatFields,'DataType','Per-ROI-NumLocalizations','FigureType','Box','StatToCalc','Median',...
                'PctLimPlot',[0 95],'PlotArgs',{'notch','on'},'DataFilterFun',clustFiltFun);

            
            
%%%%%%%%%%%%%%%%%%%%%%%%
%% ====  FLAP    ==== %%            
            
%% ---- Figure 4b-d ----- %%
%FLAP activation
%Histogram overlays and bar graphs


iStatFile = [22 20 21];
outFileNameBase = [outDir filesep 'Figure 4'];
iStatFields = 1:3;
statFieldNames = {'number','area','density'};
panelLetters = {'B','C','D'};

for iStat = 1:numel(iStatFields)


    %Make histogram overlays
    outFileName = [outFileNameBase panelLetters{iStat} ' FLAP activation histogram ' statFieldNames{iStat}];
        
    stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',statFileNames,'OutputFile',outFileName,...
                    'iStatFile',iStatFile,'iStatField',iStatFields(iStat),'DataType','Point-Weighted','FigureType','FilledHistogram')

    %Make bar graphs for inset
    outFileName = [outFileNameBase panelLetters{iStat} ' FLAP activation boxplot ' statFieldNames{iStat}];                
                    
    
    stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',statFileNames,'OutputFile',outFileName,...
                    'iStatFile',iStatFile,'iStatField',iStatFields(iStat),'DataType','Point-Weighted','FigureType','Box','StatToCalc','Median',...
                    'PctLimPlot',[0 95]);

end

%% ---- FLAP and cpla2Inhibitor pairwise ----- %%
%Pairwise comparison of 5-LO activation & cpla2 inhibitor
%Histogram overlay and bar graphs

outFileNameBase = [outDir filesep 'FLAP cpla2inh vs '];
iStatFields = 1:3;
statFieldNames = {'number','area','density'};


statFilesInd = {[21 9],...
                 [20 8],...
                 [22 10]};
statFileIndNames = {'control','C5A','C5a+CSF'};
            
             

for iStat = 1:numel(iStatFields)
    
    for iCond = 1:numel(statFileIndNames)
    
        iStatFile = statFilesInd{iCond};        

        %Make histogram overlays
        outFileName = [outFileNameBase statFileIndNames{iCond} ' histogram pairwise ' statFieldNames{iStat}];

        stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',statFileNames,'OutputFile',outFileName,...
                        'iStatFile',iStatFile,'iStatField',iStatFields(iStat),'DataType','Point-Weighted','FigureType','FilledHistogram')

        %Make bar graphs for inset
        outFileName = [outFileNameBase statFileIndNames{iCond} ' boxplot pairwise ' statFieldNames{iStat}];                


        stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',statFileNames,'OutputFile',outFileName,...
                        'iStatFile',iStatFile,'iStatField',iStatFields(iStat),'DataType','Point-Weighted','FigureType','Box','StatToCalc','Median',...
                        'PctLimPlot',[0 95]);            

    end


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ===== Knockout Data ===== %%

%% ---- Figure Y ----- %%
%FLAP activation in 5-LO KO cells
%Histogram overlays and bar graphs


iStatFile = [16 14 15];
outFileNameBase = [outDir filesep 'Figure Y 5-LO knockout '];
iStatFields = 1:3;
statFieldNames = {'number','area','density'};
panelLetters = {'B','C','D'};

for iStat = 1:numel(iStatFields)


    %Make histogram overlays
    outFileName = [outFileNameBase panelLetters{iStat} ' FLAP activation histogram ' statFieldNames{iStat}];
        
    stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',statFileNames,'OutputFile',outFileName,...
                    'iStatFile',iStatFile,'iStatField',iStatFields(iStat),'DataType','Point-Weighted','FigureType','FilledHistogram')

    %Make bar graphs for inset
    outFileName = [outFileNameBase panelLetters{iStat} ' FLAP activation boxplot ' statFieldNames{iStat}];                
                    
    
    stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',statFileNames,'OutputFile',outFileName,...
                    'iStatFile',iStatFile,'iStatField',iStatFields(iStat),'DataType','Point-Weighted','FigureType','Box','StatToCalc','Median',...
                    'PctLimPlot',[0 95]);

end

%% ---- Figure Z ----- %%
%FLAP activation in 5-LO KO cells + cPLA2 Inh
%Histogram overlays and bar graphs


iStatFile = [13 11 12];
outFileNameBase = [outDir filesep 'Figure Z 5-LO knockout and cPLA2Inh '];
iStatFields = 1:3;
statFieldNames = {'number','area','density'};
panelLetters = {'B','C','D'};

for iStat = 1:numel(iStatFields)


    %Make histogram overlays
    outFileName = [outFileNameBase panelLetters{iStat} ' FLAP activation histogram ' statFieldNames{iStat}];
        
    stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',statFileNames,'OutputFile',outFileName,...
                    'iStatFile',iStatFile,'iStatField',iStatFields(iStat),'DataType','Point-Weighted','FigureType','FilledHistogram')

    %Make bar graphs for inset
    outFileName = [outFileNameBase panelLetters{iStat} ' FLAP activation boxplot ' statFieldNames{iStat}];                
                    
    
    stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',statFileNames,'OutputFile',outFileName,...
                    'iStatFile',iStatFile,'iStatField',iStatFields(iStat),'DataType','Point-Weighted','FigureType','Box','StatToCalc','Median',...
                    'PctLimPlot',[0 95]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ===== Timecourse Data ===== %%

%% ---- Figure 5 ----- %%
%FLAP activation timecourse
%Histogram overlays and bar graphs


iStatFile = [2 4];
outFileNameBase = [outDir filesep 'Figure 5'];
iStatFields = 1:3;
statFieldNames = {'number','area','density'};
panelLetters = {'E1','E2','E3'};

for iStat = 1:numel(iStatFields)


    %Make histogram overlays
    outFileName = [outFileNameBase panelLetters{iStat} ' 2-Color FLAP timecourse histogram ' statFieldNames{iStat}];
        
    stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',statFileNames,'OutputFile',outFileName,...
                    'iStatFile',iStatFile,'iStatField',iStatFields(iStat),'DataType','Point-Weighted','FigureType','FilledHistogram')

    %Make bar graphs for inset
    outFileName = [outFileNameBase panelLetters{iStat} ' 2-Color FLAP timecourse boxplot ' statFieldNames{iStat}];                
                    
    
    stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',statFileNames,'OutputFile',outFileName,...
                    'iStatFile',iStatFile,'iStatField',iStatFields(iStat),'DataType','Point-Weighted','FigureType','Box','StatToCalc','Median',...
                    'PctLimPlot',[0 95]);

end

%% ---- Figure 5 ----- %%
%5-LO activation timecourse
%Histogram overlays and bar graphs


iStatFile = [1 3];
outFileNameBase = [outDir filesep 'Figure 5'];
iStatFields = 1:3;
statFieldNames = {'number','area','density'};
panelLetters = {'E4','E5','E6'};

for iStat = 1:numel(iStatFields)


    %Make histogram overlays
    outFileName = [outFileNameBase panelLetters{iStat} ' 2-Color 5-LO timecourse histogram ' statFieldNames{iStat}];
        
    stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',statFileNames,'OutputFile',outFileName,...
                    'iStatFile',iStatFile,'iStatField',iStatFields(iStat),'DataType','Point-Weighted','FigureType','FilledHistogram')

    %Make bar graphs for inset
    outFileName = [outFileNameBase panelLetters{iStat} ' 2-Color 5-LO timecourse boxplot ' statFieldNames{iStat}];                
                    
    
    stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',statFileNames,'OutputFile',outFileName,...
                    'iStatFile',iStatFile,'iStatField',iStatFields(iStat),'DataType','Point-Weighted','FigureType','Box','StatToCalc','Median',...
                    'PctLimPlot',[0 95]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ===== 2-Color Cluster overlay Visualizations===== %%


%Columns are channels rows are cells
%CHANGE THESE TO THE PATHS / CONDITIONS YOU WANT
clustFiles = {'full/path/to/file for control 5-LO cluster stats.mat','full/path/to/file for control FLAP cluster stats.mat'; ...
              'full/path/to/file for gmcsf 5-LO cluster stats.mat','full/path/to/file for gmcsf FLAP cluster stats.mat'};
              

cellNames = {'Control','gmcsf_c5a'};%YOU CAN CHANGE THIS TO OTHER CONDITIONS IF YOU WANT TO SHOW SOMETHING ELSE
chanNames = {'5LO','FLAP'};

nCells = size(clustFiles,1);
nChan = size(clustFiles,2);

for j = 1:nCells
    for k = 1:nChan
        cDat(j,k) = load(clustFiles{j,k});
        allDensData{j,k} = cDat(j,k).clustStats.Density;
    end
   
end

%CHANGE THIS TO THE DIRECTORY YOU WANT TO SAVE THE RESULTS TO
outDir = '/home/he19/files/CellBiology/IDAC/Hunter/Soberman/Angie/Figures/TwoColor_Overlays_11_8_2016';
outDPI = 600;

%Get limits for density display
densLims = nan(2,nChan);
for j = 1:nChan
    tmp = prctile(vertcat(allDensData{:,j}),[10 90]);
    densLims(:,j) = tmp;
end
%densLims = prctile(vertcat(allDensData{:}),[0 95]);

minAlpha = .1;%Minimuim opacity for opacity mapping
minSat = .1;%Minimum saturation for saturation mapping

chanCols = [0 1 0 ; 1 0 0];

saveFigs = true;

%% --- 2-Color overlay of points ---- %%

iCell = 2;



cf = figure;
hold on
for iChan= 1:nChan    
    plot(cDat(iCell,iChan).stormData.X,cDat(iCell,iChan).stormData.Y,'.','color',chanCols(iChan,:));
    disp(clustFiles{iCell,iChan})
end
    
axis equal

if saveFigs
    mfFigureExport(cf,[outDir filesep 'Per cluster points overlay ' cellNames{iCell} ' ' chanNames{iChan} ' only'],'DPI',outDPI)
end
    
%% ---  Convex hull overlay ---- %%




for iCell = 1:nCells

    
    % --- make and save the overlay --- %

    cf = figure;
    hold on
    for iChan = 1:nChan

        nClust = cDat(iCell,iChan).clustStats.nClust;
        for iClust = 1:nClust
            if ~isempty(cDat(iCell,iChan).clustStats.Boundary{iClust})
                patch(cDat(iCell,iChan).clustStats.Boundary{iClust}(:,1),cDat(iCell,iChan).clustStats.Boundary{iClust}(:,2),chanCols(iChan,:),'FaceAlpha',.33)
            end
        end

    end
    axis equal
    axis off
    if saveFigs
        mfFigureExport(cf,[outDir filesep 'Convex hulls overlay ' cellNames{iCell}],'DPI',outDPI)
    end

    
    % --- make and save them separately as well so we can do blending in illustrator--- %
    

    
    for iChan = 1:nChan
        
        cf = figure;
        hold on

        nClust = cDat(iCell,iChan).clustStats.nClust;
        for iClust = 1:nClust
            if ~isempty(cDat(iCell,iChan).clustStats.Boundary{iClust})
                patch(cDat(iCell,iChan).clustStats.Boundary{iClust}(:,1),cDat(iCell,iChan).clustStats.Boundary{iClust}(:,2),chanCols(iChan,:),'FaceAlpha',.33)
            end
        end
        
            axis equal
            axis off
            if saveFigs
                mfFigureExport(cf,[outDir filesep 'Convex hulls overlay ' cellNames{iCell} ' ' chanNames{iChan} ' only'],'DPI',outDPI)
            end

    end


    


end

%% ---  Convex hull overlay with density as opacity ---- %%




for iCell = 1:nCells

    
    % --- make and save the overlay --- %

    cf = figure;
    hold on
    for iChan = 1:nChan

        nClust = cDat(iCell,iChan).clustStats.nClust;
        for iClust = 1:nClust
            if ~isempty(cDat(iCell,iChan).clustStats.Boundary{iClust})
                
                densNorm = min(max(0,(cDat(iCell,iChan).clustStats.Density(iClust) - densLims(1,iChan)) / densLims(2,iChan)),1);
                densNorm = (densNorm+minAlpha) / (1+minAlpha);%Make sure the most transparent are still somehwat opaque
                
                patch(cDat(iCell,iChan).clustStats.Boundary{iClust}(:,1),cDat(iCell,iChan).clustStats.Boundary{iClust}(:,2),chanCols(iChan,:),'FaceAlpha',densNorm)
            end
        end

    end
    axis equal
    axis off
    if saveFigs
        mfFigureExport(cf,[outDir filesep 'Convex hulls overlay opacity mapped ' cellNames{iCell}],'DPI',outDPI)
    end

    
    % --- make and save them separately as well so we can do blending in illustrator--- %
    

    
    for iChan = 1:nChan
        
        cf = figure;
        hold on

        nClust = cDat(iCell,iChan).clustStats.nClust;
        for iClust = 1:nClust
            if ~isempty(cDat(iCell,iChan).clustStats.Boundary{iClust})
                densNorm = min(max(0,(cDat(iCell,iChan).clustStats.Density(iClust) - densLims(1,iChan)) / densLims(2,iChan)),1);
                densNorm = (densNorm+minAlpha) / (1+minAlpha);%Make sure the most transparent are still somehwat opaque
                
                patch(cDat(iCell,iChan).clustStats.Boundary{iClust}(:,1),cDat(iCell,iChan).clustStats.Boundary{iClust}(:,2),chanCols(iChan,:),'FaceAlpha',densNorm)
            end
        end
        
            axis equal
            axis off
            if saveFigs
                mfFigureExport(cf,[outDir filesep 'Convex hulls overlay opacity mapped ' cellNames{iCell} ' ' chanNames{iChan} ' only'],'DPI',outDPI)
            end

    end


    


end

%% ---  Convex hull overlay with density as saturation ---- %%

faceAlpha = .95;


for iCell = 1:nCells

    
    % --- make and save the overlay --- %

    cf = figure;
    hold on
    for iChan = 1:nChan

        nClust = cDat(iCell,iChan).clustStats.nClust;
        for iClust = 1:nClust
            if ~isempty(cDat(iCell,iChan).clustStats.Boundary{iClust})
                
                densNorm = min(max(0,(cDat(iCell,iChan).clustStats.Density(iClust) - densLims(1,iChan)) / densLims(2,iChan)),1);
                densNorm = (densNorm+minSat) / (1+minSat);%Make sure the least saturated still have some color
                
                %Set the saturation to the relative density
                mappedCol = rgb2hsv(chanCols(iChan,:));
                mappedCol(2) = densNorm;
                mappedCol = hsv2rgb(mappedCol);
                
                patch(cDat(iCell,iChan).clustStats.Boundary{iClust}(:,1),cDat(iCell,iChan).clustStats.Boundary{iClust}(:,2),mappedCol,'FaceAlpha',faceAlpha)
            end
        end

    end
    axis equal
    axis off
    if saveFigs
        mfFigureExport(cf,[outDir filesep 'Convex hulls overlay saturation mapped ' cellNames{iCell}],'DPI',outDPI)
    end

    
    % --- make and save them separately as well so we can do blending in illustrator--- %
    

    
    for iChan = 1:nChan
        
        cf = figure;
        hold on

        nClust = cDat(iCell,iChan).clustStats.nClust;
        for iClust = 1:nClust
            if ~isempty(cDat(iCell,iChan).clustStats.Boundary{iClust})
                densNorm = min(max(0,(cDat(iCell,iChan).clustStats.Density(iClust) - densLims(1,iChan)) / densLims(2,iChan)),1);
                densNorm = (densNorm+minSat) / (1+minSat);%Make sure the least saturated still have some color
                
                %Set the saturation to the relative density
                mappedCol = rgb2hsv(chanCols(iChan,:));
                mappedCol(2) = densNorm;
                mappedCol = hsv2rgb(mappedCol);
                
                patch(cDat(iCell,iChan).clustStats.Boundary{iClust}(:,1),cDat(iCell,iChan).clustStats.Boundary{iClust}(:,2),mappedCol,'FaceAlpha',faceAlpha)
            end
        end
        
            axis equal
            axis off
            if saveFigs
                mfFigureExport(cf,[outDir filesep 'Convex hulls overlay saturation mapped ' cellNames{iCell} ' ' chanNames{iChan} ' only'],'DPI',outDPI)
            end

    end


    


end
