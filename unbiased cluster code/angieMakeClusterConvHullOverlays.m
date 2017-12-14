
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
