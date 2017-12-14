%% --- Parameters ---- %%

%outDir = '/home/he19/files/CellBiology/IDAC/Hunter/Soberman/Figures/Figure Panels';%Parent output directory

outFigDir = '/home/he19/files/CellBiology/IDAC/Hunter/Soberman/Angie/Figures/Figs_6_24_2016'

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% ----- 5LO TimeCourse Figures ------- %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %Load the condition info
% tmp = load('/home/he19/files/CellBiology/IDAC/Hunter/Soberman/Data_3_26_2014/condition_info_file_paths_11_6_2014.mat');%This info is good for the timecourse (some other conditions need fixing in it though)
% condFiles = tmp.condStatFiles;
% 
% %-----5LO Cells -----%
% 
% condFilesCells = condFiles([35  30 29 28])'; %5LO Cells
% tData = [0 2 5 10];
% currOut = [outDir filesep '5LO Cells Time Course'];
% if ~exist(currOut,'dir')
%     mkdir(currOut);
% end
% stormTimeCourseAnalysis(condFilesCells,tData,'TimeUnits','min','OutputDirectory',currOut)
% 
% %-----5LO Membranes-----%
% 
% condFilesMem = condFiles([23 18 17 16 ])';
% tData = [0 2 5 10];
% currOut = [outDir filesep '5LO Membrane Time Course'];
% if ~exist(currOut,'dir')
%     mkdir(currOut);
% end
% stormTimeCourseAnalysis(condFilesMem,tData,'TimeUnits','min','OutputDirectory',currOut)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----- Clustering example "Proof of Principle" Figures ----- %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exampleClustFile = '/home/he19/files/CellBiology/IDAC/Hunter/Soberman/Angie/Angie_All_Data for Hunter_6_21_2016/ROI Data For Proof of Principle Fig/0'' 5LO ROI/roi cluster stats.mat';

[exampleClusterDir,~,~] = fileparts(exampleClustFile);

load(exampleClustFile)


saveFigs = true;


cropLimsXX_YY = 1e4 * [1.3700    1.6613 ; 2.0212    2.3177]; %For 0' 5LO
%cropLimsImXX_YY = 1e3 * [ .3958 1.0343; .4480 1.0865];

%% --

[N,C] = hist3([stormData.X stormData.Y],[1e3 1e3]);

N = filterGauss3D(N,1.2);

cf = figure;
imshow(N',[],'XData',[min(stormData.X) max(stormData.X)],'YData',[min(stormData.Y) max(stormData.Y)])
saturateImageColormap(gca,1)
colormap hot

h = plotScaleBar(1e3,'Location','SouthWest');



if saveFigs
    mfFigureExport(cf,[outFigDir filesep 'example STORM pseudoimage multiple cells'],'DPI',300)
end

%Show ROI
patch(cropLimsXX_YY(1,[1 2 2 1]),cropLimsXX_YY(2,[1 1 2 2]),'w','FaceColor','none','EdgeColor','w','LineWidth',5)

if saveFigs
    mfFigureExport(cf,[outFigDir filesep 'example STORM pseudoimage multiple cells ROI outline'],'DPI',300)
end


% show a crop
xlim(cropLimsXX_YY(1,:)),ylim(cropLimsXX_YY(2,:))
saturateImageColormap(gca,1)
h = plotScaleBar(1e3,1e2,'Location','SouthWest');

if saveFigs
    mfFigureExport(cf,[outFigDir filesep 'example STORM pseudoimage multiple cells ROI'],'DPI',600)
end


%% --- example clustering by hulls

cf = open([exampleClusterDir filesep 'cluster convex hulls plot.fig'])


axis ij
axis off
title('')
xlabel('')
ylabel('')

h = plotScaleBar(1e3);
h.FaceColor = 0 * ones(1,3);

if saveFigs
    mfFigureExport(cf,[outFigDir filesep 'example STORM cluster hulls multiple cells'],'DPI',600)
end

xlim(cropLimsXX_YY(1,:)),ylim(cropLimsXX_YY(2,:))

h = plotScaleBar(1e3,1e2,'Location','SouthWest','Color',[0 0 0]);

if saveFigs
    mfFigureExport(cf,[outFigDir filesep 'example STORM cluster hulls multiple cells ROI'],'DPI',600)
end

%% --- example clustering by color

cf = open([exampleClusterDir filesep 'cluster assignment plot.fig'])


axis ij
axis off
title('')
xlabel('')
ylabel('')

pts = cf.CurrentAxes.Children;
pts(2).Marker = '.';
pts(1).Marker = '.';
pts(1).MarkerEdgeColor = 'k';

cMap = randomColormap(numel(unique(pointInd)),44);
colormap(cMap)


h = plotScaleBar(1e3);
h.FaceColor = 0 * ones(1,3);


if saveFigs
    mfFigureExport(cf,[outFigDir filesep 'example STORM cluster assignment multiple cells'],'DPI',600)
end

xlim(cropLimsXX_YY(1,:)),ylim(cropLimsXX_YY(2,:))
h = plotScaleBar(1e3,1e2,'Location','SouthWest','Color',[0 0 0]);

if saveFigs
    mfFigureExport(cf,[outFigDir filesep 'example STORM cluster assignment multiple cells ROI'],'DPI',600)
end


%% --- example clustering by both

cf = open([exampleClusterDir filesep 'cluster convex hulls plot.fig'])


axis ij
axis off
title('')
xlabel('')
ylabel('')

plot(stormData.X,stormData.Y,'.','color',.5 * ones(1,3))


h = plotScaleBar(1e3);
h.FaceColor = 0 * ones(1,3);

if saveFigs
    mfFigureExport(cf,[outFigDir filesep 'example STORM cluster hulls multiple cells'],'DPI',600)
end

xlim(cropLimsXX_YY(1,:)),ylim(cropLimsXX_YY(2,:))
h = plotScaleBar(1e3,1e2,'Location','SouthWest','Color',[0 0 0]);

if saveFigs
    mfFigureExport(cf,[outFigDir filesep 'example STORM cluster hulls multiple cells ROI'],'DPI',600)
end

%% --- example points only

%cf = open([exampleClusterDir filesep 'cluster convex hulls plot.fig'])
cf = figure;

plot(stormData.X,stormData.Y,'.','color',.5 * ones(1,3))
axis equal
axis ij
axis off
title('')
xlabel('')
ylabel('')


h = plotScaleBar(1e3);
h.FaceColor = 0 * ones(1,3);


if saveFigs
    mfFigureExport(cf,[outFigDir filesep 'example STORM localizations multiple cells'],'DPI',600)
end

xlim(cropLimsXX_YY(1,:)),ylim(cropLimsXX_YY(2,:))
h = plotScaleBar(1e3,1e2,'Location','SouthWest','Color',[0 0 0]);

if saveFigs
    mfFigureExport(cf,[outFigDir filesep 'example STORM localizations multiple cells ROI'],'DPI',600)
end
