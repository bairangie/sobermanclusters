

exampleClustFile = '/home/he19/files/CellBiology/IDAC/Hunter/Soberman/New_Data_2_11_2014/121023-5lo-034_list-2012-12-07-14-30-37_S30_clustering_wholeImageStats/roi cluster stats.mat';
exampleClusterDir = '/home/he19/files/CellBiology/IDAC/Hunter/Soberman/New_Data_2_11_2014/121023-5lo-034_list-2012-12-07-14-30-37_S30_clustering_wholeImageStats';

load(exampleClustFile)


saveFigs = true;
outFigDir = '/home/he19/files/CellBiology/IDAC/Hunter/P41/VBMS_Prelim/Figures';

cropLimsXX_YY = 1e4 * [.7997 1.2683 ; .38468 .87215];
cropLimsImXX_YY = 1e3 * [ .3958 1.0343; .4480 1.0865];

%% --

[N,C] = hist3([stormData.Xc stormData.Yc],[5e3 5e3]);

N = filterGauss3D(N,1.2);

cf = figure;
imshow(N',[],'XData',[min(stormData.Xc) max(stormData.Xc)],'YData',[min(stormData.Yc) max(stormData.Yc)])
saturateImageColormap(gca,3.5)
colormap hot

h = plotScaleBar(5e3);

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
h = plotScaleBar(1e3,1e2,'Location','NorthWest');

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

h = plotScaleBar(5e3);
h.FaceColor = 0 * ones(1,3);

if saveFigs
    mfFigureExport(cf,[outFigDir filesep 'example STORM cluster hulls multiple cells'],'DPI',600)
end

xlim(cropLimsXX_YY(1,:)),ylim(cropLimsXX_YY(2,:))

h = plotScaleBar(1e3,1e2,'Location','NorthWest','Color',[0 0 0]);

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
pts.Marker = '.';

cMap = randomColormap(numel(unique(pointInd)),44);
colormap(cMap)


h = plotScaleBar(5e3);
h.FaceColor = 0 * ones(1,3);


if saveFigs
    mfFigureExport(cf,[outFigDir filesep 'example STORM cluster assignment multiple cells'],'DPI',600)
end

xlim(cropLimsXX_YY(1,:)),ylim(cropLimsXX_YY(2,:))
h = plotScaleBar(1e3,1e2,'Location','NorthWest','Color',[0 0 0]);

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

plot(stormData.Xc,stormData.Yc,'.','color',.5 * ones(1,3))


h = plotScaleBar(5e3);
h.FaceColor = 0 * ones(1,3);

if saveFigs
    mfFigureExport(cf,[outFigDir filesep 'example STORM cluster hulls multiple cells'],'DPI',600)
end

xlim(cropLimsXX_YY(1,:)),ylim(cropLimsXX_YY(2,:))
h = plotScaleBar(1e3,1e2,'Location','NorthWest','Color',[0 0 0]);

if saveFigs
    mfFigureExport(cf,[outFigDir filesep 'example STORM cluster hulls multiple cells ROI'],'DPI',600)
end

%% --- example points only

%cf = open([exampleClusterDir filesep 'cluster convex hulls plot.fig'])
cf = figure;

plot(stormData.Xc,stormData.Yc,'.','color',.5 * ones(1,3))
axis equal
axis ij
axis off
title('')
xlabel('')
ylabel('')


h = plotScaleBar(5e3);
h.FaceColor = 0 * ones(1,3);


if saveFigs
    mfFigureExport(cf,[outFigDir filesep 'example STORM localizations multiple cells'],'DPI',600)
end

xlim(cropLimsXX_YY(1,:)),ylim(cropLimsXX_YY(2,:))
h = plotScaleBar(1e3,1e2,'Location','NorthWest','Color',[0 0 0]);

if saveFigs
    mfFigureExport(cf,[outFigDir filesep 'example STORM localizations multiple cells ROI'],'DPI',600)
end

