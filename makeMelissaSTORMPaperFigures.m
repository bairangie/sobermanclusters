%Instructions:
%
%You can run a section by selecting it (clicking anywhere within it) and
%then pressint Ctrl+Enter, or by clicking the "Run Section" button in the
%Editor tab.
%
%First, run the "Setup" section. It will ask you to select a directory that
%contains all the combined_stats.mat files for each condition.
%
%Then, set the iStatFields number to choose your statistic (number,
%area or density) - see instructions within section.
%
%Then, run one of the other sections below to make either histogram
%overlays or whisker plots.
%
%You will be asked to select one or more of the combined_stats.mat from one
%or more condition to compare. Hold down Ctrl and click on multiple
%conditions to select the conditions to compare.
%
%



%% ===SETUP: RUN THIS CELL FIRST TO SET THINGS UP ====%%

%Specify the directory containing all the combined_stats.mat files for each
%condition. This will be the contents of the
%All_Melissa_Combined_Cluster_Stats_for_Angie_4_17_2017.zip file I sent you
allCondDir = uigetdir(pwd,'Select a directory containing all of the combined cluster stat files:')

condStr   = {'5LO','FLAP','Control','cPLA2 Inh','_Cell','Membrane','0 min','2 min','5 min','10 min','Ag','DMSO','IgE'};  %Strings to search for in the file names
statFileNames = {'5LO','FLAP','Control','+ cPLA2 Inh','Cell','Membrane','0 min','2 min','5 min','10 min','Ag','DMSO','IgE'}; %Strings for figure / output naming



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ===== RUN THIS SECTION TO MAKE Histogram overlays ==== %%

%RUN THIS CELL TO CREATE POINT-WEIGHTED HISTOGRAM OVERLAYS


iStatFields = 1;%<======Set to 1 for cluster number of localizations, 2 for cluster area and 3 for cluster density


%Make histogram overlays
stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',statFileNames,...
                'iStatField',iStatFields,'DataType','Point-Weighted','FigureType','FilledHistogram')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ===== RUN THIS SECTION TO MAKE WHISKER PLOTS ==== %%
            
        
%RUN THIS CELL TO CREATE POINT-WEIGHTED WHISKER PLOTS

iStatFields = 1;%<======Set to 1 for cluster number of localizations, 2 for cluster area and 3 for cluster density


stormCombinedStatsAndConditionComparisonFigures(allCondDir,'condStr',condStr,'condNames',statFileNames,...
                'iStatField',iStatFields,'DataType','Per-ROI-Mean','FigureType','Box','StatToCalc','Median',...
                'PctLimPlot',[0 95],'PlotArgs',{'notch','on'});






