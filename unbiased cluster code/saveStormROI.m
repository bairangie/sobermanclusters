function saveStormROI(varargin)
%SAVESTORMROI allows manual creation and saving of an ROI for a STORM file
%
% saveStormROI
% saveStormROI('StormFile',stormFileName)
% saveStormROI('OptionName1',optionValue1,...)
%
% This lets the user manually select areas of interest in an input storm
% file. The resulting ROI is then saved in a format which can be read by
% clusterStormROI for further processing.
%
% Note: Do not rename or re-arrange the storm .txt files and ROIs after ROI
% creation - this can break the link between the roi and corresponding
% storm .txt data file
%

%Hunter Elliott
%Sometime in 2013?


%% ----------- Input ------------ %%

oExt = 'sroi';

ip = inputParser;
ip.addParamValue('StormFile','',@(x)(isempty(x) || exist(x,'file')));%Name of STORM file to select ROI for. Default is to ask user.
ip.addParamValue('ROI',[],@(x)(isempty(x) || (size(x,2) == 2 && size(x,1) >= 3)));%Can directly input an ROI polygon.
ip.addParamValue('FramesUse',[],@(x)(isempty(x) || all(isposint(x))));%Frames to use data from. Empty means use all frames
ip.addParamValue('OutputFile','',@ischar);%Name of file to save ROI to. Default is to ask user.
ip.addParamValue('ExistingROIFile','',@(x)(isempty(x) || (exist(x,'file')>0 && strcmp(x(end-4:end),'.sroi')) ) );%Can optionally specify a Whether to use a previously saved ROI polygon but apply to a different storm .txt file.
ip.addParamValue('UseExistingROI',false,@islogical);%If true, an ROI can be loaded from an existing file. If 'ExistingROIFile' is set that file will be used, if not the user will be asked to select. If this is false then ExistingROIFile should not be specified.
ip.parse(varargin{:});

p = ip.Results;

[filePath,fileName] = optionalFileInput(p.StormFile,'*.txt','Select a STORM file:');

if filePath == 0
    return
end

    
if ~isempty(p.ExistingROIFile) && ~p.UseExistingROI
    error('Conflicting options!! - The UseExistingROI option was set to false but an ExistingROIFile was input!') %Avoid user error/confusion.
elseif p.UseExistingROI
    
    [roiPath,roiFile] = optionalFileInput(p.ExistingROIFile,'*.sroi','Select an ROI file to apply:');
    if filePath == 0
        return
    end
    
    p.ExistingROIFile = [roiPath roiFile];
    
    %Load the ROI polygon
    roiDat = load(p.ExistingROIFile,'-mat');
    p.ROI = roiDat.cropPoly;
    
end    
    
    
p.StormFile = [filePath fileName];

outVars = {'p'};%List of variables to write to file

outVars = [outVars { 'fileName', 'filePath'}];

%% ----- Loading ---- %%

disp('Reading STORM data file...')
stormData = readSTORMTxtFile(p.StormFile);
if isempty(stormData)
    %If the user clicked "cancel"
    return
end

%Display only selected frames for cropping
if ~isempty(p.FramesUse)
    stFields = fieldnames(stormData);
    nPt = numel(stormData.X);
    goodFrame = ismember(stormData.Frame,p.FramesUse);
    nFields = numel(stFields);
    for j = 1:nFields
        if numel(stormData.(stFields{j})) == nPt            
            stormData.(stFields{j}) = stormData.(stFields{j})(goodFrame);
        end
    end
end

%% ----- Cropping and Output ---- %%

[~,cropPoly] = cropStormData(stormData,p.ROI);     %#ok<NASGU>

outVars = [outVars {'cropPoly'}];

dateTime = datestr(now); %#ok<NASGU>

if isempty(p.OutputFile)
    [outFileName,outFilePath] = uiputfile(['*.' oExt],'Select where to save the ROI:',[p.StormFile(1:end-4), '_ROI_1.', oExt]);
else
    [outFilePath,outFileName] = fileparts(p.OutputFile);
    if isempty(outFilePath)
        outFilePath = pwd;
    end
    outFilePath = [outFilePath filesep];
    outFileName = [outFileName '.' oExt];
end

%TEMP - file separator consistency?

outVars = [outVars {'dateTime', 'outFilePath','outFileName'}];

save([outFilePath outFileName],outVars{:});