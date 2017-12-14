function [allErrs,times] = batchClusterStormData(parentDir)

% 
% ip = inputParser;
% ip.addParamValue('CropData',false,@islogical);
% 
% ip.parse(varargin{:});
% p = ip.Results;


if nargin < 1 || isempty(parentDir)    
    parentDir = uigetdir(pwd,'Select a directory with files to cluster:');    
end

allFiles = dir([parentDir filesep '*.txt']);

nFiles = numel(allFiles);

disp([num2str(nFiles) ' STORM files found...'])


allErrs = cell(nFiles,1);

disp('Running clustering on all files, please wait...')

times = nan(nFiles,1);
for j = 1:nFiles
        
   try
    
        currFile = [parentDir filesep allFiles(j).name];
        tic;        
        clusterStormData(currFile);
        times(j) = toc;
   catch err
       
       disp(['error on file ' num2str(j) ' : ' err.message])
       
       allErrs{j} = err;
   end
end