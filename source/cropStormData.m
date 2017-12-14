function [stormData,maskPoly,isIn] = cropStormData(stormData,maskPoly)
%CROPSTORMDATA lets a user select a polygon ROI to crop out a region of the input storm data
%
% [stormData,maskPoly,isIn] = cropStormData(stormData)
% [stormData,maskPoly,isIn] = cropStormData(stormData,maskPoly)
%
%Hunter Elliott
%Sometime in late 2012?

%If an ROI wasn't input, let the user create one
if nargin < 2 || isempty(maskPoly)
    
    maskFig = fsFigure(.5);
    plot(stormData.X,stormData.Y,'.k')
    hold on    
    axis image    
    title({'Please select the localizations to include:','**Click OK when finished**'})
    uicontrol('Style','pushbutton','String','OK','Position',[10 10 40 20],'Callback','uiresume(gcbf)');    
    pHan = impoly(get(maskFig,'CurrentAxes'));          
    uiwait(maskFig);
    maskPoly = pHan.getPosition;
            
end


%Find point indices within this polygon
isIn = inpolygon(stormData.X,stormData.Y,maskPoly(:,1),maskPoly(:,2));

sdFields = fieldnames(stormData);
nFields = numel(sdFields);
for j = 1:nFields
    stormData.(sdFields{j}) = stormData.(sdFields{j})(isIn);
end
    
if exist('maskFig','var') && ishandle(maskFig)
    close(maskFig)
end
    
    


