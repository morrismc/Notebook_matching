%Notebook loader

%  This script will load and format all of the Taubeneck noteobok data
%  and then look to save that data as a .mat file.  

% The script will also load geographic data compiled from ArcGIS (streams
% and geographic names)
clear all
clc
%%
cd('/Users/Matthew/GitHub/Notebook_matching');
%%

%load dike data
[curvPts curvTxt] = xlsread('CurvedXYElev.xlsx');
[linPts linTxt] = xlsread('LinearXYElev.xlsx');
[nonLinPts nonLinTxt] = xlsread('NonLinearXYElev.xlsx');

linPts(linPts < 0) = NaN; % turn all those pesky negatives into some NaNs.

%load geographic feature names
[geoPts geoTxt] = xlsread('GeographicFeatures.xlsx');
geoTxt = geoTxt(2:end,1); % delete the header

%load creak points
[strPts strmTxt] = xlsread('RiverCreekPtsFINAL.xlsx');
strmTxt = strmTxt(2:end,1); % delete header


%% format all of these features into cell arrays?

% at this point, all of the dike data can remain as matrices; however, the
% stream data and geographic name data needs to be reconciled with the
% strings and numeric data.

% All mapped dike data formated in the following way:
% Cl 1|| Col 2   || COl 3   ||    COL 4       ||    COL 5

% FID || Point X || Point Y || Elevation (FT) || Elevation (m) 



%now store the geographic name data



for i = 1:length(geoPts)

    geoData{i,1} = geoPts(i,1); %x data
    geoData{i,2} = geoPts(i,2); % y data
    geoData{i,3} = geoPts(i,3); % elevation data (m)
    geoData{i,4} = geoTxt(i,1); % name of the geographic point
end

% GeoData now contains all of the name, elevation, y and x data for
% goegraphic points



% now store the stream name data
strmData = cell(length(strPts),1);

for i = 1:length(strPts)
    strmData{i,1} = strPts(i,2); %x data
    strmData{i,2} = strPts(i,3); % y data
    strmData{i,3} = strPts(i,1); % elevation data (m)
    strmData{i,4} = strmTxt(i,1); % name of the geographic point
end


% All stream and geographic data is formatted the following way:
% Col 1   || COl 2   ||    COL 3       ||    COL 4

% Point X || Point Y || Elevation (m)  ||    NAME 

clearvars -except strmData geoData curvPts linPts nonLinPts 

save geographicData.mat

%%  Import Journal Data

% Import data from the Taubeneck Journal Data, Storing only that
% information that is important for this analysis, deleting empty rows etc.

[~, ~, raw] = xlsread('Journal_data_Full_copy.xlsx');

% raw now contains all of the data from the excel spreadsheet.  

% some preprocessing that needs to take place before this data is ready to
% be handled
% 1) if there is no elevation data for a dike, set row to NaN
% 2) If there is no geographic name, set row to NaN


for i = 2:length(raw) %from 2 to avoid header
    if isnan(raw{i,8}) %test if the elevation is empty
        raw(i,:)  = {NaN};
    elseif ~any(raw{i,8}) % test if the elevation is zero
        raw(i,:)  = {NaN}; % if the elevation is zero, set the row to NaN
    elseif isnan(raw{i,15})
        raw(i,:)  = {NaN};      
    end
end

jrnData = raw;
% Raw now contains an edited version of the table of dike measurements from
% the notebooks, Only those dikes that have both an elevation and a
% geographic name. 

save('preppedNotebookData.mat','jrnData')

%%
% Now all of the data needed to complete the correlation betwee dikes and
% notebooks is sorted and arrnanged

clear all
clc