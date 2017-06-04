%% This script is going to try to correlate dikes with notebooks
% clear workspace

%%%%% GOAL %%%%%
%%% THIS SCRIPT WILL, try to find the geographic names that correspond with
%%% the journal names.  
clc
clear all
close all

cd('/Users/Matthew/GitHub/Notebook_matching');
load preppedNotebookData.mat
load geographicData.mat

%% STEP 1, look for common names between geoData and raw


% Create a series of loops through each row of jrnData, then looping
% through all of the geographic names, looking for a match

% %set the number of characters to extract out of the journal local and use
% to determin the similarity between geographic name and journal name
matchThresh = 4; 

geoMatch = zeros(length(jrnData)*10,2); 
% first column: index within Journal Data
% second column: index within Geographic Names
a = 1;

h = waitbar(0, 'Please Be patient... wait .... keep waiting...');

for i = 2:length(jrnData) %Step a jump into the journal Data
    locale = char(jrnData{i,15}); %General Location Name stored in local variable
    
    waitbar(a/numel(jrnData));
%     first test if there is an Nan, length = 1
    if length(locale) < matchThresh
        flag = 1;
    else
        % create a subset of the full locale name
        localeSub = locale(1:matchThresh);
        flag = 0;
    end
        
    % create a logical statement to test for NaN conditional flag
    if flag == 1
        continue
        fprintf('oops');
    else
        for j = 1:length(geoData) %Step into the geographic names
            geoName = char(geoData{j,4});g
            %If statement comparing all geographic names to the stored local.


    %         This step compares a subset of the locale text
            if contains(geoName, locale,'IgnoreCase',true)

                geoMatch(a,:) = [i j];
                a = a + 1;
            else
            end




        end
    end
    
end

close(h)

geoMatch( ~any(geoMatch,2), : ) = [];  %delete zeroes from rows 
save('possibleGeographicMatches.mat','geoMatch')


%%%% As of June 2, this takes about an hour
%% Step next.  elevation relationship between mapped dike and geographic named place
%  There is an elevation for each notebook dike.  There is an elevation for
%  each of the geographic named features.

% Set a threshold elevation relationship so that they are close enough to
% to each other that if they're outside a 300 or 500 m elevation range I
% should throw out that match as not helpful.

% declare an elevation threshold for how close the journal and geographic
% features should be to each other.
elevThresh = 200; %200 m 
geoMatchClose = zeros(3500,2);

h = waitbar(0, 'Please Be patient... wait .... keep waiting...');
a = 1;
for i = 1:size(geoMatch,1);
    geoElev = geoData{geoMatch(i,2),3};
    jrnElev = jrnData{geoMatch(i,1),8};
    
    dist = abs(geoElev - jrnElev);
    
    
    waitbar(i/size(geoMatch,1));
    if dist <= elevThresh
        geoMatchClose(a,:) = geoMatch(i,:);
        a = a + 1;
        
    else
        continue
    end
end
geoMatchClose( ~any(geoMatchClose,2), : ) = [];  %rows
close(h)

save('GeographicNameMatchesSimilarElevations.mat','geoMatchClose');
%% create a cell array for the matched names

%           Col 1              ||           COL 2
% geographic journal name data ||  geographic database names



for i = 1:length(geoMatchClose)
    
    
    geoNames(i,:) = {jrnData{geoMatchClose(i,1),15} char(geoData{geoMatchClose(i,2),4})};
    
    
end

% Part of the problem at this stage is the lack of uniqueness. There are
% many repetitions in both the names of the locations and there are
% multiple places within the geograhic dataset that have the same name 

% So now there are nearly 35,000 matches between the geographic name
% dataset and the notebooks.
save('possibleGeographicMatchesNames.mat','geoNames');



%%
%clean up the workspace
clearvars -except jrnData geoData geoMatchClose strmData nonLinPts curvPts linPts
