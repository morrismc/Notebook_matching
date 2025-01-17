%%% Matthew Morriss June, 2017

%%%% GOAL FOR THIS SCRIPT %%%%%
%%%% This script will take the exported data from Notebook_match_1 and then
%%%% look to correlate the actual mapped dikes now with the journal
%%%% observations THROUGH the geographic names


%%% LOGIC %%%%
    %%% --> Notebook_Match_1 identified the journal entries that were close in
    %%% elevation and name to geogrpahic named points. 
    
    %%% --> STEP 1 Mapped Dikes have elevation points along them.  Identify the
    %%% elevation range for the mapped dike points and find each dike that
    %%% is within an elevation threshold of each geographic named point.
    %%% (CONNECT Mapped DIKES TO geographic elevation)
    
    %%% -->  STEP 2 Then, look to minimize the X,Y distance between those dikes that
    %%% are within this threshold of the elevation.  
    %%% (CONNECT mapped DIKES to geographic X,Y)
    
    %%% --> STEP 3 There will be some subset of the mapped dikes that are close in
    %%% elevation and X,Y distance to geographic named places.  Then cross
    %%% correlate between those mapped dikes, geographic points, and
    %%% journal observed points.
    
%% Find dikes that are within the elevation range of geographic points
%  IMPORTANT POINT, to simplify the numbers involved only look for dikes
%  that are close to geographic points that have been matched with journal
%  observations, or column 2 of geoMatchClose


% All mapped dike data formated in the following way:
% Cl 1|| Col 2   || COl 3   ||    COL 4       ||    COL 5

% FID || Point X || Point Y || Elevation (FT) || Elevation (m) 


possibleDikeGeographicPtMatch = zeros(size(linPts,1)*1000,2);

% set elevation threshold
elevThresh = 100; % 100 m

h = waitbar(0, 'Please Be patient... wait .... keep waiting...');

a = 1;

for i = 1:size(linPts,1) % enter the linear dike points array
    dikePtElev = linPts(i,5); %store each points elevation
    
    
    waitbar(i/size(linPts,1));
    
    %enter the geographicName Match array
    for j = 1:size(geoMatchClose,1)
        geoPtElev = geoData{geoMatch(j,2),3};
        
        dist = abs(dikePtElev - geoPtElev);
        
        if dist <= elevThresh
            possibleDikeGeographicPtMatch(a,:) = [i j];
            a = a + 1;
        else
        end
        
        
    end
    
    
end

close(h)

%clip off zeros
possibleDikeGeographicPtMatch( ~any(possibleDikeGeographicPtMatch,2), : ) = [];  %rows


% Format for Possible Dike Geographic Match Point (based on elevation)
%   Col 1   ||     Col 2       ||

% dikePt IX || Geographic Pt IX|| 
save('PossibleDikeGeographicPointMatch.mat','possibleDikeGeographicPtMatch');
%%  Minimize distance between the geographic points and the mapped dikes

%%% This section will use a euclidean distance to minimize the distance
%%% between the mapped dike segments and geographic points. 
