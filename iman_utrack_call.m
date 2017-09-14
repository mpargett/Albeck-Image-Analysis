%% Gap closing parameters and call for U-Track
%
% Copyright (C) 2013 LCCB 
%
% This file is part of U-Track.
% 
% U-Track is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% U-Track is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with U-Track.  If not, see <http://www.gnu.org/licenses/>.
% 
%   Modified by MPargett, 5.8.15 
%

function tracksFinal = iman_utrack_call(movieInfo, p, op)
%Version check provision
if strcmpi(movieInfo,'version'); tracksFinal = 'v1.0'; return; end
%   Assign provided parameters
pdef.mergeSplit = 3;
params = fieldnames(pdef);
for s = 1:numel(params)
   if ~isfield(op.trk, params{s})
       op.trk.(params{s}) = pdef.(params{s});
   end
end

%INFORMATION NEEDED:
%   Time Step
%   Image Size / Pixel size

%Global Parameters to employ throughout this procedure
tWin = op.trk.linkwin;    %Time window of interest about a frame
                              % Scaled on typical sampling time of 5 min
iRad = op.trk.movrad/sqrt(p.PixSizeX.^2 + p.PixSizeY.^2);   %Radius of interest
%   About a particle (cell centroid), scaled for 25 um


%Time Window to consider
gapCloseParam.timeWindow = tWin; 
%def: 7. maximum allowed time gap (in frames) between a track segment end
%   and a track segment start that allows linking them. 
%   This parameter is used for all time-based constraints

%Merging/Splitting Flag
gapCloseParam.mergeSplit = op.trk.mergeSplit;
%1 if merging and splitting are to be considered, 2 if only merging is to
%   be considered, 3 if only splitting is to be considered, 0 if no merging
%   or splitting are to be considered.  

%Minimum Track Length
gapCloseParam.minTrackLen = 3;
%def: 2. minimum length of track segments from linking to be used in gap
%   closing. 

%optional input:
%1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.
gapCloseParam.diagnostics = 0; 


%% cost matrix for frame-to-frame linking
%Function name
costMatrices(1).funcName = 'costMatRandomDirectedSwitchingMotionLink';

%Parameters
parameters.linearMotion = 0; %use linear motion Kalman filter.
parameters.minSearchRadius = 2;  %Was 0
%def: 2. minimum allowed search radius. The search radius is calculated on
%   the spot in the code given a feature's motion parameters. If it happens
%   to be smaller than this minimum, it will be increased to the minimum.  
parameters.maxSearchRadius = iRad; %Was 40
%def: 7. maximum allowed search radius. Again, if a feature's calculated
%   search radius is larger than this maximum, it will be reduced to this
%   maximum.  
%   From typical 1280x1080 image (0.65um pixel), 40 is ~26 um

parameters.brownStdMult = 6;  %Was 10
%def: 3.multiplication factor to calculate search radius from standard
%   deviation.

parameters.useLocalDensity = 1; 
%1 if you want to expand the search radius of isolated features in the
%   linking (initial tracking) step. 

parameters.nnWindow = tWin; 
%number of frames before the current one where you want to look to see a
%   feature's nearest neighbor in order to decide how isolated it is (in
%   the initial linking step).  

%Kalman filter initialization parameters.
parameters.kalmanInitParam = []; 

%Kalman filter initialization parameters.
% parameters.kalmanInitParam.searchRadiusFirstIteration = 10; 

%optional input
parameters.diagnostics = []; 
%if you want to plot the histogram of linking distances up to certain
%   frames, indicate their numbers; 0 or empty otherwise. Does not work for
%   the first or last frame of a movie.  

costMatrices(1).parameters = parameters;
clear parameters

%% cost matrix for gap closing
%function name
costMatrices(2).funcName = 'costMatRandomDirectedSwitchingMotionCloseGaps';

%parameters
parameters.linearMotion = 0;  %0 Brownian, 1 Linear, 2 Linear no-reversing
parameters.minSearchRadius = 2;     %def: 2. Was 0.  Min search radius.
parameters.maxSearchRadius = iRad;  %def: 7, Was 40. Max  search radius.
parameters.brownStdMult = 6*ones(tWin,1);  %Was 3*...
%multiplication factor to calculate Brownian search radius from st.dev.

parameters.brownScaling = [0.25 0.25]; 
%def: [0.25 0.01].power for scaling the Brownian search radius with time,
%   before and after timeReachConfB (next parameter). 
%   MP:  I see no reason to have behavior shift like this
parameters.timeReachConfB = 1; %Was 3, then tWin (unknown why)
%before timeReachConfB, the search radius grows with time with the power in
%   brownScaling(1); after timeReachConfB it grows with the power in
%   brownScaling(2).    

parameters.ampRatioLimit = [0.4 3];
%def: [0.7 4].for merging and splitting. Minimum and maximum ratios between
%   the intensity of a feature after merging/before splitting and the sum
%   of the intensities of the 2 features that merge/split.  

parameters.lenForClassify = 5; 
%minimum track segment length to classify it as linear or random.
parameters.useLocalDensity = 1; 
%1 if you want to expand the search radius of isolated features in the gap
%   closing and merging/splitting step. 
parameters.nnWindow = tWin; 
%number of frames before/after the current one where you want to look for a
%   track's nearest neighbor at its end/start (in the gap closing step). 
parameters.linStdMult = 6;  %Was 1*ones(tWin,1); 
%multiplication factor to calculate linear search radius from st.dev.

parameters.linScaling = [0.25 0.25]; 
%def: [0.25 0.01].power for scaling the linear search radius with time
parameters.timeReachConfL = 1;  %Was 4, then tWin (as with Brownian)
%similar to timeReachConfB, but for the linear part of the motion.

parameters.maxAngleVV = 360;  %No constraint on motion directions
%def: 30.maximum angle between the directions of motion of two tracks that
%   allows linking them (and thus closing a gap). Think of it as the
%   equivalent of a searchRadius but for angles.  

%optional; if not input, 1 will be used (i.e. no penalty)
parameters.gapPenalty = 1; 
%def: 1.5. penalty for increasing temporary disappearance time
%   (disappearing for n frames gets a penalty of gapPenalty^n). 
%   MP: Uncertain what this penalizes (fitness to link for gap closure?)

%optional; to calculate MS search radius
%if not input, MS search radius will be the same as gap closing search radius
parameters.resLimit = []; 
%resolution limit, which is generally equal to 3 * point spread function
%   sigma. 

costMatrices(2).parameters = parameters;
clear parameters

%% Kalman filter function names
kalmanFunctions.reserveMem  = 'kalmanResMemLM';
kalmanFunctions.initialize  = 'kalmanInitLinearMotion';
kalmanFunctions.calcGain    = 'kalmanGainLinearMotion';
kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';

%% additional input
%   saveResults.filename = 'tracksTest1DetectionAll1.mat'; 
%   name of file where input and output are saved
saveResults = 0;    %don't save results
verbose = 0;        %verbose state
probDim = 2;        %problem dimension (2D vs. 3D images)


%% tracking function call
[tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo,...
    costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);

end

