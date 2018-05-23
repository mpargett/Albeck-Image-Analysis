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
%   To be used only with uTrack 2.1.0a - MPargett 3.28.18

function [tracksFinal,kalmanInfoLink,p] = iman_utrack_call(movieInfo, pin)
%Version check provision
if strcmpi(movieInfo,'version'); tracksFinal = 'v2.1'; return; end
%   Set default parameters (units in frames and pixels)
p = struct('linkrad', 15,	'nnwin', 10,      ...
        'gaprad', 23,    	'gapwin', 5,        'mergeSplit', 3, ...
        'minTrkLength', 5, 	'motionType', 1,	'bStdMult', 4, ...
        'lStdMult', 4,      'gapPenalty', 1.5,  'msRadMult', 1,  ...
        'ampgap', 5, 'ampwin', 5, 'ampRat', [0.5, 1.5], ...
        'run', true);    

%Apply provided parameters
for s = fieldnames(pin)'
    if ~isempty(pin.(s{1}));	p.(s{1}) = pin.(s{1});  end
end
%   Backward compatibility horrors
if isfield(pin, 'movrad') && ~isempty(pin.movrad)
    p.linkrad = pin.movrad; p.gaprad = pin.movrad; end
if isfield(pin, 'linkwin') && ~isempty(pin.linkwin)
    p.nnwin = pin.linkwin; p.gapwin = pin.linkwin; end

%Time Window to consider
gapCloseParam.timeWindow = p.gapwin; 
%def: 7. maximum allowed time gap (in frames) between a track segment end
%   and a track segment start that allows linking them. 
%   This parameter is used for all time-based constraints

%Merging/Splitting Flag
gapCloseParam.mergeSplit = p.mergeSplit;
%1 if merging and splitting are to be considered, 2 if only merging is to
%   be considered, 3 if only splitting is to be considered, 0 if no merging
%   or splitting are to be considered.  

%Merging/Splitting search radius multiplier
gapCloseParam.msRadMult = p.msRadMult;

%Minimum Track Length
gapCloseParam.minTrackLen = p.minTrkLength;
%def: 2. minimum length of track segments from linking to be used in gap
%   closing. 

%Gap Closing amplitude measurement window
gapCloseParam.ampgap = p.ampgap;
gapCloseParam.ampwin = p.ampwin;

%optional input:
%1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.
gapCloseParam.diagnostics = 0; 


%% cost matrix for frame-to-frame linking
%Function name
costMatrices(1).funcName = 'costMatRandomDirectedSwitchingMotionLink_mp';

%Parameters
parameters.linearMotion = p.motionType; %
parameters.minSearchRadius = 2;  %Was 0
%def: 2. minimum allowed search radius. The search radius is calculated on
%   the spot in the code given a feature's motion parameters. If it happens
%   to be smaller than this minimum, it will be increased to the minimum.  
parameters.maxSearchRadius = p.linkrad; %Was 40
%def: 7. maximum allowed search radius. Again, if a feature's calculated
%   search radius is larger than this maximum, it will be reduced to this
%   maximum.  
%   From typical 1280x1080 image (0.65um pixel), 40 is ~26 um

parameters.brownStdMult = p.bStdMult;  %Was 10
%def: 3.multiplication factor to calculate search radius from standard
%   deviation.

parameters.useLocalDensity = 1; 
%1 if you want to expand the search radius of isolated features in the
%   linking (initial tracking) step. 

parameters.nnWindow = p.nnwin; 
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
costMatrices(2).funcName = 'costMatRandomDirectedSwitchingMotionCloseGaps_mp';

%parameters
parameters.linearMotion = p.motionType;  %0 Brownian, 1 Linear, 2 Linear no-reversing
parameters.minSearchRadius = 2;     %def: 2. Was 0.  Min search radius.
parameters.maxSearchRadius = p.gaprad;  %def: 7, Was 40. Max  search radius.
parameters.brownStdMult = p.bStdMult*ones(p.gapwin,1);  %Was 3*...
%multiplication factor to calculate Brownian search radius from st.dev.

parameters.brownScaling = [0.5 0.5]; 
%def: [0.25 0.01].power for scaling the Brownian search radius with time,
%   before and after timeReachConfB (next parameter). 
%   MP:  I see no reason to have behavior shift like this.  0.5 represents
%   a proper brownian model, with square root time scaling (of the standard
%   deviation of the displacement distribution)
parameters.timeReachConfB = 1; %Was 3, then p.linkwin (unknown why)
%before timeReachConfB, the search radius grows with time with the power in
%   brownScaling(1); after timeReachConfB it grows with the power in
%   brownScaling(2).    

parameters.ampRatioLimit = p.ampRat;
%def: [0.7 4].for merging and splitting. Minimum and maximum ratios between
%   the intensity of a feature after merging/before splitting and the sum
%   of the intensities of the 2 features that merge/split.  

parameters.lenForClassify = 5; 
%minimum track segment length to classify it as linear or random.
parameters.useLocalDensity = 1;  %***Consider disabling
%1 if you want to expand the search radius of isolated features in the gap
%   closing and merging/splitting step. 
parameters.nnWindow = p.nnwin; 
%number of frames before/after the current one where you want to look for a
%   track's nearest neighbor at its end/start (in the gap closing step). 
parameters.linStdMult = p.lStdMult;  %Was 1*ones(p.linkwin,1); 
%multiplication factor to calculate linear search radius from st.dev.

parameters.linScaling = [0.5 0.5]; 
%def: [0.25 0.01].power for scaling the linear search radius with time
parameters.timeReachConfL = 1;  %Was 4, then p.linkwin (as with Brownian)
%similar to timeReachConfB, but for the linear part of the motion.

parameters.maxAngleVV = 360;  %No constraint on motion directions
%def: 30.maximum angle between the directions of motion of two tracks that
%   allows linking them (and thus closing a gap). Think of it as the
%   equivalent of a searchRadius but for angles.  

%optional; if not input, 1 will be used (i.e. no penalty)
parameters.gapPenalty = p.gapPenalty; 
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
[tracksFinal,kalmanInfoLink,errFlag,p.svo] = trackCloseGapsKalmanSparse_mp(movieInfo,...
    costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);

if nargout > 2
    p.p.costMatrices = costMatrices;
    p.p.gapCloseParam = gapCloseParam;
    p.p.kalmanFunctions = kalmanFunctions;
    p.p.probDim = probDim;
    p.p.saveResults = saveResults;
    p.p.verbose = verbose;
end
end

