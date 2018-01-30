%IMAN_CELLID
%   Identify (segment) cell nuclei based on size and shape (roundness)
%   over a series of intensity thresholded binary masks.
%
%   [movieInfo, enuc] = iman_cellid(im, op, bkg)
%       return movieInfo (structure formatted for input into uTrack, with
%   fields xCoord, yCoord, amp [sum of signal], are [area of nucleus])
%   and enuc (a binary mask of nuclei).  The procedure takes as inputs im,
%   the image to be segmented, op, a paramter structure, and bkg background
%   values for each channel of the image.
%
%   The op structure may contain the following field:
%       chan        -Index of the channel to use for segmentation
%       cyt         -Logical - TRUE if segmentation channel is cytosolic
%       usegrad     -Logical - TRUE to use gradients (depricated, for
%                    cytoplasmic images, consider localscale)
%       localscale  -Logical - TRUE to scale by local intensity
%       minD        -Minimum Diameter of acceptable nucleus shapes (pixels)
%       maxD        -Maximum Diameter of acceptable nucleus shapes (pixels)
%       maxEcc      -Maximum Eccentricity of acceptable nucleus shapes
%                   (Eccentricity is deviation from circular shape, [0-1] 
%                    0 for a circle, 1 for a line)
%       Extent      -Minimum and Maximum 'Extent' values of acceptable
%                    nucleus shapes (Extent is the ratio of the shape area
%                    to the area of its bounding box, rotated to match the
%                    shape's orientation, [0-1], value for a circle or
%                    ellipse is PI/4.
%       sigthresh   -Signal value to use as a foreground threshold
%       hardsnr     -Logical - TRUE is sigthresh is a hard cutoff
%       nth         -Number of thresholds to use for segmentation
%

%FIXME consider better nuclear spot ID procedures
%   Find inner and outer spots (for nuc, and inside of cyto)?
%   Use multiple channels!!!

function [movieInfo, enuc] = iman_cellid(im, pin, bkg)
%Version check provision
if strcmpi(im,'version'); movieInfo = 'v2.1'; return; end

%% Operation parameters
%Default parameters
p = struct('chan', 1, 'cyt', false, 'minD', 8, 'maxD', 25, ...
           'maxEcc', 0.9, 'Extent', [0.65,0.9], 'sigthresh', [], ...
           'hardsnr', false, 'nth', 20, 'nrode', 0, ...
           'localscale', false, 'usegrad', false);

%Apply provided parameters
for s = fieldnames(pin)'; 	p.(s{1}) = pin.(s{1});  end

bkg_rat = 0.5;          %Min foreground to background ratio (~SNR-1)
%   Get image background intensity
if ~isempty(bkg); stim_bkg = bkg(p.chan); else stim_bkg = 0; end
%   Get background ratio for a threshold on fore- vs. back-ground
if ~isempty(p.sigthresh)
    %   Use a conservative estimate of 1/2 the SNR as threshold
    %   unless directed to use hard SNR threshold (op.seg.hardsnr)
    bkg_rat = max(0, (p.sigthresh./stim_bkg - 1)*(1 + p.hardsnr)/2 );
end

%% Preliminary calculations
%Calculations used for cell sizing
maxNucArea = round( (pi*(p.maxD - p.nrode*2)^2)/4 );
minNucArea = round( (pi*(p.minD - p.nrode*2)^2)/4 );
flts = max(1,floor(p.maxD/10));   %Filter size, based on expected nucleus size

%Prefiltering now depricated 2018-01-18
    %Define a Gaussian filter (smoothes noise)
    %gaussianFilter = fspecial('gaussian', flts.*[1, 1], flts/2);
    %Filter only the target segmentation channel
    % stim = double( imfilter(im(:,:,p.chan), gaussianFilter, 'replicate') );
%   All further operations consider this segmentation target image (stim)
stim = double( im(:,:,p.chan) );
%Set 'foreground' region based on background levels (previously removed)
stim_fore = stim > stim_bkg.*(bkg_rat);   %Foreground area mask
%Define binary structuring elements
morphflt = ceil(flts/4);
st1 = strel('disk', morphflt);  
st2 = strel('disk', 2*morphflt + p.nrode);


%% Image filtering (IF necessary, i.e. using edge filters for cyto)
%Gradient filtering only if requested - 2018-01-18
if p.usegrad
        %   Filter out background by setting it to uniform level
        stim(~stim_fore) = stim_bkg.*(bkg_rat);
        %Create a predefined 2d filter. "Sobel edge-emphasizing filter"
        hy = fspecial('sobel');   hx = hy';
        %Filter image for edges (gradients, via the Sobel filter)
        tstim = sqrt(imfilter(stim, hx, 'replicate').^2 ...  %Across x axis.
            + imfilter(stim, hy, 'replicate').^2);     %Across y axis.
        %   Add image to gradients (penalizes 'high' regions for cyto)
        tstim = stim + tstim;
    else    tstim = stim;  %Otherwise, just use image
end

%Cytoplasmic segmentation inverts the image (after scaling if req'd)
if p.cyt 
    %   IF using local scaling, dilate to overlap nucleus and scale
    if p.localscale
        stl = strel('disk',floor(p.maxD./2));
        tstim = tstim ./ imerode(imdilate(stim,stl),stl); 
    end    
    %Invert to make image/gradients 'low'
    stim = max(tstim(:)) - tstim;
else %IF nuclear segmentation, perform scaling very locally
    if p.localscale; stim = tstim ./ imdilate(stim,st2); end
    stim(~stim_fore) = 0;   %Neglect non-foreground aggresively
end;    clear tstim;


%% Image thresholding search
%Define thresholds at linearly spaced quantiles
thresholds = prctile(stim(stim_fore), linspace(5, 95, p.nth) );
clear stim_fore;

%Revision:  For each threshold, store nucleus estimates
%   Compare after all are generated and select 'best' for each spot

%FOR each threshold, evalute appearance of nucleus-like shapes
enuc = false(size(stim)); 
for s = 1:length(thresholds)
    %Threshold the gradient image, making a temporary image (tim)
    tim = stim > thresholds(s); %Binary, based on threshold value
    
    %Dilate and erode to eliminate salt-and-pepper noise in binary
    tim = imdilate(tim, st1); %Dilate (remove holes)
    tim = imerode(tim, st2);  %Erode (remove spots, shrink if p.nrode)
    
    %Get properties of each region in image (Area, Perimeter, ...)
    S = regionprops(tim, 'Area', 'Eccentricity', ...
        'Extent', 'Orientation', 'PixelIdxList', ...
        'MajorAxisLength', 'MinorAxisLength');
    nucArea  = cat(1,S.Area);           %Vector of Areas
    nucOr  = cat(1,S.Orientation).*pi./180; %Vector of Orientations    
    a2 = cat(1,S.MajorAxisLength).^2;       %Vector of Major Axes^2
    b2 = cat(1,S.MinorAxisLength).^2;       %Vector of Minor Axes^2
    
    %   Calculate adjusted nuclear extent, based on orientation
    nucExt = cat(1,S.Extent) .* sqrt( 1 + (a2./b2 + b2./a2 - 2)...
        .*(sin(nucOr).*cos(nucOr)).^2 );
    
    %Determine objects passing size/shape fiters
    scorePass  = nucArea > minNucArea & nucArea < maxNucArea &...
        cat(1,S.Eccentricity) < p.maxEcc & ...
        nucExt > p.Extent(1) & nucExt < p.Extent(2);
    
    %Store passing spots (overwrite older, which are always smaller)
    for ss = find(scorePass)'; enuc(S(ss).PixelIdxList) = true; end

end

%Restore nuclear size from hole removal 
%   (does not restore full size if p.nrode > 0)
enuc = imdilate(enuc, st1);


%% Store estimated nuclei data
%Re-segment stored nuclei spot (re-orders spots)
S = regionprops(enuc, 'Centroid', 'Area', 'PixelIdxList');
nCtr = cat(1, S.Centroid);  nAre = cat(1, S.Area);  
im_op = im(:,:,p.chan);
nAmp = arrayfun(@(x)sum(im_op(x.PixelIdxList)), S);

%Store coordinates IF any spots were found
if ~isempty(nCtr)
    movieInfo.xCoord = nCtr(:,1);
    movieInfo.yCoord = nCtr(:,2);
    movieInfo.amp = nAmp(:,1);
    movieInfo.are = nAre(:,1);
else  %Store empty output if no spots found
    movieInfo = struct('xCoord',[],'yCoord',[],'amp',[], 'are',[]);
end
% 
