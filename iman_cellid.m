%IMAN_CELLID
%   Image Analysis Cell Identification procedure.
%
%FIXME get a real header here

%FIXME consider better nuclear spot ID procedures
%   Find inner and outer spots (for nuc, and inside of cyto)?
%   Use multiple channels!!!

function [movieInfo, enuc] = iman_cellid(im, op, bkg)
%Version check provision
if strcmpi(im,'version'); movieInfo = 'v2.0'; return; end

%% Operation parameters
nth = 20;               %Number of Thresholds for binary operations
erode_grad = floor(op.seg.minD./4); %Erosion depth for gradients (cyto segment)
bkg_rat = 0.5;          %Min foreground to background ratio (~SNR-1)
%   Get flag for hard (or soft) SNR thresholding
if isfield(op.seg, 'hardsnr') && ~isempty(op.seg.hardsnr); hsnr = true;
else hsnr = false;  end

%   Get image background intensity
if ~isempty(bkg); stim_bkg = bkg(op.seg.chan); else stim_bkg = 0; end
%   Get background ratio for a threshold on fore- vs. back-ground
if isfield(op.seg, 'sigthresh') && ~isempty(op.seg.sigthresh)
    %   Use a conservative estimate of 1/2 the SNR as threshold
    %   unless directed to use hard SNR threshold (op.seg.hardsnr)
    bkg_rat = max(0, (op.seg.sigthresh./stim_bkg - 1)*(1+hsnr)/2 );
end

%% Preliminary calculations
%Calculations used for cell sizing
maxNucArea = round(pi*op.seg.maxD^2/4);
minNucArea = round(pi*op.seg.minD^2/4);
flts = max(1,floor(op.seg.maxD/10));   %Filter size, based on expected nucleus size

%Define a Gaussian filter (smoothes noise)
gaussianFilter = fspecial('gaussian', flts.*[1, 1], flts/2);
%Filter only the target segmentation channel
%   All further operations consider this segmentation target image (stim)
stim = double( imfilter(im(:,:,op.seg.chan), gaussianFilter, 'replicate') );
%Set 'foreground' region based on background levels (previously removed)
stim_fore = stim > stim_bkg.*(bkg_rat);   %Foreground area mask
%Define binary structuring elements
st1 = strel('disk', ceil(flts/4));  st2 = strel('disk', 2*ceil(flts/4));


%% Image filtering (IF necessary, i.e. using edge filters for cyto)
if op.seg.cyt
    %Create a predefined 2d filter. "Sobel edge-emphasizing filter"
    hy = fspecial('sobel');   hx = hy';
    %Filter image for edges (gradients, via the Sobel filter)
    stim = sqrt(imfilter(stim, hx, 'replicate').^2 ...  %Across x axis.
             + imfilter(stim, hy, 'replicate').^2);     %Across y axis.
    stim = max(stim(:)) - stim;         %Invert to make gradients 'low' 
end

%% Image thresholding search
%Define thresholds at linearly spaced quantiles
thresholds = quantile(stim(stim_fore), linspace(0.05, 0.95, nth) );
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
    tim = imerode(tim, st2);  %Erode (remove spots)
    tim = imdilate(tim, st1); %Redilate (restore size)
    
    %Get properties of each region in image (Area, Perimeter, ...)
    S = regionprops(tim, 'Area', 'Perimeter', 'PixelIdxList');
    nucArea  = cat(1,S.Area);                      %Vector of Areas
    nucPerim = cat(1,S.Perimeter);                 %Vector of Perimeters
    nucFormfactor = 4*pi*nucArea./(nucPerim.^2);   %Area/Perim ratio
    
    %Calculate continuous scores to rank spots
    %   Size target is between min and max Areas
    szScore = (nucArea - minNucArea).*(maxNucArea - nucArea);
    %   Shape target is most circular (with minimum allowable value) 
    shScore = nucFormfactor - op.seg.minF;
    %Evaluate binary decision from size and shape scores (Perimeter)
    scorePass  = szScore >= 0 & shScore >= 0;
    
    %Store passing spots (overwrite older, which are always smaller)
    for ss = find(scorePass)'; enuc(S(ss).PixelIdxList) = true; end

end

%Erode further if a cytoplasmic segmentation
%   Gradients in cyto view tend to overestimate the nucleus
%   Erode depth scaled by minimum nuc diameter (prevent overerosion)
if op.seg.cyt  
    enuc = imerode(enuc, strel( 'disk', erode_grad ));    
end


%% Store estimated nuclei data
%Re-segment stored nuclei spot (re-orders spots)
S = regionprops(enuc, 'Centroid', 'Area', 'PixelIdxList');
nCtr = cat(1, S.Centroid);  nAre = cat(1, S.Area);  
im_op = im(:,:,op.seg.chan);
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
% %Pack (compress) nuclear mask image
% enuc = bwpack(enuc);

