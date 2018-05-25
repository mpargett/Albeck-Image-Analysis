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
%       maxSmooth   -Maximum fraction of a cell's area that can be convex
%       sigthresh   -Signal value to use as a foreground threshold
%       hardsnr     -Logical - TRUE is sigthresh is a hard cutoff
%       nth         -Number of thresholds to use for segmentation
%


function [movieInfo, enuc] = iman_cellid(im, pin, bkg)
%Version check provision
if strcmpi(im,'version'); movieInfo = 'v2.1'; return; end

%% Operation parameters
%Default parameters
p = struct('chan', 1, 'cyt', false, 'minD', 8, 'maxD', 25, ...
           'maxEcc', 0.9, 'Extent', [0.65,0.9], 'maxSmooth', 0.5,...
           'mqw', [1,1,1], 'sigthresh', [], ...
           'hardsnr', false, 'nth', 20, 'nrode', 0, ...
           'localscale', false, 'usegrad', false);

%Apply provided parameters
for s = fieldnames(pin)'; 	p.(s{1}) = pin.(s{1});  end
p.maxSmooth = p.maxSmooth.*4;  %Adjust max smoothness for 1/4 scale
mqw = p.mqw(:)./sum(p.mqw);    %Scale quality weights
p.chan = p.chan(1);     %Use only main segmentation channel, if multiple

bkg_rat = 0.5;          %Min foreground to background ratio (~SNR-1)
%   Get image background intensity
if ~isempty(bkg); otim_bkg = bkg(p.chan); else otim_bkg = 0; end
%   Get background ratio for a threshold on fore- vs. back-ground
if ~isempty(p.sigthresh)
    %   Use a conservative estimate of 1/2 the SNR as threshold
    %   unless directed to use hard SNR threshold (op.seg.hardsnr)
    bkg_rat = max(0, (p.sigthresh(1)./otim_bkg - 1)*(1 + p.hardsnr)/2 );
end

%% Preliminary calculations
%Calculations used for cell sizing
maxNucArea = round( (pi*(p.maxD - p.nrode*2)^2)/4 );
minNucArea = round( (pi*(p.minD - p.nrode*2)^2)/4 );

%Define binary structuring elements
mflt = max(1, floor(p.maxD/40));
st0 = strel('disk', 1);  
st1 = strel('disk', mflt);  
st2 = strel('disk', p.nrode);
    
%   All further operations consider this original target image (otim)
otim = double( im(:,:,p.chan) );
%Set 'foreground' region based on background levels (previously removed)
%   Erode foreground mask to remove background adjacent pixels as well
%       (often salt noise in background regions)
otim_fore = imerode(otim > otim_bkg.*(bkg_rat), st1);   %Foreground mask
%Get background variance for metrics
if nnz(~otim_fore) > 3;  bkvar = var(otim(~otim_fore));  
else bkvar = 0; end


%% Image filtering (used for )
%Get gradient magnitude image
%   Filter out background by setting it to uniform level
gim = otim;     gim(~otim_fore) = otim_bkg.*(bkg_rat);
%   Create a predefined 2d filter. "Sobel edge-emphasizing filter"
hy = fspecial('sobel');   hx = hy';
%   Filter image for edges (gradients, via the Sobel filter)
gim = sqrt(imfilter(gim, hx, 'replicate').^2 ...  %Across x axis.
    + imfilter(gim, hy, 'replicate').^2)./8;     %Across y axis.

%Use Laplacian of Gaussian filter for thresholding
logf = fspecial('log', 5, 0.8);
%   Generates the Segmentation Target Image (stim)
stim = double(otim - imfilter(otim, logf, 'replicate'));

%Gradient in segmentation only if requested - 2018-01-18
%   Add image to gradients (penalizes 'high' regions for cyto)
if p.usegrad;   tstim = stim + gim;   
else            tstim = stim;  %Otherwise, just use image
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
    if p.localscale; stim = tstim./imdilate(stim,strel('disk',2*mflt)); end
end


%% Image thresholding search
%Define thresholds at linearly spaced quantiles
thresholds = prctile(stim(otim_fore), linspace(90, 5, p.nth) );
clear otim_fore;

%FOR each threshold, evalute appearance of nucleus-like shapes
enuc = inf(size(stim)); 
for s = 1:length(thresholds)
    %Threshold the gradient image, making a temporary image (tim)
    tim = stim > thresholds(s); %Binary, based on threshold value
    
    %Erode and dilate to eliminate salt noise in binary
    tim = bwlabel(imdilate(imerode(tim, st1),st1));
    %Dilate to remove holes AFTER labeling, to prevent region joining
    tim = noConMop(noConMop(tim, @imdilate, st1), @imerode, st1);
    %IF extra nuclear erosion is called, perform now
    if p.nrode > 0;  tim = imerode(tim, st2);   end
    
    %Get properties of each region in image (Area, Perimeter, ...)
    S = regionprops(tim, 'Area', 'Eccentricity', ...
        'Extent', 'Orientation', 'PixelList', 'PixelIdxList',...
        'MajorAxisLength', 'MinorAxisLength');
    nucArea  = cat(1,S.Area);           %Vector of Areas
    nucOr  = cat(1,S.Orientation).*pi./180; %Vector of Orientations  
    a = cat(1,S.MajorAxisLength);       %Vector of Major Axes^2
    b = cat(1,S.MinorAxisLength);       %Vector of Minor Axes^2
    
    %Perform primary shape filtering
    %   Calculate adjusted nuclear extent, based on orientation
    %       Solved by estimating size of bounding box for the rotated
    %       hypothetical ellipse, and dividing by unrotated area
    nucExt = cat(1,S.Extent) .* sqrt( 1 + ((a./b).^2 + (b./a).^2 - 2)...
        .*(sin(nucOr).*cos(nucOr)).^2 );
    
    %Determine objects passing size/shape fiters
    scorePass  = nucArea > minNucArea & nucArea < maxNucArea &...
        cat(1,S.Eccentricity) < p.maxEcc & ...
        nucExt > p.Extent(1) & nucExt < p.Extent(2);
    if ~any(scorePass); continue; end %Short circuit if all fail
    sp = find(scorePass)';  keeper = true(size(sp));
    
    %Evaluate mask quality metrics   
    %----------------------------------------------------------------------
    %Get outer boundary around each mask
    dtim = noConMop(tim, @imdilate, st0); %Expand masks by minimum amount
    dtim(logical(tim)) = 0; %Remove original mask to give boundaries
    %   Get properties for expanded boundaries
    S2 = regionprops(dtim, 'Area', 'PixelIdxList', 'PixelList', 'Centroid');
    
    %Boundary gradient: Maximize the gradient on the boundary
    %   Get average boundary gradient, scale by peak value and subtract
    %       from 1 to get cost. (Actually rescaled to map 5% of peak to 1
    %       "very bad" and 25% "good" to 0) 
    mq_grad = arrayfun(@(x,y)1 - 5*( mean(gim(x.PixelIdxList)) ...
                ./ max(otim(y.PixelIdxList)) - 0.05 ), S2(sp), S(sp));
        mq_grad(mq_grad < 0) = 0;  %Avoid cost on negative values
    
    %Smoothness: Minimize concave deviations on boundary
    %   Calculate the convex hull area per region (by removal)
    %       Subtract half of boundary to estimate internal area
        cha = arrayfun(@(x)sub_chull(x.PixelList, x.Centroid), S2(sp)) ...
            - cat(1,S2(sp).Area)/2;
    %   Scale smoothness cost by convex hull area
    %       Pre-weighted for 1/4 of area being concave to be "very bad"
    mq_smooth = 4*(cha - nucArea(sp))./cha;
     
    %Internal variance: Minimize CV (std/mean) of nuc
    %       Neglect typical variance of background (~dark noise)
    mq_cv = arrayfun(@(x)sqrt(var(otim(x.PixelIdxList)) - bkvar)...
         ./mean(otim(x.PixelIdxList)), S(sp));
    mq_cv = real(mq_cv);  %Imaginaries arise if var less than bkvar
    
    %Exclude regions based on quality metric bounds (secondary filter)
    keeper(mq_smooth > p.maxSmooth) = false;
    
    %Aggregate quality metrics (weighted sum of squares)
    %   Squaring weights improvement of a bad metric over a good one
    mqm = ([mq_grad, mq_smooth, mq_cv].^2)*mqw;
    %----------------------------------------------------------------------
       
    %Check mqm against past, and store passing spots with better mqm
    %   Can safely overwrite, new masks will always be larger
    for ss = 1:numel(sp);
        if keeper(ss) && mqm(ss) < min(enuc(S(sp(ss)).PixelIdxList))
            enuc(S(sp(ss)).PixelIdxList) = mqm(ss); 
        end
    end    
end

%Retain quality metric values per mask
qmask = enuc;
%Restore nuclear size from hole removal 
%   (does not restore full size if p.nrode > 0)
enuc = imdilate(bwlabel(~isinf(enuc)), st1);


%% Store estimated nuclei data
%Re-segment stored nuclei spot (re-orders spots)
S = regionprops(enuc, 'Centroid', 'Area', 'PixelIdxList');
nCtr = cat(1, S.Centroid);  nAre = cat(1, S.Area);  
im_op = im(:,:,p.chan);
nAmp = arrayfun(@(x)sum(im_op(x.PixelIdxList)), S);
nAmpStd = arrayfun(@(x)std(im_op(x.PixelIdxList)),S);

%Store coordinates IF any spots were found
if ~isempty(nCtr);
    movieInfo.xCoord = nCtr(:,1);
    movieInfo.yCoord = nCtr(:,2);
    movieInfo.amp = nAmp(:,1);
    movieInfo.ampstd = nAmpStd(:,1);
    movieInfo.are = nAre(:,1);
    movieInfo.qual = qmask(sub2ind(size(qmask), round(nCtr(:,2)), round(nCtr(:,1))));
else  %Store empty output if no spots found
    movieInfo = struct('xCoord',[],'yCoord',[],'amp',[], 'ampstd', [],...
        'are',[], 'qual', []);
end
% 
end


%% Subfunction to get convex hull area from boundary points
function cha = sub_chull(xy, ctr)
%  Get vectors from centroid to boundaries
bv = bsxfun(@minus, xy, ctr);
%   Sort by angle and get indices (counter-clockwise)
[~, bi] = sort(atan2(bv(:,2), bv(:,1)), 1, 'ascend');

%Recursively remove points in concavities to leave only the hull
cp = true;
while any(cp)   %WHILE any concavities were found on the last run
    %Get vectors along boundary points
    db = diff(xy([bi;bi(1)],:),1,1);
    %Use cross-product of adjacent vectors, which is negative if the second
    %   point is concave (when ordered counter-clockwise)
    cp = (db([end,1:end-1],1).*db(:,2) - db([end,1:end-1],2).*db(:,1)) < 0;
    %Remove concavities from boundary
    bi = bi(~cp);
end

%Calculate area of convex hull
db = xy([bi;bi(1)],:);  %Get xy coordinates of convex hull points
%   Calculate area of polygon (shoelace formula/Green's theorem)
cha = sum(db([end,1:end-1],1).*db(:,2) - db([end,1:end-1],2).*db(:,1))/2;
end


%% No conflict grayscale dilation/erosion
function im = noConMop(im, mop, stl)
%Prepare inverted image (for reverse operation)
iim = im;                   iim(iim == 0) = Inf; 
%Perform primary operation (dilation or erosion)
im = mop(im, stl);
%Repeat with inverted label values to observe overlaps
%   Inverted operation expands labels in the opposite direction
iim = -mop(-iim, stl);      iim = iim.*(iim ~= Inf);
%Remove the overlaps from both boundaries to avoid joining
im = im.*(iim == im);
end

%% Residual notes
%Other quality metrics considered
% %Size metric
% %   Compare size to average of bounds, and scale
% mq_size = (nucArea(sp) - avgNucArea)/maxNucArea;
% %Gradient Variation: Minimize the CV of the gradient on the boundary
% %   To penalize non-uniform boundaries (bright spots, etc.)
% mq_gcv = arrayfun(@(x)(std(gim(x.PixelIdxList))) ...
%     ./mean(gim(x.PixelIdxList)), S2(sp));
%
%
%     %TEMP testing
%     %spt = sp; sp = 1:1:numel(S); sp = spt;
%     m{s} = [mq_grad, mq_smooth, mq_cv];
%     imbah = im./prctile(im(:),99); [imR,imG,imB] = deal(imbah);
%     [nr,nc] = size(im);  aa = [-2,-1,0,1,2]; aa = [aa, aa+nr];
%     for sq = 1:numel(sp)
%         x = S2(sp(sq)).PixelIdxList;
%         imB(x) = 1; imR(x) = 0; imG(x) = 0;
%         xy = round(S2(sp(sq)).Centroid);
%         x = sub2ind(size(im),xy(2),xy(1)) + aa - 4*nr;  
%         if any(x<0 | x>(nr*(nc-8))); continue; end
%             imG(x) = mq_grad(sq);imB(x) = 0; imR(x) = 0;
%             x = x + 2*nr;
%             imG(x) = mq_smooth(sq);imB(x) = 0; imR(x) = 0;
%             x = x + 2*nr;
%             imG(x) = mq_cv(sq);imB(x) = 0; imR(x) = 0;
%     end
%     figure; imshow(cat(3,imR,imG,imB),[]); hold on;
%     xc = cat(1,S2.Centroid);
%     xca = bsxfun(@times, a./2, ([cos(-nucOr), sin(-nucOr)]));
%     xcb = bsxfun(@times, b./2, ([cos(pi/2-nucOr), sin(pi/2-nucOr)]));
%     plot([xc(:,1), xc(:,1)+xca(:,1)]',[xc(:,2), xc(:,2)+xca(:,2)]','r-');
%     plot([xc(:,1), xc(:,1)+xcb(:,1)]',[xc(:,2), xc(:,2)+xcb(:,2)]','r-');
%     for sq = 1:numel(sp)
%         text(xc(sp(sq),1), xc(sp(sq),2), num2str(mq_cv(sq)), 'Color', 'b');
%     end 
    