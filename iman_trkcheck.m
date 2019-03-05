%IMAN_TRKCHECK
%   Check tracking behavior with uTrack and current parameters.  This
%   function will show the search regions implemented to link points
%   together across frames, to close gaps and to identify merge/split
%   conditions.
%
%   d = iman_trkcheck(ip,op,d,...)
%       Calls segmentation (via iman_segcheck) and uTrack procedures to
%   track cells based on the parameters provided in ip and op, which should
%   be identical to those used for an iman_celltracer call.
%
%   In the plots produces, the following conventions are used:
%       * for current position of starting/ending track
%       -- for its path in time
%       x for current position of other tracks
%       - for their paths in time
%       o for position of tracks ending/starting in past/future, color
%           coded (parula), for nearness to the current time
%       [Large circles/rectangles] for search regions.  Linking regions are
%           shown as magenta color coded, gap-closing regions are color
%           coded as with the o markers above.
%       [Yellow lines (-)] for merge/split points, linking the two tracks.
%
%   The output, and optional input, d, is a data structure containing the
%   processed data, allowing for rapidly changing plotting parameters and
%   regenerating figures.  If the ip or op parameters have been changed,
%   provide an empty d (d = []). 
%
%   Plotting parameters may be provided as Name/Value pairs after the 3rd
%   input (d).  They include:
%   
%   xy      -   "XY" index to process and plot.
%   t       -   Time indices to process (must be long enough to track)
%   starts  -   Logical, set TRUE to plot gap closing for new track starts.
%                   Set FALSE to plot for track ends (merges).
%   npth    -   Number of frames of tracked path to plot.
%   plotlinks - Logical, set TRUE to show the linking regions plot.
%   plotgaps -  Logical, set TRUE to show the gap-closing regions plot.
%   frm     -   Frame(s) to plot, must be within the range of p.t. May also
%                   indicate 'splits' or 'merges' to show frames before and
%                   after each event.
%   imthresh -  Threshold percentage of image intensity to plot as white.
%


function [d] = iman_trkcheck(ip, op, d, varargin)
%% Parse Inputs
p.xy  = 1;      %Default xy to use
p.t   = 1:10;   %Default time range
p.starts = true;    %Show gap closing from perspective of starts (not ends)
p.npth = 5;     %Number of frames of track path to plot around current
p.plotlinks = true;
p.plotgaps = true;
p.printid = true;
p.frm = 5;
p.imthresh = 99;    %Threshold percentage of image intensity for plotting

%Input option pair parsing:
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end
%Splits pairs to a structure
for s = 1:2:nin;   p.(lower(varargin{s})) = varargin{s+1};   end

%Frame selection overrides p.starts
if ~isnumeric(p.frm) && ischar(p.frm)
    switch lower(p.frm)
        case 'splits';  p.starts = true;  p.frm = NaN;
        case 'merges';  p.starts = false; p.frm = NaN;
        otherwise; error('Unknown parameter in p.frm');
    end
end


%% Run segmentation for requested set of frames
if ~exist('d','var') || isempty(d)
    pst = [];
    for st = 1:numel(p.t)
        [pst, d.imo(st)] = iman_segcheck(ip, op, 't',p.t(st), 'xy',1, ...
            'pastinfo',pst, 'display', false);
    end
    d.GMD = pst.GMD; %Retain GMD for any short circuits
end

%   Re-run uTrack if tF field is not populated
if ~isfield(d, 'tF') || isempty(d.tF)
    minfo = [d.imo.m];
    %Prepare for uTrack - movieInfo fields must include a 2nd column,
    %   for Std Deviation of positions and amplitude
    mif = intersect({'xCoord','yCoord','zCoord','amp'},...
        fieldnames(minfo));  %Relevant movieInfo Fieldnames
    for sf = 1:numel(mif);  %FOR each relevant name
        if size(minfo(1).(mif{sf}),2) < 2 %Skip if already filled
            %Check if coordinate or amplitude and calculate estimated
            %   St. Dev. if associated info is present
            if ~strcmpi('amp',mif{sf}) && isfield(minfo, 'are')
                mifun = @(x)[x.(mif{sf}), sqrt(x.are./pi)/6];
            elseif strcmpi('amp',mif{sf}) && ...
                    all(isfield(minfo, {'ampstd','are'}))
                mifun = @(x)[x.amp, x.ampstd./sqrt(x.are)];
            else  %Or just append zeros for the St. Dev.
                mifun = @(x)[x.(mif{sf}), zeros(size(x.(mif{sf})))];
            end
        end
        mtmp = arrayfun(mifun, minfo, 'Un', 0); %Append St Dev
        [minfo.(mif{sf})] = deal(mtmp{:});      %Replace in mI
    end
    
    %Pre-process radii, if running uTrack
    pxscl = (d.GMD.cam.PixSizeX + d.GMD.cam.PixSizeY)./2;
    op.trk.linkrad = op.trk.linkrad./pxscl;
    op.trk.gaprad = op.trk.gaprad./pxscl;
    %Run separate, but equivalent uTrack procedures to extract search areas
    [d.tF, d.ki, d.ptrk] = iman_utrack_call(minfo, op.trk);
end

%Get coordinates for points involved in splits or merges
smxy = cell(size(d.imo));
for s = find(arrayfun(@(x)size(x.seqOfEvents,1) > 2, d.tF))';
    %   Find any splits/merges in this compound track
    sm = ~isnan(d.tF(s).seqOfEvents(:,end)) & ...
        d.tF(s).seqOfEvents(:,2) == 2-p.starts;
    if any(sm) %Get current xy and connecting xy for any found
        cfrms = d.tF(s).seqOfEvents(sm,1);
        %   Get xys for starting/ending track in this frame
        ctrk = d.tF(s).seqOfEvents(sm,3:4);
        cst = min(d.tF(s).seqOfEvents(:,1));
        for st = 1:numel(cfrms)
            %Get xys for tracks for this frame
            xys = d.tF(s).tracksCoordAmpCG(...
                ctrk(st,:),(cfrms(st)-cst)*8 + [1,2]);
            %Store
            smxy{cfrms(st)} = cat(1,smxy{cfrms(st)}, reshape(xys,1,4));
        end
    end
end

if all(isnan(p.frm))
    p.frm = find(~cellfun(@isempty, smxy));
    p.frm = unique( bsxfun(@plus, p.frm, [0;1] - p.starts) );
end

%Get final track coordinates if IDs are to be printed
[c, linfo] = iman_trackcoords(d.tF);
ids = cellfun(@num2str, num2cell(1:size(c,1)), 'Un', 0);


%% Plot each frame of interest with search areas desired
subp = p.plotlinks && p.plotgaps;
pp.starts = p.starts;  pp.npth = p.npth;  pp.imthresh = p.imthresh;
for s = 1:numel(p.frm)
    
    %Generate figure and axes
    fh(s) = figure;  set(fh(s), 'Position', [200, 300, 1400, 600]); %#ok<AGROW>
    if ~subp;  ah = axes(fh(s));  end %#ok<LAXES>
    
    %   Plot linking search regions
    if subp; ah = subplot(1,2,1); end
    if p.plotlinks
        plot_linking(ah, d.imo(p.frm(s)).im(:,:,op.seg.chan(1)), ...
            d.ki((p.frm(s)-1):p.frm(s)), p.frm(s), pp);
    end
    
    %   Plot gap closing search regions
    if subp; ah = subplot(1,2,2); end
    if p.plotgaps
        plot_gapclosing(ah, d.imo(p.frm(s)).im(:,:,op.seg.chan(1)), d.ptrk, ...
            smxy{p.frm(s)}, p.frm(s), pp);
    end
    
    %   Print ID numbers per cell
    if p.printid
        text(ah, c(:,p.frm(s),1), c(:,p.frm(s),2), ids, 'Color', 'y');
    end
end


end



%--------------------------------------------------------------------------
%Subfunctions
%--------------------------------------------------------------------------
%% Function to plot linking search regions
function [] = plot_linking(ah, im, kin, frm, pin)
p.imthresh = 99;

%Apply provided parameters
for s = fieldnames(pin)'; 	p.(s{1}) = pin.(s{1});  end

%Define generic circular plotting
x = cos(-pi:pi/8:pi);  y = sin(-pi:pi/8:pi);  

%Initialize figure and Plot current image
if isempty(ah); figure; ah = axes; end
    imshow(im, prctile(im(:),[0,p.imthresh]), 'Parent', ah);  hold on; 
title(['Linking Search Regions, Frame ',num2str(frm)]);
%Plot projected search regions (from past frame)
for s = 1:size(kin(1).sr,1)
    plot(ah, kin(1).stateVec(s,1) + x.*kin(1).sr(s), ...
        kin(1).stateVec(s,2) + y.*kin(1).sr(s), 'm-');
end

%Plot current feature centroids (from current frame)
scatter(ah, kin(2).stateVec(:,1), kin(2).stateVec(:,2), 'yx');
end


%% Function to plot gap closing search regions
function [] = plot_gapclosing(ah, im, pt, smxy, frm, pin)
p.starts = false;   %Search from starts (not ends)
p.npth = 5;         %Number of frames for which to plot paths
p.imthresh = 99;

%Apply provided parameters
for s = fieldnames(pin)'; 	p.(s{1}) = pin.(s{1});  end

%Set up positions to search from starts or ends
if p.starts;   	psch = pt.svo.poss;  ptgt = pt.svo.pose; %Positions
                lv = pt.svo.lvs;  sv = pt.svo.svs;       %Search vectors
else            psch = pt.svo.pose;  ptgt = pt.svo.poss; %Positions
                lv = pt.svo.lve;  sv = pt.svo.sve;       %Search vectors
end

%Get colormap
figure; cmap = colormap('parula'); close(gcf);

%Initialize figure
if isempty(ah); figure; ah = axes; end
    imshow(im, prctile(im(:),[0,p.imthresh]), 'Parent', ah); hold on; 
title(['Gap Closing Search Regions, Frame ',num2str(frm)]);

%Plot paths for tracks starting/ending in this frame
te = find(psch(:,3) == frm);    %Tracks starting/ending indices
nf = numel(te);                   %Number of features
%   Get only features starting/ending in this frame
fts = psch(te, :);
%Plot paths
pt.svo.trkhist  = full(pt.svo.trkhist);       %Convert to full matrix
pt.svo.trkhist(pt.svo.trkhist == 0) = NaN;    %NaNify zeros
%   Plot all, same color as current position
plot(ah, pt.svo.trkhist(te,1:2:end)', pt.svo.trkhist(te,2:2:end)', 'b-');    

%Plot select search regions for frame in window
%   Subsample frames in window
ffrm = floor(linspace(1,pt.gapwin,3));
%   Sample colormap to match
cc = interp1(linspace(0,ffrm(end),size(cmap,1))', cmap, ffrm');
%   Set up line types
lts = {'-','--'};
for s = 1:numel(ffrm)
    %Plot search boxes/disks for each feature
    for sf = 1:nf
        %   Get box definitions (centered on origin)
        [bb,isr] = get_searchregion(lv(:,ffrm(s),te(sf)), ...
            sv(:,ffrm(s),te(sf)), pt.svo.type(te(sf)),...
            pt.svo.maxDisp); 
        %Plot box (translated to feature position)
        plot(ah, fts(sf,1) + bb(1,:), fts(sf,2) + bb(2,:), ...
            lts{1+isr}, 'Color', cc(s,:));
    end
end

%Plot merges/splits
if ~isempty(smxy)
    %Plot line to split/merge partner
    plot(ah, smxy(:,1:2)', smxy(:,3:4)', 'y-', 'LineWidth', 2);
end

%Plot centroids for tracks ending before/starting after this frame
%   Resample colormap for all frames
cc = interp1(linspace(0,ffrm(end),size(cmap,1))', cmap, (1:pt.gapwin)');
for s = 1:pt.gapwin     %Consider all frames in window
    ts = ptgt(:,3) == frm + s*(1-2*p.starts);    %Ending/starting indices
    scatter(ah, ptgt(ts,1), ptgt(ts,2), [], cc(s,:), 'o');  %Plot
end

%Plot tracks that exist in the frame, with short paths
itrk = ~isnan(pt.svo.trkhist(:,2*frm));
%   Get time range to plot path
tr = max(1, frm - p.npth) : min(size(pt.svo.trkhist,2)/2, frm + p.npth);
trix = 2*(tr-1)+1; triy = 2*tr; %Get x and y coordinate indices
plot(ah, pt.svo.trkhist(itrk,trix)', pt.svo.trkhist(itrk,triy)', ...
    '-', 'Color', cc(1,:));
plot(ah, pt.svo.trkhist(itrk,2*frm-1)', pt.svo.trkhist(itrk,2*frm)', ...
    'x', 'Color', cc(1,:));

%Plot current positions (Last to ensure on top)
scatter(ah, fts(:,1), fts(:,2), 'b*');  

end


%% Function: Define plotting for search regions, from vectors
%   This to be added to the (forward time) centroid
function [bb,isr] = get_searchregion(vl,vs,islin,maxd)
isr = false;
%Restrict all vectors to max displacement
r = abs(vl) > maxd; if any(r); isr = true;  vl(r) = sign(vl(r))*maxd; end
r = abs(vs) > maxd; if any(r); isr = true;  vs(r) = sign(vs(r))*maxd; end
%Define plot points to outline the search region,
%   depending on the motion type (linear/brownian)
if isnan(islin) || ~islin
    %IF a brownian motion point, define bounding disk
    bb = [cos(-pi:pi/8:pi); sin(-pi:pi/8:pi)].*norm(vl);
else     %IF a linear motion point, define bounding box
    bb = [vl+vs, vl-vs, -vl-vs, vs-vl, vl+vs];
end
end



