%IMAN_CELLMASK
%   Segment live cell images to generation nuclear and cytoplasmic masks,
%   and mean Cell Trace values over masked images.  Use this function after
%   a cell track has been established (e.g. via u-Track).
%
%FIXME make a real header
%FIXME scale EfA averaging with FRET concentration

%Updated 2.26.17 to improve indexing into image via masks for speed.

function [valcube, mask] = iman_cellmask(imin, m, op, nmasks)
%Version check provision
if strcmpi(imin,'version'); valcube = 'v2.1'; return; end

%% Operation parameters
ncgap  = floor((0.05*op.seg.maxD)) + 2;  %Size of gap between nuc mask and cyto
ncring = floor((0.05*op.seg.maxD)) + 1;  %Desired cyto mask thickness
ncexpand = ncgap + ncring;  %Size to expand nuclear mask for cyto
outpct = [20,80];                   %Percentiles to reject for outliers

%Identify input data size
if numel(imin) < 10;    sz = imin;  else    sz = size(imin);    end
if numel(sz) > 2; nchan = sz(3);  sz = sz(1:2); else nchan = 1; end	

%Check number of cross-channel ratios requested
nrt = numel(op.msk.rt);


%% IF coordinate/track information not provided, return valcube order
%   In this case, imin may contain only the size of the image data
if isempty(m)
    valcube = cell(2*(nchan+nrt), 1);
    %Get Channel names (or fill with defaults)
    if ~all(isfield(op, {'cname','cind'})) || ...
            min(numel(op.cname), numel(op.cind)) < nchan
        %Default channel names are Ch1, Ch2, Ch3...
        cns = cellfun(@(x)['Ch',num2str(x)], num2cell(1:nchan),...
            'UniformOutput', false);
    else    cns = op.cname(op.cind);
    end
    %Filling order is Nucs,Cytos,Ratios,Coords,Median
    lcn = {'_Nuc','_Cyto'};
    valcube(1:nchan) = cellfun(@(x)[x,lcn{1}], cns, ...
        'UniformOutput', false);
    valcube(nchan+(1:nchan)) = cellfun(@(x)[x,lcn{2}], cns, ...
        'UniformOutput', false);
    for sr = 1:nrt    %Ratios based on definitions above
        valcube{2*nchan + 2*sr-1} = [cns{op.msk.rt{sr}(1)},lcn{1},'/',...
            cns{op.msk.rt{sr}(2)},lcn{1}];
        valcube{2*nchan + 2*sr} = [cns{op.msk.rt{sr}(1)},lcn{2},'/',...
            cns{op.msk.rt{sr}(2)},lcn{2}];
    end
    %Standard appendices -      coordinates    nuc area
    valcube(end + (1:3)) = {'XCoord', 'YCoord', 'nArea'};
    valcube{end + 1} = ['Note: All values are background subtracted.  ',...
        'Cross-channel ratios are taken prior to averaging.'];
    
    %IF any FRET channels are unmixed, replace channel names in output set
    if op.unmix && ~isempty(op.msk.fret) && isstruct(op.msk.fret)
        fn = fieldnames(op.msk.fret);
        for s = 1:numel(fn)
            valcube = regexprep(valcube, op.msk.fret.(fn{s}){1}, ...
                [fn{s},'_efa'], 'ignorecase');
            valcube = regexprep(valcube, op.msk.fret.(fn{s}){2}, ...
                [fn{s},'_c'], 'ignorecase');
        end
    end
    mask = [];    return
end


%% Prepare and label masks
%Split out image channels (for indexing ease when masking)
imin = mat2cell( imin, sz(1), sz(2), ones(1,nchan) );

%Unpack and label mask image
nuclm = bwlabel(nmasks);

%Get tracked labels
lbl = nuclm( sub2ind(sz, round(m.yCoord), round(m.xCoord)) );
nnuc = numel(lbl);  %Number of segmented coordinates


%% Extend to cytoplasm masks
%Thicken nuclear region (to make a donut for cytoplasm)
%   Can simply dilate because nuclm is a label matrix (different values for
%       different regions, so none will merge)
expdisk = strel('disk', ncexpand);      %Expansion depth definition
%   Expand cytoplasmic mask from final nuclear mask
cytlm = imdilate(nuclm, expdisk);
%   Repeat with inverted label values to observe overlaps (regions where
%       two cytoplasm estimates collide)
ctemp = nuclm; ctemp(ctemp == 0) = Inf;  %Inverted dilation expands labels
ctemp = -imdilate(-ctemp, expdisk);     %   in the opposite direction
%   Remove the overlaps from both cyto estimates to preclude bias/conflict
ctemp(ctemp == Inf) = 0;  cytlm(~(ctemp==cytlm)) = 0; 

%Remove Nuclear pixels from Cytoplasm region
cytlm = cytlm.*~(imdilate(nuclm, strel('disk', ncgap))>0);

%Store mask matrices (binary for whole image masking)
mask.nuc = logical(nuclm);  mask.cyt = logical(cytlm);


%% Calculate and fill outputs (the 'valcube' matrix)
%Generate pre-average cross-channel ratio images (background subtracted)
rim = cell(nrt,1);  nvc = nchan*2;
for sr = 1:nrt
    %Calculate pre-average ratio image
    rim{sr} = imin{op.msk.rt{sr}(1)}./imin{op.msk.rt{sr}(2)};
    %   Exclude pixels with values below zero
    rim{sr}( imin{op.msk.rt{sr}(1)} < 0 | imin{op.msk.rt{sr}(2)} <= 0 ) = NaN;
end

%Store values (individual and ratios)
%   Initialize valcube
valcube = nan(nnuc, 1, nvc+2*nrt+3);

%Labels and images are expanded to vectors by masking of all elements
%   (i.e. nuclei of cytoplasms) at once, then per-label values are
%   aggregated and averaged (with tails trimmed).  This process improves
%   speed over per-label masking on the whole image.
%Collect label vectors
lv.nuc = nuclm(mask.nuc);  lv.cyt = cytlm(mask.cyt);
lvn = fieldnames(lv);       %Label names
%Sort label vectors and determine indices bounding each label
[lvs.nuc,lvsi.nuc] = sort(lv.nuc);  [lvs.cyt,lvsi.cyt] = sort(lv.cyt);
ind.nuc = [0;find(diff(lvs.nuc));length(lvs.nuc)];

%Inject empties into cytoplasm index list, as needed for proper matching
%   Identify missing cytoplasm indices (lost in mask expansion)
uni = unique(lvs.nuc);  uci = unique(lvs.cyt);
%   Determine intersect (valid cyto labels)
[~, ia] = intersect(uni,uci);
%   Determine difference (missing cyto labels)
injs = uni;  injs(ia) = [];
%   Built cyto index list
ind.cyt = [0;nan(nnuc,1)];  %Initialize with start index
%   Fill with valid indices, including end index of last point
ind.cyt(ia+1) = [find(diff(lvs.cyt));length(lvs.cyt)]; 
%   Inject missing indices (iterative for sequential gaps)
while ~isempty(injs)
    ind.cyt(injs+1) = ind.cyt(injs);
    injs = injs(isnan(ind.cyt(injs+1)));
end

lblidx = find(lbl)';        %Pre-define index of valid labels
for sc = 1:nchan + nrt      %Includes cross-channel ratios
    for sl = 1:numel(lvn)
        %Get image vector
        sr = sc - nchan;  isr = sr > 0;
        if isr; 	imv = rim{sr}(mask.(lvn{sl}));
        else       	imv = imin{sc}(mask.(lvn{sl}));    end
        imv = imv(lvsi.(lvn{sl}));  %Sort image vector by appropriate label
        
        %Determine output 'channel' / valcube 'slice'
        vcs = (~isr)*(sc + nchan*(sl-1)) + isr*(nvc + 2*(sr-1) + sl);
        
        %Process label-by-label values
        for s = lblidx
            %Collect labelled values
            vals = imv(ind.(lvn{sl})(lbl(s))+1 : ind.(lvn{sl})(lbl(s)+1));
            %Trim value distribution tails for robustness
            %   i.e. outlier removal by exclusion of a percent of extrema
            pctb = prctile(vals,outpct);
            %Take mean of core values
            valcube(s, 1, vcs) = mean(vals(vals>pctb(1) & vals<pctb(2)));
            %   Note:  Percentile outlier exclusion also removes NaNs. 
        end
    end
end

%Pack mask matrices for compression
mask.nuc = bwpack(mask.nuc);  mask.cyt = bwpack(mask.cyt);

%Store coordinate values, appended
valcube(:,:,nvc+2*nrt + (1:3)) = cat(3, m.xCoord, m.yCoord, m.are);


