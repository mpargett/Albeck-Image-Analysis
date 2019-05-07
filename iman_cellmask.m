%IMAN_CELLMASK
%   Segment live cell images to generation nuclear and cytoplasmic masks,
%   and mean Cell Trace values over masked images.  Use this function after
%   a cell track has been established (e.g. via u-Track).
%
%   [valcube, mask] = iman_cellmask(imin, m, op, nmasks)
%       returns valcube (nCells x nTime x nChannels matrix of values) and
%   mask (a structure with the binary masks used - packed via bwpack).  The
%   procedure must be provided the input image (imin), a coordinate
%   structure (m - same as uTrack's movieInfo), operation parameters (op)
%   and nuclear masks from which to start (nmasks, as a logical matrix).
%
%   If called without the coordinate info (m), iman_cellmask will return
%   the naming structure (as a cell array) of the valcubes that it will
%   generate given the corresponding image size and operation parameters.
%
%   The operation parameters structure, op, may contain fields:
%       cname       -List of channel names (cell array)
%       cind        -List of channel indices used
%       unmix       -Logical, TRUE is using spectral unmixing
%       seg.maxD    -Maximum nucleus diameter allowed from Segmentation
%       msk.rt      -Indices of channel to divide prior to averaging,
%                       as 1 x nRatios cell array with 2 indices per ratio
%       msk.fret    -Definition of any FRET channels for spectral unmixing,
%                       as a structure, with fieldnames declaring new FRET
%                       channels, and values being cell array with the
%                       channel names of the donor and emitter
%       msk.aggfun  -Definition of additional aggregate functions to use,
%                       as a structure with fields:
%               name    -Name of new function (prepends channel names)
%               chan    -Channel indices on which to use
%               loc     -Localization(s) on which to use (1: Nuc, 2:Cyt)
%               fun     -Function handle to be applied
%       msk.cytadj  -Adjustments to cytoplasmic masking distances, positive
%                       values increase 'ring' radius and width, negative
%                       values decrease (pixels).  Default: 0
%

%FIXME scale EfA averaging with FRET concentration

%Updated 2.26.17 to improve indexing into image via masks for speed.
%Updated 9.15.17 to include alternate aggregation functions. MP
%Updated 4.19.18 with cyto mask refinement. MP
%Updated 4.08.19 to pass shape metrics. MP

function [valcube, mask] = iman_cellmask(imin, m, op, nuclm, bkg)
%Version check provision
if strcmpi(imin,'version'); valcube = 'v2.4'; return; end

%% Operation parameters
if ~exist('bkg','var'); bkg = []; end
ncgap  = round((0.05*op.seg.maxD)) + 2;  %Size of gap between nuc mask and cyto
ncring = round((0.05*op.seg.maxD)) + 1;  %Desired cyto mask thickness

%Backward compatibility for scripts missing parameters
bc = {'nfilt',false; 'cfilt',false; 'nrode',0; 'cgap',0; 'cwidth',0; ...
    'appendshape',false};
for s = 1:size(bc,1)
    if ~isfield(op.msk, bc{s,1}); op.msk.(bc{s,1}) = bc{s,2}; end
end

%IF masking adjustments provided, apply (allowing backward compatibility)
%   Nuclear Mask re-dilation from segmentation erosion
if isfield(op.seg, 'nrode'); nrode = op.seg.nrode; else nrode = 0; end
nrode = nrode - op.msk.nrode;       %Nuclear Mask erosion
ncgap = ncgap + op.msk.cgap;        %Nuc-Cyto Mask gap
ncring = ncring + op.msk.cwidth;    %Cytoplasm Mask ring width

ncexpand = ncgap + ncring;  %Size to expand nuclear mask for cyto
outpct = [20,80];          	%Percentiles to reject for outliers

%Identify input data size
if numel(imin) < 10;    sz = imin;  else    sz = size(imin);    end
if numel(sz) > 2; nchan = sz(3);  sz = sz(1:2); else nchan = 1; end	

%Check number of cross-channel ratios requested
nrt = numel(op.msk.rt);
%Check for additional aggregation functions
if isfield(op.msk,'aggfun')
    nagg = numel(op.msk.aggfun); 
    nagt = sum(arrayfun(@(x)numel(x.chan)*numel(x.loc), op.msk.aggfun));
else nagg = 0; nagt = 0;
end


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
    %Filling order is Nucs,Cyts,Ratios,Coords,NuclearArea
    lcn = {'_Nuc','_Cyt'};
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
    %Custom aggregations (if requested)
    if nagg > 0
        ca = cell(numel(op.msk.aggfun),1);
        for s = 1:numel(op.msk.aggfun)
            %   Get numbers of channel and localizations involved
            nc = numel(op.msk.aggfun(s).chan);   %1 - nchan
            nl = numel(op.msk.aggfun(s).loc);    %1 - 2
            %   Expand name components and align to form combinations
            ca{s} = [repmat({[op.msk.aggfun.name,'_']}, nc*nl, 1), ...
                reshape(repmat(reshape(op.cname(op.msk.aggfun.chan),1,nc),...
                nl,1),nc*nl,1), repmat(lcn(op.msk.aggfun.loc)', nc, 1)];
            %   Pack into final names
            ca{s} = cellfun(@(x)[x{:}], num2cell(ca{s}, 2), 'Un', 0);
        end
        %   Append to output (vcorder, called valcube here)
        valcube = cat(1, valcube, ca{:});
    end
    %Standard appendices -      coordinates    nuc area
    valcube(end + (1:3)) = {'XCoord', 'YCoord', 'nArea'};
    %Shape metric appendices
    if op.msk.appendshape
    valcube(end + (1:5)) = {'nEccentricity',...
        'nOrientation', 'nExtent', 'nSolidity', 'nCV'};
    end
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
imin = num2cell( imin, [1,2] );

%If nuclear masks were eroded in processing, re-dilate, accomodating
%   desired extra erosion in final masks (op.msk.nrode)
%   (Performed here, on label matrices, preserves cell identity)
if      nrode > 0; 	nuclm = imdilate(nuclm, strel('disk', nrode)); 
elseif  nrode < 0; 	nuclm = imerode(nuclm, strel('disk', -nrode)); 
end


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


%% Clean up cytoplasm masks, based on segmentation channel(s) 
%Determine available segmentation channels
if numel(op.seg.chan) < 2; op.seg.chan = [op.seg.chan, NaN]; end
nuchan = op.seg.chan(1+op.seg.cyt);  %Get Nuclear channel
cchan = op.seg.chan(2-op.seg.cyt);  %Get Cytoplasm channel
%Clean Nuclear signal from cyto mask, if channel is available
if op.msk.nfilt && ~isnan(nuchan)
    nst = prctile(imin{nuchan}(logical(nuclm)), 5); %Nuc signal threshold
    cytlm(imin{nuchan} > nst) = 0; %Remove ~Nuclear signal from Cyto mask
end
%Enforce foreground levels of Cytoplasm signal in mask, if available
if op.msk.cfilt && ~isnan(cchan) && ~isempty(bkg);  
    bkg_rat = 0.5; %Default SNR - 1 is 0.5 (1.5 SNR)
    %Handle sigthresh as with cellid, if provided
    if ~isempty(op.seg.sigthresh) ... % But check for 2nd term if ~cyt
            && (op.seg.cyt || numel(op.seg.sigthresh > 1))
        bkg_rat = max(0, (op.seg.sigthresh(2-op.seg.cyt)./bkg(cchan)-1)...
            *(1+op.seg.hardsnr)/2 );
    end
    cytlm(imin{cchan} <= bkg(cchan).*bkg_rat) = 0; %Remove bkg pixels
end

%Store mask matrices (binary for whole image masking)
mask.nuc = logical(nuclm);  mask.cyt = logical(cytlm);


%% Prepare label/mask indices to store channel values
%Labels and images are expanded to vectors by masking of all elements
%   (i.e. nuclei of cytoplasms) at once, then per-label values are
%   aggregated and averaged (with tails trimmed).  This process improves
%   speed over per-label masking on the whole image.
%Collect label vectors
lv.nuc = nuclm(mask.nuc);  lv.cyt = cytlm(mask.cyt);
lvn = fieldnames(lv);  nlv = numel(lvn);       %Label names
%Sort label vectors and determine indices bounding each label
[lvs.nuc,lvsi.nuc] = sort(lv.nuc);  [lvs.cyt,lvsi.cyt] = sort(lv.cyt);

%Assemble bounding indices for label vectors
for sl = 1:nlv  %For both Nuc and Cyt indices
    %   Get bounding indices and number of labels incremented
    dd = diff([1;lvs.(lvn{sl})]); id = find(dd);  dd = dd(id);
    %   Assign bounding with duplication where labels skipped
    ind.(lvn{sl}) = [0; cell2mat(arrayfun(@(x,y)repmat(y,x,1), ...
        dd, id-1, 'Un', 0)); length(lvs.(lvn{sl}))];
end
%   Append ends for any lost cyto masks of highest label number
%       Only if necessary (robust to no masks)
if length(ind.cyt) < length(ind.nuc)
    ind.cyt(end+1:lvs.nuc(end)+1) = length(lvs.cyt);
end

%% Calculate and fill outputs (the 'valcube' matrix)
%Generate pre-average cross-channel ratio images (background subtracted)
rim = cell(nrt,1);
for sr = 1:nrt
    %Calculate pre-average ratio image
    rim{sr} = imin{op.msk.rt{sr}(1)}./imin{op.msk.rt{sr}(2)};
    %   Exclude pixels with values below zero
    rim{sr}( imin{op.msk.rt{sr}(1)} < 0 | imin{op.msk.rt{sr}(2)} <= 0 ) = NaN;
end

%Get labels (label value is 1st dim index of m, per iman_cellid)
lbl = unique(nuclm(:)); lbl = lbl(lbl>0);
nnuc = numel(lbl);  %Number of segmented coordinates

%   Initialize valcube
valcube = nan(nnuc, 1, nlv*(nchan+nrt) + nagt + 2 + 6*op.msk.appendshape);
lblidx = find(lbl)';        %Pre-define index of valid labels
%Fill standard channels
for sc = 1:nchan
    for sl = 1:nlv
        %Get image vector, and Sort by appropriate label
        imv = imin{sc}(mask.(lvn{sl})); imv = imv(lvsi.(lvn{sl}));
        %Get current valcube slice (channel) index
        vcs = sc+nchan*(sl-1);
        %Fill valcube slice (channel)
        valcube(:, 1, sc+nchan*(sl-1)) = fillchan(valcube(:, 1, vcs), ...
            imv, lbl, lblidx, ind.(lvn{sl}), outpct, @mean);        
    end
end

%Fill ratio channels
for sc = 1:nrt
    for sl = 1:nlv
        %Get image vector, and Sort by appropriate label
        imv = rim{sc}(mask.(lvn{sl})); imv = imv(lvsi.(lvn{sl}));
        %Get current valcube slice (channel) index
        vcs = nchan*nlv + sl+2*(sc-1);
        %Fill valcube slice (channel)
        valcube(:, 1, vcs) = fillchan(valcube(:, 1, vcs), ...
            imv, lbl, lblidx, ind.(lvn{sl}), outpct, @mean);        
    end
end

%Fill additional aggregate function channels
vcs = nlv*(nchan+nrt);
for sa = 1:nagg
    for sc = 1:numel(op.msk.aggfun(sa).chan)
        for sl = numel(op.msk.aggfun(sa).loc)
            %Get image vector, and Sort by appropriate label
            imv = imin{ op.msk.aggfun(sa).chan(sc) }(mask.(lvn{sl})); 
            imv = imv(lvsi.(lvn{ op.msk.aggfun(sa).loc(sl) }));
            %Get current valcube slice (channel) index
            vcs = vcs + 1;
            %Fill valcube slice (channel)
            valcube(:, 1, vcs) = fillchan(valcube(:, 1, vcs), imv, ...
                lbl, lblidx, ind.(lvn{ op.msk.aggfun(sa).loc(sl) }),...
                outpct, op.msk.aggfun(sa).fun);
        end
    end
end

%Store coordinate values and nuclear metrics, appended
valcube(:,:, vcs + (1:3)) = cat(3, m.xCoord, m.yCoord, m.are);
%Append shape metrics if called
if op.msk.appendshape
    valcube(:,:, vcs + (4:8)) = cat(3, m.nEcc, m.nOrient, ...
        m.nExt, m.nSold, m.nCV);
end

%Pack mask matrices for compression
mask.nuc = bwpack(mask.nuc);  mask.cyt = bwpack(mask.cyt);


end


%% Subfunction: Fill an output channel
function vc = fillchan(vc, imv, lbl, lblidx, lidx, outpct, aggfun)
%Process label-by-label values
for s = lblidx
    %Collect labelled values
    vals = imv(lidx(lbl(s))+1 : lidx(lbl(s)+1));
    %Trim value distribution tails for robustness
    %   i.e. outlier removal by exclusion of a percent of extrema
    pctb = prctile(vals,outpct); 
    %Take mean of core values
    vc(s) = aggfun(vals(vals>pctb(1) & vals<pctb(2)));
    %   Note:  Percentile outlier exclusion also removes NaNs.
end
end


