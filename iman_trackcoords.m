%IMAN_TRACKCOORDS
%   Extract tracked coordinates from u-Track output ('tracksFinal')
%
%   [c, wmm] = iman_trackcoords(tF)
%       returns a nTracks x nTime x 2 (X,Y) matrix of coordinates (c) and
%   (optionally) a nTracks x 2 index of lineage (wmm).  tF must be the
%   'tracksFinal' output from u-Track.
%
%   Optional lineage index wmm contains the index of the track from which
%   the current track split (1st column) or that which it merged into (2nd
%   column).  If the track initiated without splitting (or ended without
%   merging), the correspoding value is NaN.


function [c, whosmymommy] = iman_trackcoords(tF)
%Version check provision
if strcmpi(tF,'version'); c = 'v2.0'; return; end

%% Get track dimensions
nT = max( arrayfun(@(x)x.seqOfEvents(end,1), tF) ); % # of Times
nct = length( tF );     %Number of compound tracks
ntk = arrayfun(@(x)size(x.seqOfEvents,1)./2, tF);   %# of Cells / compound
nC = sum( ntk );        % # of Cells

c = nan(nC, nT, 2);     %Initialize Coordinate matrix
whosmymommy = nan(nC, 2);      %Initialize Lineage matrix
%% Extract each tracks coordinates and lineage info
for s = 1:nct
    %Note: seqOfEvents gives: 1-frame index of event, 2-(1=start, 2=end),
    %   3-local track index, 
    %   4-(NaN=true start/end, value=local index of interacting track)
    
    %Get start and end indices of compound track (global indices)
    cst = min(tF(s).seqOfEvents(:,1));
    cnd = max(tF(s).seqOfEvents(:,1));
    
    %For each component track in this compound
    pstk = sum(ntk(1:s-1));     %Number of tracks already done
    for st = 1:ntk(s)
        %Get indices for this track
        tk = pstk + st;     %Global track index
        
        %Fill Coordinates
        c(tk, cst:cnd, 1) = tF(s).tracksCoordAmpCG(st, 1:8:end); %x-Coord
        c(tk, cst:cnd, 2) = tF(s).tracksCoordAmpCG(st, 2:8:end); %y-Coord

        %Record lineages (global index values)
        for se = 1:2    %1:Start (Split from), 2:End (Merged into)
        whosmymommy(tk, se)	= pstk + tF(s).seqOfEvents( ...
            tF(s).seqOfEvents(:,3)==st & tF(s).seqOfEvents(:,2)==se, 4 );
        end
    end
end

%NaNify zero values coordinates
c(c==0) = NaN;

%Sort tracks, first by start time, then by length
cg = ~isnan(c(:,:,1));  %Logical matrix for good indices
%   Get indices for starts (pattern = FALSE, TRUE), working along tracks
[tst, ci] = find( (cg & ~[false(nC,1), cg(:,1:end-1)])' );
%   Get first start index per cell track (ignore starts after gaps)
[~,st1,~] = unique(ci);  tst = tst(st1);
%   Get length of tracks
tlength = sum(cg, 2);   tlength(tlength == 0) = NaN;
%   Ensure alignment of starts if any tracks are empty
tstart = nan(nC,1);     tstart(tlength > 0) = tst;

%Sort lineage indices to match
[~,sli] = sort(tlength, 'descend');     %First on length (becomes secondary)
[~,sti] = sort(tstart(sli), 'ascend'); 	%Last on start time (primary)
%   Apply final sort order
si = sli(sti);  c = c(si,:,:);  whosmymommy = whosmymommy(si,:);
%   Invert sorting order for lineage mapping
[~,isi] = sort(si);
%   Apply sorting to lineage values
for s = unique(whosmymommy(~isnan(whosmymommy)))'	
    whosmymommy(whosmymommy == s) = isi(s);  
end


%% Interpolate gaps in coordinate vectors
%   Fill bounded NaNs with interpolated values (NaNs outside stay NaN)
index = 1:nT;
for s = 1:nC
    %Get failed coordinates
    nanD = any( isnan(c(s,:,:)), 3 );
    if ~any(nanD); continue; end    %Skip if no failed coordinates
    %Skip any track with < 2 valid coordinates
    if nnz(~nanD) < 2; c(s, ~nanD, :) = NaN;  continue;   end
    %Interpolate between valid coordinates
    c(s, nanD, :) = interp1(index(~nanD), shiftdim(c(s, ~nanD, :),1), ...
        index(nanD), 'linear', NaN); 
end


end