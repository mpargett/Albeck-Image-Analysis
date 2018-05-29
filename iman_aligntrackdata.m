%IMAN_ALIGNTRACKDATA
%   Align tracked coordinates to segmented centroids, to allocate raw
%   masked data to tracks.  Used after both segmentation and tracking in
%   the Celltracer procedure.
%
%   [vc_out, vi, viv] = iman_aligntrackdata(coord, valcube)
%       selects the nearest centroid from valcube for each tracked
%       coordinate in coord.  Overlaps (centroids linked to two tracks) and
%       out-of-range centroids (greater than estimated cell radius) are
%       rejected.
%
%   Outputs:
%   vc_out  - valcube array (see celltracer)
%   vi      - raw indices aligned
%   viv     - index of valid indices, vi(viv) yields valid indices from raw
%               (input) valcube array



function [vc_out, vi, viv] = iman_aligntrackdata(coord, valcube)
%Version check provision
if strcmpi(coord,'version'); vc_out = 'v1.0'; return; end

%% Curate input coordinates
%Get scale of dataset
[ncells, ntime, ~] = size(coord);  nfield = size(valcube{1}, 3);

%Align tracked coordinates to nearest centroid for each time point
vc_out = nan(ncells, ntime, nfield);     %Initialize output dataset
%   Initialize coordinates as tracked estimates
vc_out(:,:,end-2:end-1) = coord;
for st = 1:size(coord,2)    %FOR each time point
    %   Reshape by pre-allocating is faster than squeeze
    vcc = nan(size(valcube{st},1),2);  
    %Get centroid coordinates from valcube
    vcc(:,:) =  valcube{st}(:,1,end-2:end-1);
    %Short circuit for empty fields
    if all(any(isnan(vcc),2));  continue; 	end
    %Get estimate radius of cell (from segmentation)
    erad = sqrt(valcube{st}(:,1,end)./pi); 
    
    %Get pairwise distance between tracked coordinate and centroids
    vdist = sqrt( bsxfun(@minus, coord(:,st,1), vcc(:,1)').^2 + ...
                  bsxfun(@minus, coord(:,st,2), vcc(:,2)').^2   );
    %   Exclude distances outside the radius of the nucleus
    vdist(bsxfun(@gt,vdist, erad')) = NaN;
    
    %Align coordinates to closest valcube centroids
    [vm, vi] = min(vdist, [], 2);
    %   Find any overlaps (centroids linked to two coords)
    ovr = bsxfun(@eq, vi, vi');
    ovr(1:(numel(vi)+1):end) = false;  ovr = any(ovr,2);
    %   Neglect entries for overlapped coordinates and failed location
    viv = ~ovr & ~isnan(vm);    %Map to valid indices
    
    %Rearrange valcube entries per tracked coordinates
    vc_out(viv,st,:) = valcube{st}(vi(viv),1,:);
end

end


%OBSOLETE - retained for debugging

% % coord = round(coord);   %Ensure integer coordinates
% %Enforce Minimum Values on coordinates (xy)
% coord( repmat(any(coord < 1, 3), [1,1,2]) ) = NaN;
% %Enforce Maximum Values on coordinates (xy)
% coord( repmat(any(bsxfun(@lt, coord, shiftdim(imsz,-2)), 3), [1,1,2]) ) = NaN;
% % coord( (coord(:,:,2) > sz(1) | coord(:,:,1) > sz(2)), :, : ) = NaN;

% %Convert coordinates to linear indices, from subscripts
% lind = sub2ind( sz, coord(:,:,2), coord(:,:,1) );
% gi = find(lind>0);      %Store index of 'good' coordinates
% %Get filtered linear index (nonzero value tracks)
% lind = lind( gi );      %Note: gi is retained to fill outputs properly