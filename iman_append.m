%IMAN_APPEND
%   Prepare secondary imaging data as an appendix to celltracer data.
%   Aligns the end frame of the first image data to the first frame of the
%   latter data.  
%
%   vc_app = iman_append(in1, in2)
%
%   in1 and in2 must both be processed imaging data, as one of:
%       valcube array [num]
%       filename of mat-file containing 'valcube' [char]
%       cell array of either of the above [cell]
%
%   vc_app, the output, is a cell array containing the aligned data (from
%       the 2nd input) such that each row is matched with the corresponding
%       row in the 1st input data set.  Rows that do not match are
%       discarded.  Cell array elements are ordered matching the input cell
%       arrays.
%
%   Version 1.0

%   M. Pargett, Albeck Lab, UC Davis, MCB.  2016

function vc_app = iman_append(in1, in2, varargin)
%Version check provision
if strcmpi(in1,'version'); vc_app = 'v1.0'; return; end

p.xyind = [];     %Default XY indices are empty (found as end-2:end-1)
p.trim  = {};     %Default is to not trim (empty)

%Input option pair parsing:
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end
%Splits pairs to a structure
for s = 1:2:nin;   p.(lower(varargin{s})) = varargin{s+1};   end

%% Get processed data
%   Check input when receiving cell arrays
if iscell(in1) || iscell(in2); assert(iscell(in1) && iscell(in2) ...
        && numel(in1) == numel(in2), 'IMAN:AppendCheck', ...
        ['When submitting multiple data sets for appending, they must ',...
        'be in cell arrays and the number of data sets in inputs 1 ',...
        'and 2 must match.']);
end

%Parse inputs
vc1 = sub_parse_data(in1);  vc2 = sub_parse_data(in2);

%Run registration and alignment
vc_app = cellfun(@(x,y)sub_align(x, y, p.xyind, p.trim), ...
    vc1, vc2, 'UniformOutput', false);

end

% --- Sub-Function: Parse Input Data --- 
function d = sub_parse_data(d)
%Extract valcube data from files, as needed
if ischar(d);        a = load(d);   d = {a.valcube};
elseif iscell(d)    
    if ischar(d{1}) %IF multiple names, load all, place in cell array
        a = cellfun(@(x)load(x,'valcube'), d);  d = {a.valcube};
    end
elseif isnumeric(d); d = {d};  %Ensure arrays are in a cell
end

end


% --- Sub-Function: Perform Registration and Alignment --- 
function vc_app = sub_align(vc1, vc2, xycs, t)
%% Register coordinates from valcube arrays
%   Bring coordinates to the same orientation (that of the former set)

%Determine data sizes
nv1 = size(vc1,3);
nv2 = size(vc2,3);    nt2 = size(vc2,2);
if isempty(xycs); xyc1 = nv1-2:nv1-1; xyc2 = nv2-2:nv2-1;
else  xyc1 = xycs{1}; xyc2 = xycs{2};
end

%Get coordinate sets for Scene (scn) and Model (mdl)
scn = squeeze(vc1(:, end, xyc1));    mdl = squeeze(vc2(:, 1  , xyc2));

%Trim Data as needed (may improve problem smoothness for large mismatches)
if ~isempty(t);
    scn = scn(  (scn(:,1) > t{1}(1) & scn(:,1) < t{1}(2)) & ...
                (scn(:,2) > t{1}(3) & scn(:,2) < t{1}(4)), : );
    mdl = mdl(  (mdl(:,1) > t{2}(1) & mdl(:,1) < t{2}(2)) & ...
                (mdl(:,2) > t{2}(3) & mdl(:,2) < t{2}(4)), : ); 
end

%Register Model coordinates to Scene
[mdl, tp, ctf] = iman_kcreg( scn, mdl );

%Re-place coordinates with those transformed
if isempty(t);   vc2(:, 1, xyc2) = mdl;
else   nc2 = size(vc2,1);  
    vc2(:, 1, xyc2) = ctf([squeeze(vc2(:, 1, xyc2)), ones(nc2,1)], tp);
end
%Replace time dimension in scene coordinates
scn = shiftdim( shiftdim( shiftdim(scn,1), -1), 2);


%% Align coordinates across gaps
[vc_app, vi, viv] = iman_aligntrackdata(scn, {vc2(:,1,:)});

if nt2 > 1  %IF latter data have more than one time point
    %Transform all latter coordinate sets
    for s = 2:nt2
        vc2(:,s,xyc2) = ctf([squeeze(vc2(:, s, xyc2)), ones(nc2,1)], tp); 
    end 
    %Align all latter data (neglects new traces in latter data)
    vc_app(:, 2:nt2, :) = vc2(vi(viv), 2:nt2, :);
end
end

