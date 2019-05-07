%IMAN_APPEND
%   Prepare secondary imaging data as an appendix to celltracer data.
%   Aligns cell traces from the end frame of the first image data to those
%   in the first frame of the second dataset.  
%
%   vc_app = iman_append(IN1, IN2, ...)
%
%   The input datasets, IN1 and IN2, must both be processed imaging data, 
%       as either a valcube array [numeric], or a filename of a mat-file
%       containing 'valcube' [char].
%
%   The output, vc_app, is a numeric 'valcube' array containing the aligned
%       data from IN2, such that each row is matched with the corresponding
%       row in IN1.  Rows from IN2 that do not match anything in IN1 are
%       discarded.   
%
%   Additional parameters may be passed as Name/Value pairs:
%   xyind - Cell array of indices of the X and Y coordinates in the input
%               valcube data (the "slices" of the 3rd dimension), with a
%               separate cell array element for input 1 and 2, such as
%               {[Xind1,Yind1], [Xind2,Yind2]}.  Typically not needed if
%               providing file names of processed data to load.     
%   trim  - Cell array of coordinates defining rectangular regions of the
%               frames to trim prior to alignment. For each input frame,
%               the region is provided by defining the upper left and lower
%               right corner, e.g. {[Xul,Yul; Xlr,Ylr]}. Optional, but
%               likely to improve alignment if large deviations exist in
%               frame centering (or size).
%
%   Example usage:
%   
%   vc_app = iman_append('path\filename1.mat', 'path\filename2.mat', ...
%       'xyind', {[5,6],[2,3]}, 'trim', {[100,400; 700,800],[300,200; 900,600]});
%
%   Version 2.0

%   M. Pargett, Albeck Lab, UC Davis, MCB.  2016
%   Version 2.0 update to enforce provision of vcorder to
%       iman_aligntrackdata, and removes input as cell arrays

function vc_app = iman_append(in1, in2, varargin)
%Version check provision
if strcmpi(in1,'version'); vc_app = 'v2.0'; return; end

p.xyind = cell(1,2);  %Default XY indices are empty (found as end-2:end-1)
p.trim  = {};         %Default is to not trim (empty)

%Input option pair parsing:
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end
%Splits pairs to a structure
for s = 1:2:nin;   p.(lower(varargin{s})) = varargin{s+1};   end

%% Get processed data
%Parse inputs
[vc1, xycs{1}] = sub_parse_data(in1, p.xyind{1}); 
[vc2, xycs{2}] = sub_parse_data(in2, p.xyind{2});

%Run registration and alignment
vc_app = sub_align(vc1, vc2, xycs, p.trim);

end


% --- Sub-Function: Parse Input Data --- 
function [d, xyi] = sub_parse_data(d, pxy)
if ischar(d) %IF multiple names, load all
    a = load(d,'valcube','vcorder'); %Load
    d = a.valcube;    %Pack up data arrays
    %Get X and Y channel indices
    if isfield(a, 'vcorder')
        xyi = [ find(strcmpi(a.vcorder, 'XCoord')),... %X index
                find(strcmpi(a.vcorder, 'YCoord')) ];  %Y index
        xyi = [xyi, find(strcmpi(a.vcorder, 'nArea'))];
    else xyi = []; %IF no vcorder, no xyi
    end
elseif isnumeric(d);  xyi = []; %IF no file, no xyi
else
    error('IMAN:Append:BadInput', ['Input data must be either a ',...
        'valcube array (numeric) or a string consisting of a ',...
        'processed data file name.'])
end

%Check XY indices
if isempty(xyi)
    if isempty(pxy);  error('IMAN:Append:NoXYind', ...
            ['If data are provided without xy indices (e.g. as numeric ',...
            'arrays rather than processed data file names), X and Y ',...
            'indices must be provided via the parameter xyind.']);
    else xyi = pxy; %Fill from parameters if necessary
    end
elseif ~isempty(pxy) %IF both parameters and files give indices, check
    if ~all( xyi(:) == pxy(:) )
    error('IMAN:Append:BadXYind', ['X and Y coordinate indices provided',...
        ' in the parameter xyind do not match those indicated by the ',...
        'vcorder variable found in the processed data file(s).']); 
    end
end

end


% --- Sub-Function: Perform Registration and Alignment --- 
function vc_app = sub_align(vc1, vc2, xycs, t)
%% Register coordinates from valcube arrays
%   Bring coordinates to the same orientation (that of the former set)

%Determine data sizes
nt2 = size(vc2,2);
%   Get XY coordinate indices
xyc1 = xycs{1}(1:2); xyc2 = xycs{2}(1:2);

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
nc2 = size(vc2,1);  %Number of cells in original Model (Input 2)
if isempty(t);   vc2(:, 1, xyc2) = mdl;
else vc2(:, 1, xyc2) = ctf([squeeze(vc2(:, 1, xyc2)), ones(nc2,1)], tp);
end
%Replace time dimension in scene coordinates
scn = shiftdim( shiftdim( shiftdim(scn,1), -1), 2);


%% Align coordinates across gaps
[vc_app, vi, viv] = iman_aligntrackdata(scn, {vc2(:,1,:)}, xycs{2});
%   Replace Model coordinates into output (trackdata retains scence coords)
vc_app(viv, 1, xyc2) = vc2(vi(viv), 1, xyc2);

if nt2 > 1  %IF latter data have more than one time point
    vc_app(:,2:nt2,:) = NaN; %Initialize as NaN
    %Transform all latter coordinate sets
    for s = 2:nt2
        vc2(:,s,xyc2) = ctf([squeeze(vc2(:, s, xyc2)), ones(nc2,1)], tp); 
    end 
    %Align all latter data (neglects new traces in latter data)
    vc_app(viv, 2:nt2, :) = vc2(vi(viv), 2:nt2, :);
end
end

