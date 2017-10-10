%CELLTRACER_V2_PARTNER
%   Accompaniment functions for Celltracer run scripts.

function [varargout] = celltracer_v2_partner(flag, varargin)
switch lower(flag)
    case 'validate'
        [varargout{1:2}] = celltracer_validation(varargin{1:2});
    case 'initialize'        
        [varargout{1:2}] = celltracer_initialize_parameters();
    case 'match_names'
        [varargout{1}] = match_names(varargin{1:2});
end
end

%% Pre-Run validation
function [ip,op] = celltracer_validation(ip,op)
% --- Version verification --- 
%Call for version of main celltracer function
[ctv, rqv] = iman_celltracer('version');
fprintf(['Checking code versions for iman_celltracer ', ctv,' ...']);
fnn = fieldnames(rqv);  %Get list of functions with version requirements
%Call for version of each image analysis function
for s = 1:numel(fnn)  %FOR each function to be checked
    vs = regexpi(eval(['iman_',fnn{s},'(''version'')']), ...
        'v(?<main>\d+)\.(?<sub>\d+)', 'names');  %Get and format versions
    vmain.(fnn{s}) = str2double(vs.main);   %Main version numbers
    vsub.(fnn{s}) = str2double(vs.sub); %#ok<STRNU> %Sub version numbers
end;  vcheck = struct2cell(vmain);
%Compare required versions with those available and error if incompatible
assert(all( cellfun(@(x,y)any(x==y), vcheck, struct2cell(rqv)) ), ...
    'IMAN:VersionCheck', ['Code versions do not match.  ',...
    'The main versions self-reported by each function must match ',...
    'those required by iman_celltracer.  See the rqv and vmain ',...
    'structures for requirements and currently available versions.']);
fprintf(' OK\n');
% --- ---------------------------------------------------------------- ---

fprintf('Validating image and run parameter definitions ...');
% --- Parameter validation ---
%Match/check channel names/indices to use (op.cind, allow string input)
if isempty(op.cind);    op.cind = 1:numel(op.cname);	end
[op.cind] = name_map(op.cind, op.cname, 'Index');
assert(isnumeric(op.cind) && max(op.cind) <= numel(op.cname), ...
    'IMAN:cindCheck', ['Validation failed. Channel indices to use ',...
    'exceed number of declared channels. Verify op.cind against ',...
    'op.cname/ip.bkmd.exp.Channel']);

%Match/check segmentation channel (op.seg.chan, allow string input)
[op.seg.chan] = name_map(op.seg.chan, op.cname, 'Segmentation', op.cind);

%Match/check aggregation channels (op.seg.chan, allow string input)
if ~isempty(op.msk.aggfun)
[op.msk.aggfun.chan] = name_map(op.msk.aggfun.chan, op.cname, ...
    'AggFun', op.cind);
[op.msk.aggfun.loc] = name_map(op.msk.aggfun.loc, {'Nuc','Cyt'}, ...
    'AggFunLoc');
end

%Match/check requested pre-averaging ratios
for s = 1:numel(op.msk.rt)  %Check each element
    if iscell(op.msk.rt{s})     %Check if is cell (typical)
        if ischar(op.msk.rt{s}{1})  %Check if char names
            [op.msk.rt{s}] = name_map(op.msk.rt{s}, op.cname,...
                'Ratio', op.cind);  %Map to channel indices
        elseif isnumeric(op.msk.rt{s}{1}); op.msk.rt{s} = [op.msk.rt{s}{:}];
        end
    elseif isempty(op.msk.rt{s}); op.msk.rt(s) = [];  %Remove empty entries
    end
end
%Check index ranges requested
if isempty(op.trng); op.trng  = [1, ip.indsz.t]; end
assert( max(op.cind) <= ip.indsz.c && max(op.xypos) <= ip.indsz.xy && ...
    max(op.trng) <= ip.indsz.t, 'IMAN:MaxIndexCheck', ['Validation failed.  ',...
    'An index requested for channel, xy or time (in op) exceeds the ',...
    'maximum index defined for the file (in ip).']);
% --- ---------------------------------------------------------------- ---

% --- Background review ---
ip.bkg = ip.bkg(op.xypos);   	%Restrict background to XYs used
%Assert proper sizes of background region definitions
bkc = arrayfun(@(x)~isempty(x.altxy) || (x.fix && numel(x.reg) == ip.indsz.c) || ...
    (numel(x.reg) == 4 && x.reg(1)<x.reg(2) && x.reg(3)<x.reg(4)), ip.bkg);
assert(all(bkc), 'IMAN:bkgCheck', sprintf(['Validation failed. ',...
    'Background regions for XY position(s) ', repmat('%d, ', 1, ...
    nnz(~bkc)-1), repmat('and ', 1, nnz(~bkc)>1), '%d',' are improperly ',...
    'sized or defined. Review ip.bkg for compliance.'], op.xypos(~bkc)));
%Rescrict any fixed values to current channels
for sb = find([ip.bkg.fix]); ip.bkg(sb).reg = ip.bkg(sb).reg(op.cind); end
% --- ---------------------------------------------------------------- ---

% --- Backup MetaData review ---
%   Ensure filled sample time in MetaData
if isfield(ip,'tsamp'); 
    ip.bkmd.cam.tsamp = ip.tsamp; % keep for historical reasons
       op.trk.linkwin = ceil(op.trk.linkwin/ip.tsamp); % change from minutes to frames
elseif ~isfield(ip.bkmd.cam,'tsamp') || isempty(ip.bkmd.cam.tsamp)
    error('IMAN:tsamp', ['Either ip.tsamp or ip.bkmd.cam.tsamp must ',...
        'be defined as a scalar value.']);
end

%Review for Spectral Unmixing procedures
if op.unmix
    %   Transfer Channel IDs to MetaData for passage to unmixing
    if numel(op.cind)<numel(op.cname); ip.bkmd.exp.cind=op.cind; end
    %Validate naming for Fluorophores and Filters
    nmsc = [ip.bkmd.exp.Filter, ip.bkmd.exp.FPhore];
    nv = iman_naming('validate', nmsc{:});
    assert(all(nv), 'IMAN:bkmdCheck', ['Validation failed. The following ',...
        'names provided in MetaData do not match allowed names: ', ...
        repmat('%s, ', nnz(~nv)-1), '%s.  Call [fpn, ftn] = iman_naming ',...
        'to see allowed names.'], nmsc{~nv});
    %Validate consistent number of Channels
    assert(isequal(numel(ip.bkmd.exp.Channel), numel(ip.bkmd.exp.Filter), ...
        numel(ip.bkmd.exp.FPhore)), 'IMAN:bkmdCheck', ['Validation failed. ',...
        'The number of names for Channels, Filters, and Fluorophores in ',...
        'ip.bkmd.exp do not match.']);
    %Validate FRET names match Channel names
    if ~isempty(ip.bkmd.exp.FRET)
        frtn = struct2cell(ip.bkmd.exp.FRET); frtn = [frtn{:}];
        frtv = cellfun(@(x)any(strcmp(x, ip.bkmd.exp.Channel)), frtn);
        assert(all(frtv), 'IMAN:bkmdCheck', ['Validation failed. The following ',...
            'FRET component names to not match Channel names in ip.bkmd.exp: ', ...
            repmat('%s, ', nnz(~frtv)-1), '%s.'], frtn{~frtv});
    end
end
%   Validate if specifying SPECTRAX (multi-line) light source
if any(strcmpi(ip.bkmd.exp.Light, {'SPECTRAX'}))  %IF a multi-line source
    mln = {'ExVolt', 'ExLine', 'ExWL'};
    exc = cellfun(@(x)iscell(ip.bkmd.exp.(x)), mln);
    if all(exc); exc = cellfun(@(x)cellfun(@(y)numel(y), ip.bkmd.exp.(x)),...
            mln, 'UniformOutput', false); else exc = {1,2,3}; end
    assert(isequal(exc{:}), 'IMAN:bkmdCheck', ['Validation failed. ',...
        'When specifying a multi-line light source (e.g. SPECTRAX), ',...
        'ExVolt, ExLine and ExWL must all be cell array with the same ',...
        'number of elements and the size of each entry must match among ',...
        'the three fields.']);
end
% --- ---------------------------------------------------------------- ---

% --- XY Shift definition conversion ---
%   Automated conversion of xyshift time values for indicated time range
if ~isempty(ip.xyshift)     %IF a shift is defined
    if ~isempty(op.trng)        %IF the time range is restricted
        badi = ip.xyshift.frame < op.trng(1)+1 | ...
            ip.xyshift.frame > op.trng(end);         ni = numel(badi);
        for s = {'frame','dx','dy'};    %Remove out of range shifts
            if numel(ip.xyshift.(s{1})) == ni
                ip.xyshift.(s{1})(badi) = []; end
        end
        ip.xyshift.stime = op.trng(1); %Indicate starting time point
    else    ip.xyshift.stime = 1;   %Indicate if start time is beginning
    end
    %Match a tracking channel name to proper index
    if ~isempty(ip.xyshift.trackchan) && ischar(ip.xyshift.trackchan)
        [ip.xyshift.trackchan, erm] = match_names(ip.xyshift.trackchan, ...
            op.cname);
        if ~isempty(erm); fprintf('\n');  warning('IMAN:ShiftCheck', ...
                ['Tracking channel name not found. ', erm]); end
    end
end
% --- ---------------------------------------------------------------- ---

fprintf(' OK\n');

%Path management
fprintf('Applying path definitions for Java ...');
op.pth = iman_assertjavapaths;
fprintf(' OK\n');
end


%% Parameter structure initialization
function [ip, op] = celltracer_initialize_parameters()
%Image parameters
ipfields = {'fname', 'sname', 'indsz', 'bkg', 'xyshift', 'bkmd', 'bval'};
ip = cell2struct( cell(size(ipfields)), ipfields, 2 );
    ip.indsz = struct('xy', [], 'z', [], 't', [], 'c', []);
    ip.bkg   = struct('reg',[], 'fix',[], 'dyn',[], 'altxy',[], 'bkonly',[]);
    ip.xyshift = struct('frame',[], 'dx',[], 'dy',[], ...
        'trackwell',[], 'trackchan',[], 'badtime',[], ...
        'shifttype', 'translation');
    ip.bkmd = [];

%Operation parameters
opfields = {'cname', 'cind', 'xypos', 'trng', 'nW', 'objbias', 'unmix', ...
            'fixshift', 'mdover', 'seg', 'msk', 'trk', 'disp', 'pth'};
op = cell2struct( cell(size(opfields)), opfields, 2 );
    op.seg = struct('chan',[], 'cyt',false, 'maxD',[], 'minD', [], ...
                    'maxEcc',[], 'Extent',[]);
    op.msk = struct('fret',[], 'freti',[], 'rt',{{}}, 'storemasks',false, ...
                    'saverawvals',false, 'aggfun', []);
    op.trk = struct('movrad', [], 'linkwin', []);
    op.disp = struct('meta',false, 'samples',false, 'shifts', false, ...
                     'warnings',false);
end


%% Name matching
function [ind, erm] = match_names(in, ref)
if ischar(in); in = {in}; end   %Ensure input is in a cell array

%Match names by comparing each with the reference set
nin = numel(in);  ind = nan(1,nin);  %Get size and initialize
for s = 1:nin
    tmp = find( strcmpi(in{s}, ref) );   
    if isempty(tmp); continue; end;    ind(s) = tmp;
end

%Check for validity and return error message if check failed for any
ni = isnan(ind);
if any(ni); erm = sprintf(['The following name(s) were not matched: ', ...
        repmat('%s, ',1, numel(ni)-1), '%s'], in{ni});
else erm = '';
end

end

%% Channel name checking/remapping procedure
function [cin] = name_map(cin, cref, ctype, cmap)
%Match/check segmentation channel (op.seg.chan, allow string input)
if ischar(cin) || iscell(cin)
    [cin, erm] = match_names(cin, cref); 
    assert(isnumeric(cin) && isempty(erm), ['IMAN:',ctype,'ChanCheck'], ...
        ['Validation failed. ',ctype,' channel name does not ',...
        'match declared names. ', erm]);
elseif isempty(cin) || ~isnumeric(cin)
    error(['IMAN:',ctype,'ChanCheck'], ['Validation failed. ',ctype,...
        ' channel must be numeric or a string.']);
end
%   Remap channel to reduced indices, as needed
if exist('cmap', 'var')
    [chk,cin] = ismember(cin, cmap);
    assert(all(chk), ['IMAN:',ctype,'ChanCheck'], ['Validation ',...
        'failed. ',ctype,' channel index not found within selected ',...
        'channel indices. Verify name or index provided against ',...
        'those in op.cind']);
end

end

