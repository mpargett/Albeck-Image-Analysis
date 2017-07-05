%IMAN_IMGRENAME
%   Rename image files from a defined format to the xy/z/t/c format used by
%   IMAN functions.
%
%   [nfrm, imp] = iman_imgrename(fpath, [par1, val1, ...])
%
%   Renames images found in fpath and returns the name format used (nfrm)
%   and the parsed image names (imp), which define the mapping from old
%   names to new names.
%
%   Parameters:
%   imtypes - cell array of strings defining file extensions to consider as
%       images (can be a regular expression)
%   nameformat - name format to use when parsing names (either a string
%       matching an internally defined format, or a regular expression)
%   revert - logical, if TRUE, uses the imp parameter (MUST be provided)
%       from a previous run to revert image names back to the originals
%   imp - imp output from a previous run, used for reversion

function [nfrm, imp] = iman_imgrename(fpath, varargin)
%Default parameters
p.imtypes = {'tif{1,2}', 'png', 'jpe?g'};
p.nameformat = 'collins';
p.revert = false;
p.imp = [];

%Input option pair parsing:
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end
%Splits pairs to a structure
for s = 1:2:nin;   p.(lower(varargin{s})) = varargin{s+1};   end


%% If requested, revert previously changed names
if p.revert
    if ~isempty(p.imp) && isstruct(p.imp) && isfield(p.imp, 'oldname')
        fpath = [regexprep(fpath, '[\\/]\s*$', ''),'\'];
        for s = 1:numel(p.imp)
            movefile([fpath,p.imp(s).newname], [fpath,p.imp(s).oldname]);
        end
        nfrm = []; imp = p.imp; return;
    else
        error('IMAN:BADINPUT', ['For filename reversion, "imp" ',...
            'parameter must be provided.']);
    end
end


%% Main process
%Get parse string for the current name format
pstr = get_parse_string(p.nameformat);

%Verify folder path
assert(exist(fpath, 'dir')==7, 'IMAN:BADPATH', ['Folder path provided ', ...
    'is not a valid directory (or is not accessible).']);
%   Get info on files in the folder
fpath = [regexprep(fpath, '[\\/]\s*$', ''),'\'];
dl = dir(fpath);
imn = regexpi({dl.name}, '(?<nm>[^\.]*)\.?(?<ext>[^\.]*)$', 'names');
imn = [imn{:}];
%   Determine type of image (for each file)
imi = cellfun(@(imt) ~cellfun(@isempty, regexpi({imn.ext}, imt)), ...
    p.imtypes, 'Un', 0);    imi = cat(1,imi{:})';

%Parse current image names
imind = any(imi,2);  nim = nnz(imind);
imp = regexpi({imn(imind).nm}, pstr, 'names');  imp = [imp{:}]';
[imp.oldname] = deal(dl(imind).name);

%Build new image names
[nnm{1:nim,1}] = deal('');
%   Base name
nnm = cellfun(@(x,y)[x,y], nnm, {imp.tag}', 'Un', 0);
nfrm.base = unique(nnm);

%Define mapping to new name fields
flds = {'xy',   'z',     'c',    't'    ;...
        'site', 'stack', 'chan', 'time'     };

%Map each field and new name of each image
for s = 1:size(flds,2)
    if ~isfield(imp,flds{2,s})
        fv.(flds{2,s}) = ones(nim,1); n.(flds{1,s}) = 1;  continue; 
    end
    %Get unique list of field values, and indices...
    [c, ~, ic] = unique({imp.(flds{2,s})}');
    n.(flds{1,s}) = numel(c);
    %   Retain values for each field
    fv.(flds{2,s}) = ic;
    
    if strcmpi(flds{1,s},'t')  %IF on time index
        %Time indices are mapped for each other index
        %   (Avoids complication of timestamps being different for each
        %   well/stack/channel etc.)
        %   Get all timestamps and initialize
        tstamps = {imp.(flds{2,s})}';  ic = nan(nim,1);
        for sx = 1:n.xy      %FOR each xy site
            for sz = 1:n.z       %FOR each z-stack
                for sc = 1:n.c       %FOR each channel
                    %Get current image indices
                    cii = fv.site == sx & fv.chan == sc & fv.stack == sz;
                    %Get sorted (unique) time stamps
                    [c, ~, ict] = unique(tstamps(cii));
                    %Arrange time indices in full set
                    ic(cii) = ict;
                    %Store map of timestamps over other indices
                    nfrm.(flds{1,s})(sx,sz,sc,:) = c;
                end
            end
        end
                
    else    %IF not operating on time
        %   Store ordered list of old field values
        nfrm.(flds{1,s}) = c;
    end
    %Append new field to new names
    nnm = cellfun(@(x,ci)[x,flds{1,s},num2str(ci)], ...
        nnm, num2cell(ic), 'Un', 0);
end

%Append file extensions
nnm = cellfun(@(x,y)[x,'.',y], nnm, {imn(imind).ext}', 'Un', 0);

%Store new names
[imp.newname] = deal(nnm{:});

%Rename files
for s = 1:nim
    movefile([fpath,imp(s).oldname], [fpath,imp(s).newname]);
end

%Save renaming info, for reference, or reversion
save([fpath,'imgrename_info.mat'], 'nfrm', 'imp');

end


%SUBFUNCTION: Parse string definitions
function pstring = get_parse_string(in)

switch lower(in)
    case 'collins'
        %Image name parsing strategy for Collins lab data:
        %   tag_AQID_SITE_CHANNEL_yyymmddThhmmss
        pstring = ['(?<tag>\w*)_(?<aqid>\d*)_(?<site>[a-zA-Z]*\d*)_',...
            '(?<chan>[\w-]*)_(?<time>\d*T\d*)'];
    otherwise
        %Pass input through if not recognized.  May be a parse string.
        pstring = in;
end


end