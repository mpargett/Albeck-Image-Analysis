%IMAN_IMAGEACCESS
%   Image data access procedure for celltrace processing.
%
%   [dao] = iman_imageaccess(in);
%       Returns a Data Access Object, dao, built for the string input, in.
%   in must specify a valid image file (e.g. ND2 format), or folder in
%   which images are stored (e.f. TIF format).  dao can be passed to
%   iman_getframe, with a vector of index values to retrieve image data
%   from desired frames.
%
%   When a single standalone image (TIF, etc.) is specified, the containing
%   folder is searched for  further images of the same type and naming
%   pattern.
%
%   For standalone images, image naming patterns must conform to [BaseName,
%   IDsequence, Appendix], where IDsequence is constructed of the
%   dimensions used (abbreviated as t,c,xy,z for Time,Channel,XY Position,
%   Z Position) in any order, each followed by the index value for that
%   image.  For example, TestImaget01c2xy12.tif is the 1st Time point, 2nd
%   Channel, and 12th XY position, with no Z stack.  Spacers between
%   indices may be -,_,or whitespace.  Ensure that exported images are
%   named such that all values for a given index have the same number of
%   characters, e.g. t001 for the 1st point if Time goes up to t100. 
%   


%FIXME Also need a metadata handling routine

function [dao] = iman_imageaccess(in, indsz, varargin)
%Version check provision
if strcmpi(in,'version'); dao = 'v1.0'; return; end

%% Set defaults and initialize
imtyp = {'tif','jpg','png'};  %Image types to consider
dao = struct('type',[], 'read',[], 'rinfo', [], 'info', [], 'abort', false);
if ~exist('indsz','var'); indsz = []; end

%% Determine data type
if ischar(in)       %Input is a filename
    %Check existence and type of name
    exflag = exist(in); %#ok<EXIST>
    
    %Prep info based on name type
    switch exflag
        case 2  %Input is a filename
            fdef = regexp(in, ['^(?<loc>.*(\\|/))(?<nam>[^\\/]*)',...
                '\.(?<ext>\S*)$'], 'names');
            switch lower(fdef.ext)
                case 'nd2'; dao.type = 'nd2';
                case imtyp; dao.type = 'image'; dao.info.ftype = fdef.ext;
            end
            
        case 7  %Input is a folder name
            fdef.loc = in;  dao.type = 'image';  %Assume folder = images
            %Ensure loc ends with \
            if ~any(strcmp(fdef.loc(end), {'\','/'}))
                fdef.loc = [fdef.loc,'\'];
            end
            %Type of file not yet determined
            dao.info.ftype = [];
            
        otherwise  %Error on bad input
            error('IMAN:BADFILENAME',['The input name is not a valid',...
                ' file or folder.  File names must include the extension.']);
    end
end

%% Fill DAO based on type of input
switch dao.type
    case 'nd2'  %IF ND2 file, run bfopen (modified)
        dao = access_nd2(dao, in, indsz);        
    case 'image'
        dao = access_imfile(dao, fdef, indsz, exflag, imtyp);
end

end


% ----------------------- FILE ACCESS SUBFUNCTIONS -----------------------
%% ND2 files
function [dao] = access_nd2(dao, in, indsz)
%Get Reader and associated info
[dao.read, dao.rinfo, rp] = bfopen_custom(in);
%Store extra Reader parameters (includes memopath info)
dao.info.rp = rp;

%Determine XY storage order (sometimes as Series, sometimes as z)
nd2sz = [numel(dao.rinfo), dao.rinfo(1).ZTCsize];
iford = {'xy','z','t','c'};
if isstruct(indsz)
    %Ensure incoming fields are properly ordered, robustly
    itmp = indsz; indsz = cell2struct(cell(1,4),iford,2);
    for s = intersect(iford,fieldnames(itmp))'
        indsz.(s{1}) = itmp.(s{1}); end
    %Ensure all fields are filled (assume empty = 1)
    indsz = structfun(@(x)max(x, isempty(x)), indsz, ...
        'UniformOutput', false);
    %Compare data sizes to expected sizes
    xys = nd2sz(1) == indsz.xy;
    isvalid = nd2sz([2-xys, 1+xys,3:end]) == ...
        cell2mat(struct2cell(indsz))';
    %IF any are not valid, report. Return for Abort if not fixable.
    if ~all(isvalid);
        %Build report of which fields do not match
        itmp = iford(~isvalid);  istr = itmp{1}; itmp(1) = [];
        itmp = cellfun(@(x)[', ',x], itmp,'UniformOutput',false);
        istr = [istr, itmp{:}];
        %Double-check if XY stacked in Time
        if all(nd2sz(1:2) == 1) && nd2sz(3) == indsz.t*indsz.xy;
            warning('IMAN:XYinTime',['Time dimension size is ',...
                'product of expected XY and Time sizes.  ',...
                'Treating as XY layering into Time.']);
            nt_keep = floor(dao.rinfo.ZTCsize(2)./tnxy);
            if nt_keep < indsz.t;
                warning('IMAN:shortTime', ['Maximum time shorter', ...
                    ' than expected.  Truncating run to ',...
                    num2str(nt_keep), ' time points.']);
            else nt_keep = indsz.t; %Do not exceed indsz.t
            end
            %Adjust indexing for XY laying in Time
            dao.rinfo.lblr = reshape(dao.rinfo.lblr(:,1:nt_keep*indsz.xy,:), ...
                indsz.xy, nt_keep, indsz.c);
            dao.rinfo.ZTCsize = size(dao.rinfo.lblr);
        else
            %Issue warning and Return for Abort
            warning('IMAN:ReaderDimensionMismatch', ['The dimension',...
                ' of the ', istr,' indices do not match the ',...
                'provided values.  Aborting and returning ',...
                'Reader.']);  dao.info.ThreadSafe = false;
                dao.abort = true; return;
        end
    end
else %IF no expected index sizes available
    %Assume multiple Series are XYs and if singular, z is XY
    if nd2sz(1) > 1;  xys = true;  else xys = false;  end
end
%Store index ordering and max index values
% dao.info.iord = iford;  dao.info.imax = nd2sz([2-xys, 1+xys,3:end]);
dao.info.imax = cell2struct(num2cell(nd2sz([2-xys, 1+xys,3:end])), iford, 2);
dao.info.xySeries = xys;  dao.info.ThreadSafe = false;
end


%% Image files
function dao = access_imfile(dao, fdef, indsz, exflag, imtyp)
%Get file names/types in the folder
%   Determine which files are to be processed
dirtst = dir(fdef.loc);
dirinf = regexp({dirtst.name}, ...
    ['(?<fname>[^\.]*)','\.?(?<ext>.*)'], 'names'); clear dirtst;
dirfn = cellfun(@(x)x.fname, dirinf, 'UniformOutput', false);
dirext = cellfun(@(x)x.ext, dirinf, 'UniformOutput', false);
clear dirinf;

targetim = [];
switch exflag
    case 2  %Input is a file
        imt = dao.info.ftype;  targetim = fdef.nam;
        %FIXME check on if multipage image
        
    case 7  %Input is a folder
        %Get only valid image types
        imt = intersect(unique(dirext), imtyp);
        %IF multiple image types found, select most common
        if numel(imt) > 1
            nmt = numel(imt); nper = zeros(1,nmt);
            for s = 1:nmt; nper(s) = nnz(strcmp(imt{s}, dirext)); end
            [~, mm] = max(nper);  imt = imt{mm};
        else imt = imt{1};
        end;    dao.info.ftype = imt;
end

%Get filenames of proper filetype
fns = dirfn(strcmp(dirext, imt));  clear dirfn dirext;
%IF no target file specified, use first in listing
if isempty(targetim); targetim = fns{1}; end

%Get image file info for target image
imfi = imfinfo([fdef.loc,targetim,'.',imt]);

%Store basic image properties
dao.info.color = imfi(1).ColorType;
dao.info.xydim = [imfi(1).Width, imfi(1).Height];
dao.info.nchan = imfi(1).SamplesPerPixel;
dao.info.bits = imfi(1).BitsPerSample;

%Record image storage type
dao.info.multipage = numel(imfi) > 1;
dao.info.multichan = imfi(1).SamplesPerPixel > 1;

%Establish file name structure for frame IDs
[fid,fbeg] = regexp(targetim, '([-_\s]?(t|xy|z|c)\d+)*', 'match', 'split');
%Parse non-id segments into beginning and end (if any)
if ~isempty(fid);   fend = fbeg{2};
    %Parse id segment into spacers, ids, and values
    fid = regexp(fid{1}, ...
        '(?<spc>[-_\s]?)(?<id>(t|xy|z|c))(?<val>\d+)', 'names');
else   fend = [];
end;   fbeg = fbeg{1};

 

%Construct ID matching specification, based on naming pattern
idmatch = [];   fnum = struct();
for s = 1:numel(fid)
    idmatch = [idmatch, fid(s).spc, fid(s).id,'(?<',fid(s).id,'>\d+)']; %#ok<AGROW>
end

%Parse valid indices for valid files
fstr = cell2mat(regexp(fns, ['^',fbeg,idmatch,fend], 'names')); clear fns;
for s = 1:numel(fid)
    fnum.(fid(s).id) = str2double({fstr.(fid(s).id)})';
end
%Get field widths for each index              %(below)Fixes size of empties
fw = cellfun(@numel,struct2cell(fstr(:)));  %Widths for each
%   Compare to see if pre-padded with zeros, make sprintf format value
fw = 1 + all(all(bsxfun(@eq, fw(:,1),fw),2)).*(fw(:,1)'-1);  
clear fstr; if isempty(fw); fw = []; end

%Store index order and maximum index values
dao.info.imax = structfun(@max, fnum, 'UniformOutput',false);  clear fnum;
dao.info.ThreadSafe = true;
%Assemble file name creation function
dao.read = @(x)cat( 2, fdef.loc, fbeg, cell2mat(arrayfun(...
    @(z,xx,w)[z.spc,z.id,sprintf(['%0',num2str(w),'d'],xx)], ...
    fid, x, fw, 'UniformOutput', false)), fend, '.', imt   );

%Collect image information for each image (imfi)
fdims = cell2mat(struct2cell(dao.info.imax));   %Get number of files
nfd = numel(fdims);  imfi = cell([fdims(:);1]');  %Allocate imfi
for s = 1:prod(fdims);
    [oi{1:nfd}] = ind2sub(fdims,s);  oi = [oi{:}]; %Get indices
    if isempty(oi); oi = []; end  %Fixes size of empties
    imfi{s} = imfinfo(dao.read(oi));  clear oi;
end

%IF a multi-page image, provide necessary indexing information
if dao.info.multipage 
    dao.info.imfi = imfi;  %Store full image information
    %Identify and store dimension(s) rolled in the pages
    if isfield(indsz,'mpdim');  dao.info.mpdim = mpdim;
    else 
        fdn = fieldnames(dao.info.imax); 	%File Dimension Names
        mpsz = rmfield(indsz, fdn);         %expected MultiPage Size
        mpdn = fieldnames(mpsz);            %MultiPage Dimension Names
        mpsv = cell2mat(struct2cell(mpsz)); %MultiPage Size Vector
        
        %Validate MultiPage Size
        assert(prod(mpsv) == numel(imfi{1}),'IMAN:ImagePageCheck', ...
            ['Identified dimensions (t,c,xy,z) of the apparently ',...
            'multi-page image do not account for the number of image ',...
            'pages. Recheck image structuring and provide dimension ',...
            ' size(s) as indsz.mpdim.']);
        
        %Assert expected page dimension order (t > xy > c) z??
        expord = {'t', 'xy', 'c'};
        mpord = cellfun(@(x)find(strcmpi(mpdn, x)), expord, ...
            'UniformOutput', false);  mpord = [mpord{:}];
        dao.info.mpdim = orderfields(mpsz, ...
            [mpord, setdiff(1:numel(mpdn),mpord)]);
    end
end

end


