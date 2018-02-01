%IMAN_FRAMESTATS
%   Extracts full frame statistics from and image Data Access Object.
%
%   y = iman_framestats(dao, ...)
%       returns statistics (sum, mean, var) for the imagery in the Data
%       Access Object dao (obtained via iman_imageaccess).  y is a
%       structure with fields for requested statistics, each field contains
%       a XY x Time x Channel x Z array of values, corresponding with the
%       requested frames (i.e. in the order provided).  Frames are
%       selected via additional inputs.  
%
%   Additional inputs may be provided as Name/Value pairs:
%
%   t - Time point(s) to include (vector)
%   c - Channel(s) to include (vector)
%   xy - XY position(s) to include (vector)
%   z - Z position(s) to include (vector)
%   trim - Percentiles at which to trim tails of intensity distribution 
%           (1x2 vector, e.g. [5,95])
%   sum - Locigal, true to return sum (default = TRUE)
%   mean - Locigal, true to return mean (default = TRUE)
%   var - Locigal, true to return variance (default = TRUE)
%   reg - Region of image to use (optional), 2x2 array defining a box
%       (UpperLeftX, UpperLeftY; LowerRightX, LowerRightY)

%FIXME - Non-Uniform background subtraction, when ready

function y = iman_framestats(dao, varargin)
%Default parameters
p = struct('t',1,'c',1,'xy',1,'z',1,'trim', [], ...
    'sum',true,'mean',true,'var',true,'reg',[]);

%Input option pair parsing:
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end
%Splits pairs to a structure
for s = 1:2:nin;   p.(lower(varargin{s})) = varargin{s+1};   end

%Set flag for region restriction
if isempty(p.reg); rg = false; else rg = true; end

%% Process frames
%Initialize
n = [numel(p.xy), numel(p.t), numel(p.c), numel(p.z)];  
flds = {'sum';'mean';'var'};  uf = [p.sum;p.mean;p.var];
y = cell2struct(repmat({nan(n)},nnz(uf),1), flds(uf), 1);
%Iterate through each frame
for sx = 1:n(1)
    for st = 1:n(2)
        for sc = 1:n(3) 
            for sz = 1:n(4)
            %Pull frame
            im = iman_getframe(dao,[p.t(st),p.c(sc),p.xy(sx),p.z(sz)]);
            %Restrict region, if requested
            if rg; im = im(p.reg(3):p.reg(4), p.reg(1):p.reg(2)); end
            %Trim distribution of values, if requested
            if ~isempty(p.trim)
                pt = prctile(im(:), p.trim);
                im = im > pt(1) & im < pt(2);
            end
            %Take stats
            if p.sum;   y.sum(sx,st,sc,sz) = sum(im(:));    end
            if p.mean;  y.mean(sx,st,sc,sz) = mean(im(:));  end
            if p.var;   y.var(sx,st,sc,sz) = var(im(:));    end
            end
        end
    end
end

end

