%IMAN_READDATASHEET
%   Reads properly formatted imaging datasheet and extracts metadata to 
%   a structure contaning plate-sized cell-arrays and an indexing structure
%   with matched lists of XY-position and values for Cells and Treatment
%
%   [pmd, idx] = iman_readdatasheet(f);
%       returns the plate metadata (pmd) and indexing structure (idx) for
%   datasheet referred to by filename f.  f may be a full path if the file
%   is not on the MATLAB path.
%
%   pmd contains fields available on the datasheet template, each
%   containing a cell array wherein the 3rd dimension stacks parts of
%   metadata depending on the field (as follows):
%       XY - reduced to vectors in cell array over plate
%   	Cell - {Name, Density, Units, Repeat...} per well
%       Genotype - {Name} per well
%       Pre-Tx - {Reagent, Conc., Unit, Time, Unit} per well
%       Tx - {Reagent, Conc., Unit, Time, Unit, Repeat...} per well
%
%   idx is a structure containg XY indices for different cell types and
%   treatments, with matched data for treatments (dose, timing).  Each cell
%   type or treatment component is given a new field in the structure.
%   Cell types are the concatenated 'Cell' and 'Gene' from pmd.
%
%   Example idx fields:
%   idx.Cell.mcf10Aekar3: [1x10 double]  <-List of XYs with these cells
%   idx.pTx.Gef.xy: [1x20 double]   <- List of XYs pre-treated with Gef
%              .dose: [1x20 double] <- Gef dose (value only)
%              .dunit: 'nM'         <- Gef does unit (for all XYs)
%              .time: [1x20 double] <- Time relative to imaging (value)
%              .tunit: 'tp'         <- Time units

function [pmd, idx] = iman_readdatasheet(f, varargin)
%% Input parsing
p.block = false;
%Input option pair parsing:
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end
%Splits pairs to a structure
for s = 1:2:nin;   p.(lower(varargin{s})) = varargin{s+1};   end



%% Load datasheet and Parse headers (accommodate added lines etc.)
%Read excel sheet (specified by f)
[n, t, r] = xlsread(f); %#ok<ASGLU>

%   Discard leading empty columns or rows
cfill = false;  rfill = false;
while ~cfill || ~rfill
    cfill = cfill || ~all(isbad(r(:,1)));    
    rfill = rfill || ~all(isbad(r(1,:)));
    r = r( (1 + ~rfill) : end, (1 + ~cfill) : end );
end;

%   Replace NaN cells with empty strings
br = isbad(r); [r{br}] = deal('');

%   Find first main header (Movie Name)
mnhi = find(~cellfun( @isempty, regexpi(r(1:10,1), 'movie name') ));
hdi = mnhi:(mnhi+10);
%   Assert remaining main headers
hds = {'movie name', 'fName'; 'directory', 'fDir'; ...
       'processed data [dn]ame', 'pName'; 'processed data dir', 'pDir'; ...
       'plate type', 'Plate'; 'plating cond', 'Cond'; ...
       'cell density', 'Density'; 'sampling time', 'Timing'; ...
       'movie length', 'Length'; 'experimentalist', 'Person'; ...
       'notes', 'Notes'};
   
assert( all(~cellfun( @isempty, regexpi(r(hdi,1), hds(:,1)) )), ...
    'DataSheetReader:BadHeaders', 'Main information headers do not ', ...
    'match expected labels. Check sheet against template.');
   
%   Find plate headers (XY, Cell, Geno, Pre-Tx)
    phds = {'xy', 'cell type', 'genotype', 'pre-treatment$'; ...
            'xy', 'Cell',      'Gene',     'pTx'};
    ph = cell2struct(cell(size(phds(2,:))), phds(2,:), 2);
    coltext = cellfun(@num2str, (r(hdi(end)+1:end,2)), 'Un', 0);
    for s = 1:size(phds,2)
        ph.(phds{2,s}) = hdi(end) + find(~cellfun( @isempty, ...
            regexpi( coltext, phds{1,s}) ));
    end
    
%   Find plate headers for Treatments
    tmp = regexpi( coltext, '^treatment (\d*)', 'tokens');
    tmpi = hdi(end) + find(~cellfun( @isempty, tmp ));
    for s = 1:numel(tmpi)
        ph.(['Tx',tmp{tmpi(s)-hdi(end)}{1}{1}]) = tmpi(s);
    end
    pfn = fieldnames(ph);   %Keep fieldnames handy
    

%% Allocate Metadata
%Get unit definitions
[ut, tu, cu, du] = defineunits(pfn); %#ok<ASGLU>

%Allocate main information
for s = 1:numel(hdi)
    switch lower(hds{s,2})
        case {'timing', 'length'}
            d.(hds{s,2}) = eparse(r{hdi(s),2}, 'basic', tu);
        case 'density'
            d.(hds{s,2}) = eparse(r{hdi(s),2}, 'basic', du);
        otherwise;  d.(hds{s,2}) = r{hdi(s),2};
    end
end

%Get valid plate range from XY
pr = {1:8, 2:13};  rtmp = r(ph.xy + pr{1}, pr{2});
vpr = ~cellfun( @isempty, r(ph.xy + pr{1}, pr{2}) );
%   Remove any entries with only space characters
    strchk = vpr & cellfun(@ischar,rtmp);
    vpr(strchk) = ~cellfun(@isempty, regexpi(rtmp(strchk), '\S'));
%Identify size of plate
p24 = false(size(vpr));     p24(1:4,1:6) = true;
if ~any(vpr(~p24));  vpr = reshape(vpr(p24),4,6);  pr = {1:4, 2:7};   end
p0 = cell(size(vpr));

%Pre-allocate plate metadata structure
pmd = cell2struct([{[]};repmat({p0},numel(pfn),1)], [{'src'};pfn], 1);
%   Parse filename for reporting
pmd.src = regexpi(f, '(?<path>^.*\\)?(?<fname>[^\\]*$)', 'names');

%Assign metadata for each field
%   Assign unit type for each field
for s = 1:numel(pfn)
    tmp = cellfun(@num2str, r(ph.(pfn{s}) + pr{1}, pr{2}), ...
        'UniformOutput', false);
    pmd.(pfn{s})(vpr) = eparse(tmp(vpr), pfn{s}, ut{s});
    
    %   Orient per well cells to for sliceable matrix
    %   Get max number of entries (in cells only)
    nc = max([1;cellfun(@(x)iscell(x).*numel(x), pmd.(pfn{s})(vpr))] );  
    pmd.(pfn{s})(~vpr) = cellfun(@(x)num2cell(nan(1,nc)), ...
        pmd.(pfn{s})(~vpr), 'UniformOutput', false);
                
    %Map empty spaces by upper-left corners
    %   Identify remaining empties (including empty strings in cells)
    ec = cellfun( @(x)isempty(x) || (iscell(x) && ...
            all(cellfun(@isempty,x))), pmd.(pfn{s}) );
    if all(ec(:)|~vpr(:))
        pmd = rmfield(pmd, pfn{s});  continue;  %Skip and remove if empty
    elseif any(ec(:))
        if p.block
            %Identify block structures
            cdt = vpr & ~ec;    %Candidate seeds
            %   Remove buried entried (not empty below or to right)
            cdt( [cdt(:,2:end),false(size(cdt,1),1)] & ...
                [cdt(2:end,:);false(1,size(cdt,2))] ) = false;
            %   Candidates with falses to right or down
            [cdti{1},cdti{2}] = find(cdt); cdti = cell2mat(cdti);
            ncdt = size(cdti,1);  blk = repmat({false(size(vpr))},ncdt,1);
            for sc = 1:ncdt
                %Build block associated with each seed
                blk{sc}( cdti(sc,1) + (0:find(...   %Look down
                    [cdt(cdti(sc,1)+1:end, cdti(sc,2));true],1)-1), ...
                    cdti(sc,2) + (0:find(...   %Look across
                    [cdt(cdti(sc,1), cdti(sc,2)+1:end),true],1)-1) ) = true;
                %Include only valid plate range
                blk{sc} = blk{sc} & vpr;
                %Propagate seeds
                [pmd.(pfn{s}){blk{sc}}] = ...
                    deal( pmd.(pfn{s}){cdti(sc,1), cdti(sc,2)} );
            end
            %Check for un-allocated empties
            if any(any( vpr & ~any(cat(3,blk{:}), 3) ))
                %Return error regarding plate map completeness
                error('IMAN:DSHEET', ['Unassigned datasheet entries for ',...
                    pfn{s}, ', even considering block structures. Check ',...
                    'sheet and assign a value to all wells imaged.']);
            end;  clear cdti;
        else
            %Allocate empties
            
        end
    end
    %Cat third dimension to a single cell layer
    pmd.(pfn{s}) = cellfun(@(x)reshape(x,1,1,numel(x)), ...
        pmd.(pfn{s}), 'UniformOutput', false);
    %	Expand short cell arrays to match longest (3rd dim)
    mxl = max(max(cellfun(@(x)size(x,3), pmd.(pfn{s}))));
    pmd.(pfn{s}) = cellfun(@(x)cat(3,x,cell(1,1,mxl-size(x,3))), ...
        pmd.(pfn{s}), 'UniformOutput', false);
    %   Cat nested cells to 3D array
    pmd.(pfn{s}) = reshape([pmd.(pfn{s}){:}], [size(pmd.(pfn{s})), nc]);
end
pfn = fieldnames(pmd);  %Reassign field names
%   Remove any newfound empty xy wells (e.g. if bg only well)
vpr( cellfun(@isempty, pmd.xy) ) = false;

%% Build indexing structure
idx.src = pmd.src;  %Copy source information to index structure

%FOR Cell
%   Cat celltype and genotype with "_", removing spaces/bad chars in each
%   Get number of celltypes
nct = ceil(size(pmd.Cell,3)/3);  ci = 3*((1:nct)-1) + 1;

%   Cat genes together (DOUBLE CHECK IF NECESSARY)
genecat = p0;
for s = 1: size(pmd.Gene,3)
    genecat = cellfun(@(x,y)[x,y], genecat, pmd.Gene(:,:,s), ...
        'UniformOutput', false);
end
%   Cat genes to cells
cnames = cellfun(@(x,y)[x,'_',y], pmd.Cell(:,:,ci), ...
            repmat(genecat,[1,1,numel(ci)]), ...
            'UniformOutput', false);
%   Remove any non-word characters to make usable as names
cnames = regexprep(cnames, {'\W','^[\d_]*(\w)'}, {'','$1'});
%   Make a name for each unique catted string
[ucnames] = unique(cnames(vpr));
%   Disregard invalid names and warn
gn = cellfun(@isvarname,ucnames);
if any(~gn); warning(['Invalid cell name found for idx: ', ucnames{~gn}]);
    ucnames = ucnames(gn);
end
%   Assign matching xy positions to each compiled cell name
for s = 1:numel(ucnames)
    %   Assemble list of xys with the current cname
    idx.Cell.(ucnames{s}) = [pmd.xy{strcmp(ucnames{s},cnames)}];  
end

%FOR Txs (ptx, tx1, tx2, etc.)
%   Get treatment stage names (pTx, Tx1, Tx2, etc.)
txn = pfn( ~cellfun(@isempty,regexpi(pfn, '.*tx.*')) );   
for s = 1:numel(txn);
    %Make unique list of tx (drug) names in this Tx
    %   Get indices for drug names (may be multiple), and get names
    di = 1:5:size(pmd.(txn{s}),3);    dn = pmd.(txn{s})(:,:,di);
    %   Parse dosing and timing for each name set
    dsv = cat(3,str2double( pmd.(txn{s})(:,:,di+1) ));
    dun = cat(3,pmd.(txn{s})(:,:,di+2));
    tmv = cat(3,str2double( pmd.(txn{s})(:,:,di+3) ));
    tun = cat(3,pmd.(txn{s})(:,:,di+4));
    %   Remove NaNs
    bi = cellfun(@(x)any(isnan(x)), dn);   [dn{bi}] = deal('');
    %   Make names proper for fields
	dn = regexprep(dn, '\W', '');
    %   Get unique list, and remove any empties
    udn = unique( dn );  udn = udn(~cellfun(@isempty,udn));
    %Fill all values for each treatment 
    for ss = 1:numel(udn)
        %   Get well indices with matching name
        wi = strcmp(udn{ss},dn);  nwi = nnz(wi);
        %   List xys with that name
        idx.(txn{s}).(udn{ss}).xy = [pmd.xy{any(wi,3)}];  
        %   Get mapping of xys in wells
        xym = cellfun(@(x,y)ones(1,numel(x)).*y, pmd.xy(any(wi,3)), ...
            num2cell(1:nwi)', 'Un', 0);     xym = [xym{:}];
        %Convert any doses or times with different units to minimum unit
        %   List dose (same order) and State dose unit
        [idx.(txn{s}).(udn{ss}).dose, idx.(txn{s}).(udn{ss}).dunit] = ...
            simplify_units(dsv(wi), dun(wi));
        %   List times (same order) and State time unit
        [idx.(txn{s}).(udn{ss}).time, idx.(txn{s}).(udn{ss}).tunit] = ...
            simplify_units(tmv(wi), tun(wi));
        %   Apply xy-mapping for xy-oriented list of values
        idx.(txn{s}).(udn{ss}).dose = idx.(txn{s}).(udn{ss}).dose(xym)';
        idx.(txn{s}).(udn{ss}).time = idx.(txn{s}).(udn{ss}).time(xym)';
    end
end

end


%% SUBFUNCTIONS
% - Subfunction to check if any of a cell array are filled (not NaN) -
function isn = isbad(x)
isn = cellfun(@(xx)numel(xx)==1 && isnumeric(xx) && isnan(xx), x);
end


% ----- Subfunction to parse data inputs -----
function out = eparse(x, tt, uu)
%Ensure cell wrapping of input
if ~iscell(x); x = {x}; uncell = true; else uncell = false; end
x = cellfun(@num2str, x, 'UniformOutput', false);   %Ensure characters

%Catch alternative unit phrasings
x = regexprep(x, {'d(ay)?s?','h(ou)?rs?', 'mi?n(ute)?s?', ...
    'sec(ond)?s?', 'tps'}, {'d','h','m','s','tp'});

%Look for entry separation (via semi-colons)
tmp = regexpi(x, ';', 'split');

%  Lump all tx flags to same procedure
tt = regexprep(tt, '\w*tx\d*', 'tx', 'ignorecase');

out = cell(numel(tmp), 1);
for sw = 1:numel(tmp)
    %Look for recognized units and structure, depends on type (tt)
    switch lower(tt)
        case 'basic'
            y = regexpi(tmp{sw}, ...
                ['(?<val>[\d\.,]*)\s*(?<unit>',uu,')'], 'names');
            for s = 1:numel(y)
                if isempty(y{s}); out{sw}{s} = tmp{sw}{s};
                else out{sw}{s} = {str2double(y{s}.val), y{s}.unit}; end
            end
            out{sw} = [out{sw}{:}];
            
        case 'xy'
            xyt = regexpi(tmp{sw}, '(,|\s)(?=\s*\d)', 'split'); %Split on comma, spc
            %   Get numbers, excluding wells labeled as blank or background
            y = regexpi( [xyt{:}] , '\d+(?!.*b[agkl])', 'match' ); 
            %   Pack as numbers
            out{sw}{1} = cellfun(@str2double, [y{:}]);  
            
        case 'cell'
            y = regexpi(tmp{sw}, ['\s*(?<name>[^@]*)',...
                '@(?<val>[\d\.,]*)\s*(?<unit>',uu,')'], 'names');
            for s = 1:numel(y)
                if isempty(y{s}) %Pack empties with spacing
                    out{sw}((s-1)*3 + 1:3) = {tmp{sw}{s}, '', ''}; 
                else             %Cell gives: Type, Conc., Units
                out{sw}((s-1)*3 + 1:3) = {y{s}.name, y{s}.val, y{s}.unit};
                end
            end
            
        case 'gene'
            y = regexpi(tmp{sw}, '^\s*(\S.*\S)\s*$', 'tokens'); %Trim spcs
            y = [y{:}]; out{sw} = [y{:}]; %Cat all listed names
            
        case 'tx'
            y = regexpi(tmp{sw}, ['^\s*(?<pfx>[^\d@]*)\s*',... %Prefix
                '(?<cv>[\d\.]*)\s*',...             %Conc. value
                '(?<cu>(?(cv)',uu,'))\s*',...       %Conc. units
                '(?<nm>(?(cv)[^@]*))',...           %Tx Name
                '(@?)\s*(?<tv>(?(5)[\d\.]*))',...   %Timing value (if @)
                '\s*(?<tu>(?(tv)',uu,'?))'], 'names');       %Timing units
            for s = 1:numel(y)
                if ~isfield(y{s},'tv'); %Fill missing time fields
                    [y{s}.tv] = deal(''); [y{s}.tu] = deal(''); end
                if isempty(y{s})  %Pack empties with spacing (bad data)
                    out{sw}{s} = [tmp{sw}, cell(1,4)]; continue; end  
                %   Determine if names are in nm or pfx
                if isempty(y{s}(end).nm) && ~isempty(y{s}(end).pfx)
                    nm1st = true;       else nm1st = false; 
                end
                for sy = 1:numel(y{s})
                    si = (sy-1)*5 + 1;
                    if nm1st;     out{sw}{s}(si) = {y{s}(sy).pfx};
                    else        out{sw}{s}(si) = {y{s}(sy).nm};      end
                    out{sw}{s}(si+(1:4)) = {y{s}(sy).cv, y{s}(sy).cu,...
                        y{s}(sy).tv, y{s}(sy).tu};
                end
            end
            out{sw} = [out{sw}{:}]; %Cat 's' layer
    end
end

%Remove cell wrapping if only one output
if uncell; out = out{1}; end

end


% ----- Subfunction defining units -----
function [u, tu, cu, du] = defineunits(pfn)
%List units by type (use regexp notation as convenient)
time_units = {'d', 'h', 'm', 's', 'tp'};
conc_units = {'[fpnum]?g(/[umcd]?l)?', '[fpnum]M', '%'};
dens_units = {'c(ells?)?/w(ell)?', 'c(ells?)?/[um]l'};

%Time only
tu = cell2mat(cellfun(@(x)[x,'|'], time_units, ...
    'UniformOutput', false));    tu = tu(1:end-1);
%Concentration only
cu = cell2mat(cellfun(@(x)[x,'|'], conc_units, ...
    'UniformOutput', false));    cu = cu(1:end-1);
%Cell density only
du = cell2mat(cellfun(@(x)[x,'|'], dens_units, ...
    'UniformOutput', false));    du = du(1:end-1);

%Assign units according to plate fieldnames
pfn = regexprep(pfn, '(tx)\d+', '$1', 'ignorecase'); u = cell(numel(pfn),1);
for s = 1:numel(pfn)
    switch lower(pfn{s})
        case {'xy', 'geno'};    u{s} = [];
        case 'cell';            u{s} = ['(',du,')'];
        case {'ptx', 'tx'};     u{s} = ['(', cu, '|', tu,')'];
    end
end

%Package individual unit types
tu = ['(',tu,')'];      cu = ['(',cu,')'];      du = ['(',du,')'];
end


% ----- Subfunction for unit converstion -----
function [v,u] = simplify_units(v,u)
%Short circuits
if ischar(u) || numel(u) < 2; return; end   %For only one unit
if isequal(u{:}); u = u{1}; return; end     %For uniform units

%Define reference units, and associated conversions
tr = {'s','m','h','d'};             tc = [60, 60, 24];
ur = {'f','p','n','u','m','c','d','_'}; uc = [1e3,1e3,1e3,1e3,10,10,10];

%Split numerator/deminator as needed
tmp = regexpi(u, '/', 'split');
%   Parse units for scale and unit type
uf = @(x)regexpi(x, '(?<scl>[fpnumcd]?)(?<typ>\D+)', 'names');
uu = cellfun(uf, tmp, 'UniformOutput', false);

%Assert unit consistency
nn = cellfun(@numel, uu);  badu = false; 
if all(nn == nn(1)); nn = nn(1); 
    typs = cell(1,nn); scs = typs; istime = false(1,nn);
    for ss = 1:nn  %Check each unit set (numerator/denominator)
        typs{ss} = cellfun(@(x)[x{ss}.typ,''], uu, 'UniformOutput', false);
        scs{ss}  = cellfun(@(x)[x{ss}.scl,''], uu, 'UniformOutput', false);
        %   Assert that all units are same (or are time units)
        if isequal(typs{ss}{:});  typs{ss} = typs{ss}{1};
            fu = cellfun(@isempty,scs{ss}); [scs{ss}{fu}] = deal('_');
        elseif all( cellfun(@(x)any(strcmpi(x, tr)), typs{ss}) )
            %   Setup for time units
            scs{ss} = typs{ss}; typs{ss} = []; istime(ss) = true; 
        else badu = true;
        end
    end
else badu = true;
end
%IF inconsistent units, warn and return
if badu     
    warning('IMAN:DSHEET', ['Inconsistent units found.  ',...
        'Units in index stucture may not be accurate.']); return;
end

%Identify scale values
ubuild = [];
for ss = 1:nn
    if istime(ss); reff = tr; cnv = tc; else reff = ur; cnv = uc; end
    %   Convert scale character to rank
    rnk = cellfun(@(x)find(strcmpi(x, reff)), scs{ss});
    %   Get minimum rank (smallest unit scale)
    mrnk = min(rnk);
    %   Convert unit values
    if rem(ss,2) == 1   %IF numerator
        v = arrayfun(@(vv,rr)vv.*prod(cnv(mrnk:rr)), v, rnk-1);
    else                %IF demoninator
        v = arrayfun(@(vv,rr)vv./prod(cnv(mrnk:rr)), v, rnk-1);
        ubuild = [ubuild,'/']; %#ok<AGROW>
    end
    ubuild = [ubuild, reff{mrnk}, typs{ss}]; %#ok<AGROW>
end
%Assign built unit
u = ubuild;

end


