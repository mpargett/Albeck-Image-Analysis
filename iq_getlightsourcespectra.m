%IQ_GETLIGHTSOURCESPECTRA
%   Gathers spectral data for microscopy light source(s), as assembled in
%       datasheets.
%
%   q = iq_getlightsourcespectra()
%       returns the structure q with all light source spectra

function [q, cc] = iq_getlightsourcespectra(lsn)
%Define standard light source names
lsnames = {	'SOLA',         '^.*sola(?!.*(se|sm|ii)).*$';...
            'SOLA_SM_II',   '^.*sola(\s|_)?sm(\s|_)?(2|ii)?';...
            'SOLA_SE_II',   '^.*sola(\s|_)?se(\s|_)?(2|ii)?';...
            'SPECTRAX',     '^.*spectra\s?x';                      };

%Check for input of light source name, or metadata
gmd = []; cc =[];
if exist('lsn','var') && ~isempty(lsn)
    %Check for a Metadata structure and extract Light Source Name
    if isstruct(lsn) && isfield(lsn, 'exp')
        gmd = lsn; lsn = gmd.exp.Light;
        %   Get calibration data
        ct = iq_getexcal();
    end
    %   Sanitize source name
    lsn = regexprep(lsn, lsnames(:,2), lsnames(:,1), 'ignorecase');
    if iscell(lsn) && numel(lsn) == 1; lsn = lsn{1}; end %Check input
    %   Restrict info to the relevent Light Source
    nuse = ~cellfun(@isempty, regexpi(lsn, lsnames(:,2)));
    lsnames = lsnames(nuse,:);
end

%Wavelength range
wl_lo = 300;  wl_hi = 800;      %nm
q.WaveLength = wl_lo:wl_hi;   %Expect always 1 nm step

%Spectra X Line Definitions {Light, Filter, Calibration Power, ExWL}
lnn.SPECTRAX.FName = 'sx(\s|_)filt(er)?';
lnn.SPECTRAX.Lines(1,:) = {'Violet',  'Violet', 21.6,    395}; %266
lnn.SPECTRAX.Lines(2,:) = {'Blue',    'Blue',  	22.6,    440}; %224
lnn.SPECTRAX.Lines(3,:) = {'Cyan',    'Cyan', 	17.8,    470}; %149
lnn.SPECTRAX.Lines(4,:) = {'Teal',    'Teal',  	1.9,     508}; %62
lnn.SPECTRAX.Lines(5,:) = {'Green',   'Green', 	3.95,    555}; %225
lnn.SPECTRAX.Lines(6,:) = {'Red',     'Red',   	11.6,    640}; %224
%lnn.SpectraX{???} = {'Green', 'Yellow', 271, ???};

%Define data files to use
data_light = 'LightSourceSpectra.xlsx';

%Get data from file
[d.n, d.t, d.r] = xlsread(data_light);


%% PARSE FOR HEADERS
%Main header index for data in properties file
d.mainhead = find( ~cellfun('isempty', ...
    regexpi(d.t(:,1), '^model')) );


%% COLLECT SPECTRA
%Get names used in the properties data file
names_in = d.t(d.mainhead,:);
f_ind = cellfun(...                     %Get sheet index for each name
    @(x)find( ~cellfun('isempty', regexpi(names_in, x)) ),...
        lsnames(:,2), 'UniformOutput', false);
    f_use = find( cellfun(@(x)any(x), f_ind) )';     %Remove unmatched names
fmatch = cell2struct(f_ind(f_use), lsnames(f_use,1), 1); %Structure with index
fmn = fieldnames(fmatch);

%Get indices for the wavelengths (300 to 800 nm, incremented by 1)
ind_u = d.mainhead + find(~cellfun(...
    'isempty',regexpi(d.t(d.mainhead+1:end,1), 'wave(length)?') ));
ind_wl = ind_u + ( find([d.r{ind_u+1:end,1}] == wl_lo) : ...
                    find([d.r{ind_u+1:end,1}] == wl_hi) );

%Get names of any parameters (between Model and WaveLength)
ind_p = d.mainhead + 1 : ind_u - 1;
pn = regexprep( d.r(ind_p, 1), {'\s','\W'}, {'_',''});

%IF a MultiLine light source requested, parse Line definitions
mln = fieldnames(lnn);  
isml = cellfun(@(x)strcmpi(fmn(:), x), mln, 'UniformOutput', false);
isml = cat(2,isml{:});  mln = mln( any(isml,1) );

%Get spectral data for Single Line devices
for s = find(~isml)'
    %Include any parameters in the data table
    for ss = 1:numel(pn)
        q.(fmn{s}).(pn{ss}) = d.r{ind_p(ss), fmatch.(fmn{s})};
    end
    q.(fmn{s}).spec = [d.r{ind_wl, fmatch.(fmn{s})}];
    q.(fmn{s}).unit = d.r{ind_u, fmatch.(fmn{s})};
end

%Get spectral data for Multi Line devices
for sfn = 1:numel(mln)   %FOR each MultiLine device
    %Assemble excitation source filter spectra separately
    %   Get filter names available
    ftn = regexpi(names_in, [lnn.(mln{sfn}).FName,...
        '[_\s](?<name>[^_\d\s]*$)'], 'names');
    %   Get indices and store names (of both filter family and part)
    f_i = find(~cellfun('isempty', ftn)); 
    ftn = squeeze(struct2cell([ftn{f_i}]));
    %   Include a null filter
    q.(mln{sfn}).Filter.None = ones(size(ind_wl));
    for sl = 1:numel(f_i) %For each filter
        %Collect filter spectrum 
        q.(mln{sfn}).Filter.(ftn{sl}) = [d.r{ind_wl,f_i(sl)}];
    end
    
    %Parse light sources and filters for the device
    nl = size(lnn.(mln{sfn}).Lines, 1);
    for sl = 1:nl       %FOR each Line
        %Get source index
        s_i = fmatch.(mln{sfn})( ~cellfun(@isempty, ...
            regexpi(names_in(fmatch.(mln{sfn})), ...
            [mln{sfn},'.*',lnn.(mln{sfn}).Lines{sl,1}])) );
        %Get filter index
        f_i = ~cellfun(@isempty, regexpi(names_in, ...
            [lnn.(mln{sfn}).FName,'.*',lnn.(mln{sfn}).Lines{sl,2}]));
        %Store parameters
        for ss = 1:numel(pn)
            q.(mln{sfn}).Line(sl).(pn{ss}) = d.r{ind_p(ss), f_i};
        end
                
        %Calculate filtered spectrum (mW/nm)
        tspec = [d.r{ind_wl,s_i}] .* [d.r{ind_wl,f_i}];
        %Scale spectrum to calibration power and assign
        q.(mln{sfn}).Line(sl).spec = lnn.(mln{sfn}).Lines{sl,3} * ...
            tspec./trapz(tspec);
        %Include units
        q.(mln{sfn}).Line(sl).unit = d.r{ind_u, s_i};
        %Provide raw spectra
        q.(mln{sfn}).Line(sl).rawspec = [d.r{ind_wl,s_i}];
    end
    
end

%Remove all NaNs, replacing with zeros
q = zeronan(q);


%% Calibrated Channels
%   Run if Metadata for an experiment is provided
if ~isempty(gmd)
    %Number of photons (n) may be calculated by the Planck-Einstein relation:
    %   E = n*hc/l, where h = Planck's const(J-s), c = speed of light(m/s),
    %       l = wavelength(m), n = number of photons, and E = energy(J)
    %       and Power (pw) = E/s;
    %   Pre-define relevant constants (includes 1e12 for nm-m and mW-W)
    hc = 6.62606957*10^-34 * 299792458     * 1e27;
    %    ^Planck's           ^Speed of Light ^Extra scaling for convenience
    
    %Determine Light Source and get relevant indices
    lsi = find(~cellfun(@isempty, regexpi({ct.chan.model}, lsn)));
    
    %Define Relative Flux Estimate function
    rfe = @(ci,x)interp1(ct.pwr, ct.chan(ci).pwr, x, 'linear')...
        .*ct.chan(ci).calwl/hc; %Scaled by avg. wavelength for flux
    %Get number of channels and voltage values (to match type)
    nc = numel(gmd.exp.Channel);   vlt = gmd.exp.ExVolt;
   
    %Load Filter Spectral Data, as needed
    if ~exist('sd','var') || isempty(sd) %#ok<NODEF>
        sd = iq_getspectralpar(); end
    
    %Alter protocol if Light Source is Single vs. Multi-Line
    switch lower(lsn)
        case {'sola'}; pflux = nan(nc,1);         	%Single Line
            %   Ensure array-packed voltage %
            if iscell(vlt); vlt = cell2mat(vlt); end
            %Estimate relative photon flux
            for s = 1:nc
                %Determine Filter and set Channel Index
                ci = lsi(~cellfun(@isempty, ...
                    regexpi(gmd.exp.Filter{s}, {ct.chan(lsi).filter})));
                %Calculate Relative Flux Estimate
                pflux(s) = rfe(ci, vlt(s));
                %Divide by filter-cross-source power (to match expected usage)
                pflux(s) = pflux(s)./...
                    sum(sd.(gmd.exp.Filter{s}).ex.*q.(lsn).spec);
            end; cc = pflux; 
            
        case {'spectrax'}; cc = cell(nc,1);      	%Multi-Line
            %   Ensure cell-encapsulated voltage %
            if ~iscell(vlt); vlt = num2cell(vlt); end
            fltn = cellfun(@(x)regexprep(x,'Filter_',''), ...
                gmd.exp.Filter, 'Un', 0);
            %Estimate relative photon flux (per excitation line)
            for s = 1:nc
                %Find valid calibrations for the used Filter
                uflt = lsi(ismember({ct.chan(lsi).filter}, fltn{s})); 
                %Get line IDs and match to Calibration Lines
                [ci, loc] = ismember([ct.chan(uflt).exline], gmd.exp.ExLine{s});
                %   Check that all lines are calibrated
                if isempty(loc) || ~all(loc)
                    error(['Calibrations are not available for ',...
                        'all Excitation Line and Filter combinations']); 
                end
                ci = uflt(ci);  ci = ci(find(loc));      %#ok<FNDSB>
                li = ct.chan(ci).exline;
                %Calculate Relative Flux Estimate, per Line
                pflux = arrayfun(rfe, ci, vlt{s}, 'Un', 0);
                
                %Divide by filter-cross-source power (to match expected usage)
                pflux = cellfun(@(x,L)x./sum(sd.(gmd.exp.Filter{s}).ex.*L./trapz(L)),...
                    pflux, {q.(lsn).Line(li).spec}, 'Un', 0);
                
                %Assemble calibrated channel
                cc{s} = cellfun(@(c,p,f)p.*c./trapz(c), ...
                    {q.SPECTRAX.Line(li).spec}, pflux, 'Un', 0);
                cc{s} = sum(cat(1,cc{s}{:}),1); %Sum sources     
            end            
    end
    
end



end


%% Sub-Function to set NaN to zero
function out = zeronan(in)
switch class(in)
    case 'struct'
        if numel(in) == 1
            out = structfun(@zeronan, in, 'UniformOutput', false);
        else tname = fieldnames(in); tcell = struct2cell(in);
            out = cell2struct(zeronan(tcell), tname, 1);
        end
    case 'cell';   out = cellfun(@zeronan, in, 'UniformOutput', false);
    case 'double'; out = in; out(isnan(in)) = 0;
    otherwise;     out = in; 
end
end

%% Sub-Function: Get Excitation Calibrations from table
function ct = iq_getexcal()
%Load calibration table data (from Excel sheet)
[d.n, d.t, d.r] = xlsread('LightSourceCalibration.xlsx');

%Find calibration header
%Main header index for data in calibration file
d.mainhead = find( ~cellfun('isempty', ...
    regexpi(d.t(:,1), '^mfr\.$')) );
%Find number of headers
d.numhead = find( cellfun(@(x)isequal(x,1),  d.r(:,1)) ) - d.mainhead;
%Find number of entries
d.lastchan = find(cellfun(@(x)isempty(x) || (isnumeric(x)&&isnan(x)), ...
    d.r(d.mainhead,:)), 1, 'first') -1;

%Collect Metadata fields
md = lower( d.t(d.mainhead + (0:d.numhead-1),1) );
%Sanitize Metadata field names
md = regexprep(md, '\W', '');
%   Populate Metadata
mdd = d.r(d.mainhead + (0:d.numhead-1), 2:d.lastchan);
vf = ~any(cellfun(@(x)isnumeric(x)&&isnan(x),mdd),1);
%   Include column index in sheet
md = cell2struct([mdd(:,vf);num2cell(find(vf)+1)], [md;{'col'}], 1);

%Get calibration values for all channels
vals = cell2mat(d.r(d.mainhead+d.numhead:end, [1,md.col]));
%   Clear NaN rows and cell-pack
vals = num2cell(vals(~any(isnan(vals),2),:),1);
%   Assign to MD structure
[md.pwr] = deal(vals{2:end});
%Handle the % power list somehow... ***
ct.pwr = vals{1};
ct.chan = md;    

end


% function mr = get_Cy5Mirror()
% %Original call for this:
% %             %Correct for Cy5 Mirror used in calibration
% %             mr = get_Cy5Mirror();	mrr = 1-mr.Cy5Mirror; %Get spectra
% %             sxscl = arrayfun(@(x)sum(x.rawspec)./sum(x.rawspec.*mrr), ...
% %                 q.SPECTRAX.Line); %Get scaling factor
% %             for li = 1:numel(sxscl) %Loop and scale power values
% %                 cti = cellfun(@(x)isnumeric(x) && x==li, {ct.chan.exline});
% %                 ct.chan(cti).pwr = ct.chan(cti).pwr.*sxscl(li);
% %             end        
%         
% %Load sheet data
% [~,t,r] = xlsread('FilterSpectra.xlsx');
% %Get header indices
% hi = find(~cellfun('isempty',regexpi(t(:,1), 'application')));
% li = find(~cellfun('isempty',regexpi(t(:,1), 'wavelength')));  
% %   Assemble header structure
% h = cell2struct( num2cell(hi:li-1)', t(hi:li-1,1), 1);
% %   Get number of rows
% nr = find(cellfun(@(x)isnumeric(x) && ~isnan(x), r(li+1:end,1)), 1, 'last');
% %Assign data
% wl = cell2mat(r(li+(1:nr),1));
% ri = li + find(wl == 300);
% rf = li + find(wl == 800);
% 
% 
% %Get index for Cy5 Mirror
% cyi = find(~cellfun('isempty',regexpi(t(h.Application,:), 'Cy5')) & ...
%           ~cellfun('isempty',regexpi(t(h.Type,:), 'Mirror'))        );
% %Assign data
% mr.WaveLength = wl((ri-li):(rf-li))';
% mr.Cy5Mirror = cell2mat(r(ri:rf,cyi))'; %#ok<FNDSB>
% end

%More Obsolete:
% %   Get source filter if specified (or use default)
%                 if isfield(gmd.exp, 'ExFilt')
%                     %   Define a lenient filter name match (first 3 chars)
%                     ftm = cellfun(@(x)[x(1:3),'\w*$'], ftn, 'Un', 0);
%                     %   Match filter names per channel
%                     sflt = cellfun(@(x) ftn{~cellfun('isempty', ...
%                         regexpi(x, ftm))}, gmd.exp.ExFilt{s}, 'Un', 0);
%                 else    sflt = lnn.SPECTRAX.Lines(li,2)'; %Defaults
%                 end
%                 %Calculate calibrated channel spectrum, including filter
%                 cc{s} = cellfun(@(c,p,f)p.*c.*q.SPECTRAX.Filter.(f)./trapz(c), ...
%                     {q.SPECTRAX.Line(li).rawspec}, num2cell(pflux), sflt, 'Un', 0);
%                 cc{s} = sum(cat(1,cc{s}{:}),1); %Sum sources

