%IQ_GETLIGHTSOURCESPECTRA
%   Gathers spectral data for microscopy light source(s), as assembled in
%       datasheets.
%
%   q = iq_getlightsourcespectra()
%       returns the structure q with all light source spectra

function q = iq_getlightsourcespectra(lsn)
%Light source names
lsnames = { ...
    'SOLA',         'sola(?!.*(se|sm|ii)).*$';...
    'SOLA_SM_II',   'sola(\s|_)?sm(\s|_)?(2|ii)?';...
    'SOLA_SE_II',   'sola(\s|_)?se(\s|_)?(2|ii)?';...
    'SPECTRAX',     'spectra\s?x';...
};

%If a desired source is requested, remove all others from list
if exist('lsn','var') && ~isempty(lsn)
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
names_in = lower( d.t(d.mainhead,:) );
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
            q.([mln{sfn},num2str(sl)]).(pn{ss}) = d.r{ind_p(ss), f_i};
        end
        
        %Calculate filtered spectrum (mW/nm)
        tspec = [d.r{ind_wl,s_i}] .* [d.r{ind_wl,f_i}];
        %Scale spectrum to calibration power and assign
        q.([mln{sfn},num2str(sl)]).spec = lnn.(mln{sfn}).Lines{sl,3} * ...
            tspec./trapz(tspec);
        %Include units
        q.([mln{sfn},num2str(sl)]).unit = d.r{ind_u, s_i};
    end
end

%Remove all NaNs, replacing with zeros
q = zeronan(q);


end


%% Sub-Function to set NaN to zero
function out = zeronan(in)
switch class(in)
    case 'struct'; out = structfun(@zeronan, in, 'UniformOutput', false);
    case 'cell';   out = cellfun(@zeronan, in, 'UniformOutput', false);
    case 'double'; out = in; out(isnan(in)) = 0;
    otherwise;     out = in; 
end
end