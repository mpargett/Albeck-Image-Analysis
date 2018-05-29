%IMAN_POWERRATIO_EST
%   Estimate the power ratio between two channels (i.e. for FRET pairs).
%
%   prat = iman_powerratio_est(MD, SP, CNAMES)
%       returns the power ratio using a Global Metadata structure MD and
%       Spectral Data structure SP, for the channel names specified in
%       CNAMES.
%
%   CNAMES must be a 2 element cell, with the channel name of the numerator
%   first, then denominator.  Names must match those in MD.  e.g. CNAMES =
%   {'CFP','YFP'};  (This is the default)
%
%   If SP is not provided, it is recovered (if possible) from data tables
%   via IQ_GETSPECTRALPAR.

function prat = iman_powerratio_est(MD, SP, cn)

%Default Fret pair = CFP/YFP
if ~exist('cn','var') || isempty(cn); cn = {'cfp','yfp'}; end

%Ensure spectral parameters are available
if ~exist('SP', 'var') || isempty(SP)
    %Run routine to extract spectral parameters
    SP = iq_getspectralpar;
end

%Assert propernames for MetaData on FPhores and Filters
[MD.exp.FPhore, MD.exp.Filter] = ...
    iman_naming('match',MD.exp.FPhore, MD.exp.Filter);

%Get Light source spectral power
lsPower = iq_getlightsourcespectra(MD.exp.Light);
lsn = fieldnames(lsPower); lsn = lsn(~strcmpi(lsn,'WaveLength'));
%   IF Multi Line source, assemble line list
if numel(lsn) > 1 && isfield(MD.exp, 'ExLine') && ~isempty(MD.exp.ExLine)
    %Store index of lines relevant to each field
    lines = cellfun(@str2double, regexpi(lsn, '(?<=\D*)\d+', 'match'));
    is_multiline = true;
    %Double check MetaData format
    if ~iscell(MD.exp.ExVolt); MD.exp.ExVolt = num2cell(MD.exp.ExVolt); end
else %IF single line, store spectrum and ensure numeric array MD
    is_multiline = false;  lsp = lsPower.(lsn{1}).spec;
    if iscell(MD.exp.ExVolt); MD.exp.ExVolt = [MD.exp.ExVolt{:}]; end
    if iscell(MD.exp.Exposure); MD.exp.Exposure = [MD.exp.Exposure{:}]; end
end

%Get camera QE (quantum efficiency) curve
QEcam = iq_getcameraqe(MD.cam.Desc, SP.WaveLength);


%Evaluate channel data to get intensity gains
rig = [0,0];
for s = 1:2
    %Get Channel ID for the exposure
    id = find( ~cellfun('isempty', regexpi(MD.exp.Channel, cn{s}, 'start')) );
    
    %   Get light source spectra if a Multi-Line source
    if is_multiline;  lsp = zeros(size(lsPower.WaveLength));  exv = 1;
        cln = find(ismember(lines, MD.exp.ExLine{id}));
        for ss = 1:numel(cln) %Each Line can be scaled differently
            lsp = lsp + lsPower.(lsn{cln(ss)}).spec .* MD.exp.ExVolt{id}(ss); 
        end
    else exv = MD.exp.ExVolt(id); %IF single Line, voltage is a vector
    end
    
    %Define relevant filter and fluorophore spectra
    fp = SP.(MD.exp.FPhore{id});  ft = SP.(MD.exp.Filter{id});
    %Calculate relative intesity gain
    rig(s) = sum( lsp .* ft.ex .* fp.ab ) * sum( ft.em .* fp.em .* QEcam ) ...
        * fp.mec * fp.QY * MD.exp.Exposure(id) * exv / sum(fp.em);
end

%Take ratio of RIG to get Power Ratio (e.g. Donor to Acceptor in FRET pair)
prat = rig(1) / rig(2);

%Power Ratio Calculation in brief
% prat = Abs(D_d) * QY(D) * Em(D_d) * t(d) * P(d) ./ ...
%        Abs(A_a) * QY(A) * Em(A_a) * t(a) * P(a);
% 

end