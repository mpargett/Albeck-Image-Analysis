%IMAN_UNMIX
%   Linear unmixing for image data with multiple channels.  Uses known
%   spectra to calculate the overlap among fluorophores with given filters,
%   defining a mixing weight matrix. Inverts the mixing weight matrix (as a
%   linear system), and applies it to each image pixel.  The resulting
%   image approximates a uniform Global Gain (which includes cs, a
%   convenient scaling) times the concentration of each fluorophore. 
%
%   [im, invA, cs, freti] = iman_unmix(im, MD, SP, [], FRET);
%       calculates the unmixing weight matrix and corrects the image, im,
%       using metadata MD and spectral parameters (SP, optional).  If SP is
%       not provided, it is collected from internally designated stored
%       files.  freti is the returned indices of any FRET channel names
%       provided (this may be passed to future iman_unmix calls,
%       eliminating the need to pass MD).
%
%   MD, metadata must be a structure containing the following fields (with
%       examples, all are typical of the iman procedure Global Metadata):
%       cam.Desc       =   'Zyla5';             %Description of camera
%       exp.Exposure   =   [800, 500, 800];     %Exposure times (ms)
%       exp.ExVolt     =   [15, 15, 30];        %Relative voltages used
%       exp.Filter     =   {'Filter_CFP', 'Filter_YFP', 'Filter_Cherry2'};
%       exp.FPhore     =   {'mTurq2', 'YPet', 'mCherry'};
%
%   SP is a structure containing spectra and properties of any named
%   fluorophores and filters.  See iq_getspectralpar for its generation.
%
%   FRET is a structure indicating channels that are FRET pairs
%       (single-molecule only).  The fieldnames of the structure should
%       indicate the FRET reporter name (not used here), and the value is a
%       2 element cell array with the channel names of the FRET pair.
%       e.g. FRET.EKAR = {'CFP','YFP'};
%
%   im = iman_unmix(im, MD, SP, invA, FRET);
%       applies the previously calculated unmixing weight matrix, invA, to
%       the image, im.  This is much faster than re-calculating the matrix
%       for each image to be corrected.  SP is not needed, and MD is only
%       needed if FRET channel names need to be resolved.

%FIXME:  Consider excitation chaining? e.g. excited CFP exciting GFP?

function [im, invA, cs, freti] = iman_unmix(im, MD, SP, invA, freti, mct)
%Version check provision
if strcmpi(im,'version'); im = 'v2.0'; return; end

%Ensure neglected channels do not show in MetaData for unmixing
if ~isempty(MD) && isfield(MD.exp, 'cind');     
    exn = {'Channel', 'Filter', 'FPhore', 'Exposure', 'ExVolt'};
    if isfield(MD.exp,'ExLine'); exn = [exn, {'ExLine', 'ExWL'}]; end
    for s = exn; MD.exp.(s{1}) = MD.exp.(s{1})(MD.exp.cind);  end
end

%% Parse FRET channel names to indices (if provided)
if ~exist('freti','var') || isempty(freti);   freti = {};
if ~isempty(MD) && isfield(MD.exp, 'FRET')
    if isstruct(MD.exp.FRET)
        fn = fieldnames(MD.exp.FRET);  freti = cell(1,numel(fn));
        for s = 1:numel(fn)
            for ss = 1:2;  freti{s}(ss) = find( ~cellfun(@isempty, ...
                regexpi( MD.exp.Channel, MD.exp.FRET.(fn{s}){ss} )) ); end
        end
        %   OR if indices are provided directly, just keep them
    elseif iscell(MD.exp.FRET) && isnumeric(MD.exp.FRET{1}); freti = MD.exp.FRET;
    end
end
end


%% Calculate unmixing weight matrix
%   Skip if pre-calculated matrix is provided
cs = [];
if ~exist('invA','var') || isempty(invA)
    nfp = numel(MD.exp.FPhore);  nft = numel(MD.exp.Filter);
    
    %Get spectral data for FPs and Filters
    if ~exist('sp','var') || isempty(SP); SP = iq_getspectralpar; end
    
    %Get camera QE (quantum efficiency) curve
    X.qe = iq_getcameraqe(MD.cam.Desc, SP.WaveLength);
    %Get light source spectrum
    LSS = iq_getlightsourcespectra(MD.exp.Light);
    lsn = fieldnames(LSS); lsn = lsn(~strcmpi(lsn,'WaveLength'));
    X.ls = LSS.(lsn{1}).spec;
    
    %Get weight matrix
    %   Integrate each FP spectrum with each Filter Set
    A = zeros(nft,nfp);    fpi = setdiff(1:nfp, [freti{:}]);
    for sft = 1:nft
        %Get aggregate spectrum and excitation power (for multi-line)
        if numel(lsn) > 1 && iscell(MD.exp.ExVolt)  %IF multi-line source
            ExVolt = max(MD.exp.ExVolt{sft});  %Max line power
            lst = arrayfun(@(x,y)LSS.([MD.exp.Light,num2str(x)]).spec...
                .*y./ExVolt, MD.exp.ExLine{sft}, MD.exp.ExVolt{sft}, ...
                'UniformOutput', false);  %Scales by relative voltage
            X.ls = sum(cat(1,lst{:}),1);  %Scaled sum of spectra
        else ExVolt = MD.exp.ExVolt(sft);
        end
        %Add excitation time and power to structure
        SP.(MD.exp.Filter{sft}).tim = MD.exp.Exposure(sft);  
        SP.(MD.exp.Filter{sft}).pow = ExVolt;
        for sfp = fpi
            A(sft,sfp) = rel_imaging_gain( SP.(MD.exp.FPhore{sfp}), ...
                SP.(MD.exp.Filter{sft}), X );
        end
        %   For any FRET pairs, provide FP data simultaneously
        for s = 1:numel(freti)
           A(sft, freti{s}) = rel_imaging_gain( ...
               [ SP.(MD.exp.FPhore{freti{s}(1)}), ...
                 SP.(MD.exp.FPhore{freti{s}(2)}) ], ...
               SP.(MD.exp.Filter{sft}), X );
        end
    end

    %Look for mixing calibration table, and apply if provided
    if exist('mct', 'var')
        %Align indices from Calibration Table
        mct_fti = cellfun(@(x)mct.ftn.(x), MD.exp.Filter);
        mct_fpi = cellfun(@(x)mct.fpn.(x), MD.exp.FPhore);
        %Adjust Gain Table (A) by calibration ratios
        A = A.*mct.table(mct_fti,mct_fpi);
        %Adjust for any Direct Replacement values
        if isfield(mct, 'replace');   repi = ~isnan(mct.replace);
            repv = bsxfun( @times, mct.replace, A(sub2ind( ...
                [nft,nfp], mct.repref(mct_fpi), 1:nfp )) );
            A(repi) = repv(repi);
        end
    end
        
    %Truncate weight matrix, at 0.1% Threshold, per FP
    A(bsxfun(@rdivide, A, max(A, [], 2)) < 0.001) = 0;
    
    %Get norm of weight matrix, to scale for convenience
    cs = norm(A);
    
    %Invert weight matrix to give unmixing matrix
    invA = pinv(A./cs);
end


%% Unmix
%Check if an image is provided to unmix
if isempty(im); im = A; return; end

[nY, nX, nChan] = size(im);
%Unmix each pixel independently
%   Rearrange pixels to a matrix (nChan x nPixels)
im = reshape(shiftdim(im,2), nChan, nY*nX);
%   Multiply by inverted weight matrix
im = invA*im;
%   Rearrange unmixed values to image (nY x nX x nChan)
im = shiftdim( reshape(im, nChan, nY, nX), 1);

%If a FRET pair was used, calculate the FRET image (EfA)
for s = 1:numel(freti)
    im(:,:, freti{s}(1)) = im(:,:, freti{s}(1)) ./ im(:,:, freti{s}(2));
end

end


%% SUBFUNCTION:  Calculate relative imaging gain from spectra
function rig = rel_imaging_gain(fp, ft, z)
%Evaluate the Fluorescence Imaging Model to get the imaging gain (ig)
%   ig = sum(exF*FPab)*FPmec*FPqy*sum(FPem/sum(FPem)*emF*cQE)*tE*pE*GG;
%   Global Gain (GG) comprises camera Gain, geometric effects, absolute
%       light source power function, sample path length, etc.
%   Relative image gain (rig) is scaled by the (unknown) global gain, ig/GG

switch numel(fp)
    case 1 %Single fluorophore - standard imaging model
        rig = sum(z.ls .* ft.ex .* fp.ab) * fp.mec * fp.QY * ...
            (sum(fp.em .* ft.em .* z.qe)/sum(fp.em)) * ft.tim * ft.pow;
        
    case 2 %Two fluorophores - FRET imaging model        
        %   Donor absorbance
        ad = sum(z.ls .* ft.ex .* fp(1).ab) * fp(1).mec * ft.tim * ft.pow;
        %   Acceptor absorbance
        aa = sum(z.ls .* ft.ex .* fp(2).ab) * fp(2).mec * ft.tim * ft.pow;
        %   Donor emission
        ed = (sum(fp(1).em .* ft.em .* z.qe)/sum(fp(1).em)) * fp(1).QY;
        %   Acceptor emission
        ea = (sum(fp(2).em .* ft.em .* z.qe)/sum(fp(2).em)) * fp(2).QY;
        %Relative intensity gains with [c*EfA, c]'
        rig = [ad*(ea - ed), ad*ed + aa*ea];
end
    
end