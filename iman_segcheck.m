%IMAN_SEGCHECK
%   Check segmentation on a single image.
%   
%   iman_segcheck(ip,op,...)
%       uses the ip and op parameters set in the celltracer run script to
%   print an example image from the dataset with segmented masks, per
%   current parameters.  It is recommended to run the celltracer script 
%   through calling the partner routine (e.g. celltracer_v2_partner) prior
%   to running iman_segcheck.  iman_segcheck takes additional parameters a
%   Name/Value pairs, options for which are defined below.
%
%   pastinfo = iman_segcheck(ip,op,...)
%       also returns a structure containing the Data Access Object,
%       MetaData and Background information necessary to speed up further
%       calls for the same dataset.  Provide the information as the
%       "pastinfo" parameter, as in:
%   pastinfo = iman_segcheck(ip, op, 'pastinfo', pastinfo)
%
%Optional Parameters:
%   t           - Time point to print
%   c           - Channel to plot
%   xy          - XY Position to plot
%   z           - Z Position to plot (typ. not used)
%   nclr        - Color of nuclear masks, RGB triplet
%   cclr        - Color of cytoplasm masks, RGB triplet
%   pastinfo    - Data Access info for future calls on the same dataset
%   imthresh    - Image percentile threshold for display (default 95)
%   display     - TRUE (default) to display image and masks

function [imd, imo] = iman_segcheck(ip,op,varargin)
%% Handle options
p.nclr = [0,0,1];
p.cclr = [0,1,0.5];
p.pastinfo = [];
p.imthresh = 95;
p.t  = 1;
p.c  = op.seg.chan;
p.xy = 1;
p.z  = 1;
p.display = true;

%Input option pair parsing:
nin = length(varargin);     %Check for even number of add'l inputs
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end
%Splits pairs to a structure
for s = 1:2:nin;   p.(lower(varargin{s})) = varargin{s+1};   end

%% Get Image
%Access image data
if isempty(p.pastinfo);    
    imd.dao = iman_imageaccess(ip.fname, ip.indsz);
    imd.GMD = iman_getmeta(imd.dao, ip.bkmd, op);
    imd.bkg = cell(ip.indsz.xy,1);
else    imd = p.pastinfo;
end

%Pull image
%   Initialize image
im = nan(imd.GMD.cam.PixNumY, imd.GMD.cam.PixNumX, ip.indsz.c);
%   Load image data from DAO
for sc = 1:ip.indsz.c;   
    im(:,:,sc) = iman_getframe(imd.dao, [p.t,sc,p.xy,p.z]);   
end

%Refine image
[im, e_inv] = iman_refine(im, imd.GMD, ip.bval, op.objbias);

%% Handle background
txyz = [p.t, p.xy, p.z];  %By default, point to this image
%   If alternate XY indicated, use it
if ~isempty(ip.bkg(txyz(2)).altxy);  txyz(2) = ip.bkg(txyz(2)).altxy; end
%   If not dynamic, use time point 1
if ip.bkg(txyz(2)).dyn; nbt = ip.indsz.t; else txyz(1) = 1; nbt = 1; end

%IF this XY not yet accessed for background, pre-allocate w/ NaN
if isempty(imd.bkg{txyz(2)}); imd.bkg{txyz(2)} = nan(nbt, ip.indsz.c); end
%Check if this background was already collected
if any(isnan(imd.bkg{txyz(2)}(txyz(1),:)))  %IF not collected, do so now
    if ip.bkg(txyz(2)).fix;  bkg(txyz(1),:) = ip.bkg(txyz(2)).reg - ip.bval;  %IF fixed, apply
    else 	bkg(1,:) = get_bkg(imd, ip, txyz, e_inv);  %ELSE get values
    end;    imd.bkg{txyz(2)}(txyz(1),:) = bkg;
else bkg = imd.bkg{txyz(2)}(txyz(1),:);	%IF already collected, use old
end

%Remove backround level from image
im = bsxfun( @minus, im, reshape(bkg, 1, 1, numel(bkg)) );

%Perform linear channel unmixing (optional)
if op.unmix; im(:,:,op.cind) = iman_unmix(im(:,:,op.cind),imd.GMD); end
                
%% Segment and mask image
%   Adjust size parameters for calibration
pxscl = (imd.GMD.cam.PixSizeX + imd.GMD.cam.PixSizeY)./2;
op.seg.minD = op.seg.minD./pxscl;  op.seg.maxD = op.seg.maxD./pxscl;
%Get segmentation masks
[m, nmask] = iman_cellid(im, op.seg, bkg);
[valcube, mask] = iman_cellmask(im, m, op, nmask, bkg);

%Display
if p.display
for s = p.c(:)'
    imm = iman_maskoverlay(im(:,:,s),mask.nuc,p.nclr,3, p.imthresh);
    imm = iman_maskoverlay(imm,mask.cyt,p.cclr,2, p.imthresh); 
    figure; imshow(imm); 
    title(['XY:', num2str(p.xy),' C:', num2str(s), ' T:', num2str(p.t)]);
end
end

%IF full image and data are requested, package
if nargout > 1;    imo.im = im;    imo.m = m;    imo.msk = mask;    end

end


function [bkg] = get_bkg(imd, ip, txyz, e_inv)

%Pull image for background
im = zeros(imd.GMD.cam.PixNumY, imd.GMD.cam.PixNumX, ip.indsz.c);
for sc = 1:ip.indsz.c;
    im(:,:,sc) = iman_getframe(imd.dao, [txyz(1), sc, txyz(2), txyz(3)]);
end

%Refine image
im = iman_refine(im, imd.GMD, ip.bval, e_inv);

%Get background value from image region
bkg = mean(mean( ...
    im( ip.bkg(txyz(2)).reg(3):ip.bkg(txyz(2)).reg(4), ...
    ip.bkg(txyz(2)).reg(1):ip.bkg(txyz(2)).reg(2), : ), 1 ), 2);
end
