%MaskOverlay

function im = iman_maskoverlay(im,msk,mcolor,sz,imth)
%Replicate grayscale to RGB
if size(im,3) == 1
    im = im./prctile(im(:),imth); %Scale for visisbility
    im = repmat(im,1,1,3); 
end

%Unpack msk if necessary
if isa(msk,'uint32'); msk = bwunpack(msk,size(im,1)); end

%   Dilate and subtract mask for a boundary
msk = imdilate(msk,strel('disk',sz)) & ~msk;
%   Set mask boundaries to desired color
im(repmat(msk,1,1,3)) = ones(nnz(msk),1)*mcolor;

end