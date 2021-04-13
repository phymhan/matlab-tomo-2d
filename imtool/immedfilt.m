function imOut = immedfilt(im, siz, padopt, outType)
% function immedfilt(im, siz, padopt, outType)

if ~exist('outType', 'var') || isempty(outType)
    outType = class(im);
end
if ~exist('padopt', 'var') || isempty(padopt)
    padopt = 'symmetric';
end
if ~exist('siz', 'var') || isempty(siz)
    siz = [3 3];
end
[height, width, nChnl] = size(im);
imOut = zeros(height, width, nChnl);
for iChnl = 1:nChnl
    imOut(:,:,iChnl) = medfilt2(im(:,:,iChnl), siz, padopt);
end
imOut = cast(im, outType);
