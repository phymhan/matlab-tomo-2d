function im = imblendmask(im, mask, color, alpha, outType)
% function im = imblendmask(im, mask, color, alpha)
%   alpha: layer(mask) = layer(mask)*(1-alpha)+color(k)*alpha;

thisType = class(im);
if ~exist('outType', 'var')
    outType = thisType; % or 'uint8' for images
elseif strcmp(outType, 'same') || isempty(outType)
    outType = thisType;
end
if ~exist('alpha', 'var'), alpha = 0.5; end
if size(im, 3) == 1, im = cat(3, im, im, im); end %###
if max(color) <= 1, color = color*255; end
color = color.'; % one row defines one color
for c = 1:size(im, 3)
    x = double(im(:,:,c));
    x(mask) = alpha(1).*color(c) + (1-alpha(1)).*x(mask);
    im(:,:,c) = x;
end
if numel(color) > 3
    % add boundaries
    alpha = [alpha 1];
    color = color(4:end);
    bd = edge(mask);
    for c = 1:size(im, 3)
        x = double(im(:,:,c));
        x(bd) = alpha(2).*color(c) + (1-alpha(2)).*x(bd);
        im(:,:,c) = x;
    end
end
if ~isa(thisType, outType), im = cast(im, outType); end
