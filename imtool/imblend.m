function J = imblend(I, ix)
% function J = imblend(I, ix)

if numel(I) == 1, J = I{1}; return; end
cast_ = @(x) cast(x, 'like', I{1});
I = cellfun(@single, I, 'UniformOutput', 0);
masks = cellfun(@create_mask, ix, 'UniformOutput', 0);
mask = sum(cat(3, masks{:}), 3);
masks = cellfun(@(x) x./mask, masks, 'UniformOutput', 0);
J = cellfun(@(x,y) bsxfun(@times,x,y), masks, I, 'UniformOutput', 0);
J = sum(cat(4, J{:}), 4);
J = cast_(J);

function mask = create_mask(ix)
mask = bwboundary(ix);
mask = bwdist(mask, 'city');
mask = mask.*ix;
mask = mask/max(max(mask));

function b = bwboundary(ix)
se = strel('disk', 1);
b = xor(imdilate(ix, se), ix);
