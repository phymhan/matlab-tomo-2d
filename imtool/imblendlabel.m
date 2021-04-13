function Ic = imblendlabel(I, L, C, alpha)
% I, original image
% L, segmentation result
% C, color code where i-th row defines color for i-th label
% alpha, Ic = (1-alpha)*I + alpha*L
% Ic, output image with colorized label

if ~exist('alpha', 'var') || isempty(alpha), alpha = 0.8; end
sz = [size(I,1) size(I,2)];
if size(I,3) == 1, I = repmat(I, [1 1 3]); end
L = imresize(L, sz, 'nearest');
Lc = zeros([sz 3]);
for l = setdiff(unique(L(:))', 0)
    mask = L==l;
    for j = 1:3
        Lc(:,:,j) = Lc(:,:,j) + mask*C(l,j);
    end
end
Ic = (1-alpha)*I + alpha*cast(Lc,'like',I);
