function im1 = nonuniform_sampling(im0, sampleGrid, opt)

if ~exist('opt', 'var') || isempty(opt)
    opt = 'trombone';
end
if ~exist('sampleGrid', 'var') || isempty(sampleGrid)
    imageSize = size(im0);
    [sx, sy] = create_sample_grid(imageSize, opt);
else
    imageSize = size(im0);
    [sx, sy] = deal(sampleGrid(:,:,1), sampleGrid(:,:,2));
end
im0 = double(im0);
x = linspace(-1, 1, imageSize(2));
y = linspace(-1, 1, imageSize(1));
% im1 = interp2(x, y, im0, sx, sy, 'linear', 0);
im1 = zeros(size(im0), 'like', im0);
for c = 1:size(im0, 3)
    im1_ = interp2(x, y, im0(:,:,c), sx, sy);
    im1_(isnan(im1_)) = mean(im1_(~isnan(im1_)));
    im1(:,:,c) = im1_;
end

%%
function [sx, sy] = create_sample_grid(imageSize, opt)
% % opt: 'trombone', 'zun', 'gaussian', 'quadratic'

% sample points range: [-1 1]
x = linspace(-1, 1, imageSize(2));
y = linspace(-1, 1, imageSize(1));
[X, Y] = meshgrid(x, y);

n = 2;
switch lower(opt)
    case 'zun'
        r = max(cat(3, 1-(X-1).^n, 1-(X+1).^n, 1-(Y-1).^n, 1-(Y+1).^n), [], 3);
    case 'trombone'
        r = sqrt(X.^2+Y.^2).^(1/n)/sqrt(2)^(1/n);
    case 'gaussian'
        r = 1-exp(-(X.^2+Y.^2)/0.4);
    case 'quadratic'
        r = (X.^n+Y.^n)/2;
end
[TH, R] = cart2pol(X, Y);
R = R.*r;
[sx, sy] = pol2cart(TH, R);
