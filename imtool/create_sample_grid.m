function [saGrid, sx, sy] = create_sample_grid(imageSize, opt)
% opt: 'trombone', 'zun', 'gaussian', 'quadratic'
if ~exist('opt', 'var') || isempty(opt)
    opt = 'trombone';
end

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
saGrid = cat(3, sx, sy);