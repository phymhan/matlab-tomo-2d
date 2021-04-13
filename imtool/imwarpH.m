function imgWarp = imwarpH(img, T, outRange, magnification, interpMethod, newClass)
% function outputImage = imwarpH(inputImage, T, outRange, magnification, interpMethod, newClass)
%   [x; y; 1] = T * [u; v; 1];
%   out_range = [xmin, xmax, ymin, ymax];
% 12/04/2015

oldClass = class(img);
if ~exist('newClass', 'var')
    newClass = oldClass;
end
if ~exist('interpMethod', 'var')
    interpMethod = 'linear';
end
if ~exist('magnification', 'var')
    magnification = 1;
end
if length(outRange) == 2
    outRange = [1 outRange(2) 1 outRange(1)];
end
nc = size(img, 3);
step = 1/magnification;
[x, y] = meshgrid(outRange(1):step:outRange(2), outRange(3):step:outRange(4));
[M, N] = size(x);
coordi = T\[x(:)'; y(:)'; ones(1, M*N)];
xi = reshape(coordi(1,:)./coordi(3,:), [M N]);
yi = reshape(coordi(2,:)./coordi(3,:), [M N]);
img = double(img);
imgWarp = zeros(M, N, nc);
for k = 1:size(img, 3)
    imgWarp(:,:,k) = interp2(img(:,:,k), xi, yi, interpMethod);
end
imgWarp = cast(imgWarp, newClass);
