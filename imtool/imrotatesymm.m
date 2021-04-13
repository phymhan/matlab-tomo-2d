function imOut = imrotatesymm(A, angled, method)
% function imOut = IMROTATESYMM(A, angled, method)
% July 18, 2016

if ~exist('method', 'var') || isempty(method)
    method = 'bilinear';
end
[height, width, ~] = size(A);
imLr = fliplr(A);
imUd = flipud(A);
imPi = fliplr(imUd);
imSymm = cat(1, cat(2, imPi, imUd, imPi), cat(2, imLr, A, imLr), ...
    cat(2, imPi, imUd, imPi));
imRot = imrotate(imSymm, angled, method, 'crop');
trans = [0 0];
imOut = imRot(trans(1)+height+1:trans(1)+height+height, ...
    trans(2)+width+1:trans(2)+width+width, :);
