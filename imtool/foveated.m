function im1 = foveated(im0, filterBank, opt)

if ~exist('opt', 'var') || isempty(opt)
    opt = 'imfilter';
end

if ~exist('filterBank', 'var') || isempty(filterBank)
    imageSize = size(im0);
    filterSize = floor(max(imageSize)/6.4/2)*2+1;
    filterNum = 20;
    filterBank = create_filter_bank(filterSize, filterNum);
else
    imageSize = size(im0);
    filterSize = size(filterBank, 1);
    filterNum = size(filterBank, 4);
end

im0 = double(im0);
filterBank = double(filterBank);

x = linspace(-1, 1, imageSize(2));
y = linspace(-1, 1, imageSize(1));
[X, Y] = meshgrid(x, y);
I = sqrt(X.^2+Y.^2);
I = round(I/sqrt(2)*(filterNum-1)+1);

c = size(im0, 3);
im1 = zeros([imageSize(1) imageSize(2) c]);

for j = 1:size(im0, 3)
    im0_ = im0(:,:,j);
    if strcmpi(opt, 'imfilter')
        % % use imfilter
        IM1 = zeros([imageSize(1) imageSize(2) filterNum], 'double');
        for i = 1:filterNum
            IM1(:,:,i) = imfilter(im0_, filterBank(:,:,1,i), 'symmetric');
        end
    elseif strcmpi(opt, 'vl_nnconv')
        % % use vl_nnconv
        % pad = floor((fsize-1)/2);
        % IM1 = vl_nnconv(im0, filterBank, zeros(fnum, 1), 'pad', pad);
        im0_ = single(im0_);
        filterBank = single(filterBank);
        pad = floor((filterSize-1)/2);
        IM1 = vl_nnconv(im0_, filterBank, zeros(filterNum, 1, 'single'), 'pad', pad);
        IM1 = double(IM1);
    else
        error('Not implemented error.');
    end
    % % pack
    for i = 1:filterNum
        ix = I == i;
        IM1(:,:,i) = IM1(:,:,i).*ix;
    end
    im1(:,:,j) = sum(IM1, 3);
end

%%
function fb = create_filter_bank(fsize, fnum)
sigma = linspace(0.05, (fsize-1)/4, fnum);
fb = zeros([fsize fsize 1 fnum]);
for i = 1:fnum
    fb(:,:,1,i) = fspecial('Gaussian', fsize, sigma(i));
end
