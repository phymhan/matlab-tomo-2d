function J = flow_warp(I, flow, padding, interpolation, extrapval)

if ~exist('interpolation', 'var')
    interpolation = 'bilinear';
end
if ~exist('padding', 'var')
    padding = 'symmetric';
end
sz = [size(I,1) size(I,2)];
if sz(1)~=size(flow,1) || sz(2)~=size(flow,2)
    flow = imresize(flow, sz, 'bicubic');
end
if strcmpi(padding, 'symmetric')
    I_  = padarray(double(I), floor(sz/2), 'symmetric');
    [x_, y_] = meshgrid(...
        -floor(sz(2)/2):sz(2)+floor(sz(2)/2)-1, ...
        -floor(sz(1)/2):sz(1)+floor(sz(1)/2)-1);
elseif strcmpi(padding, 'none')
    I_ = double(I);
    [x_, y_] = meshgrid(0:sz(2)-1, 0:sz(1)-1);
else
    error('Not implemented error.')
end
[x, y] = meshgrid(0:sz(2)-1, 0:sz(1)-1);
xq = x + flow(:,:,1);
yq = y + flow(:,:,2);
J = zeros(size(I), 'like', I);
for c = 1:size(I,3)
    J(:,:,c) = interp2(x_, y_, I_(:,:,c), xq, yq, interpolation, extrapval);
end
