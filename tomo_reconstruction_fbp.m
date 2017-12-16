function im_rec = tomo_reconstruction_fbp(projmat_bp,angles,show_animation)
%  Tomographic reconstruction using Filtered Back Projection (FBP) 
%   algorithm.
%
%   05-Aug-2013 17:23:35
%   hanligong@gmail.com

if nargin < 3
    show_animation = false;
else
    figure('name','FBP')
end
[n_proj,D] = size(projmat_bp);
if length(angles) ~= n_proj
    fprintf('Size of proj_mat is incorrect.\r')
    im_rec = [];
    return
end
im_rec = zeros(D);

% N = 3;
% h = zeros(1,N);
% m = floor((N+1)/2);
% for n = 1:N
%     if n == m
%         h(n) = 1/4;
%     elseif rem(n, 2) == 0
%         h(n) = 0;
%     else
%         h(n) = -1/(pi*(n-m))^2;
%     end
% end

for k = 1:n_proj
    proj_vec = hpf(projmat_bp(k,:));
    %proj_vec = filtersinc(projmat_bp(k,:));
    %proj_vec = conv(projmat_bp(k,:), h, 'same');
    im_rot = repmat(proj_vec/D,D,1);
    im_rot = imrotate(im_rot, angles(k), 'bilinear', 'crop');
    im_rec = im_rec+im_rot/n_proj;
    %DEBUG
    if show_animation
        imshow(uint8(imscale(im_rec)))
        pause(0.01)
    end
    % %
end
