function [proj_mat, D] = tomo_projection_2d(im,angles)
%TOMO_PROJECTION_2D   [proj_mat,D] = tomo_projection_2d(im,angles)
%   unit of angles: Degree;
%   projection direction: When angle == 0, X-ray passes through vertical 
%   (up-down) direction.
%
%   02-Aug-2013 14:07:06
%   hanligong@gmail.com

%Pad image
% [im_pad,D] = impad(im);
[im_pad,D] = impad(im,'diagonal',1); %set im_margin to 1
%Calculate projection
num_proj = length(angles);
proj_mat = zeros(num_proj,D);
for k = 1:num_proj
    im_rot = imrotate(im_pad,-angles(k),'bilinear','crop');
    proj_mat(k,:) = sum(im_rot,1);
end
end
