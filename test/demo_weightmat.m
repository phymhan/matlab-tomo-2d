clear;
clc;
im = imtest('phantom',[128,128]);
sz = size(im);
ang = 0:12:179.99;
%% Method 1
fprintf('The first way to call ''build_weight_matrix''\r')
[proj_mat,D] = tomo_projection_2d(im,ang);
[W, ix] = build_weight_matrix(im, ang, false, 'area');
p = proj_mat';
p = p(:);
W = sparse(W);
W(~ix,:) = [];
p(~ix) = [];
im_rec = tomo_reconstruction_lsqr(W, p, sz, 0.5*1e-3, 200);
im_rec = uint8(imscale(im_rec));
figure('name','The first method')
imshow(im_rec);
%% Method 2
fprintf('The second way to call ''build_weight_matrix''\r')
[W, p, ~, ~] = build_weight_matrix(im, ang, true, 'area');
im_rec = tomo_reconstruction_lsqr(W, p, sz, 0.5*1e-3, 200);
im_rec = uint8(imscale(im_rec));
figure('name','The second method')
imshow(im_rec);
