%TEST 08: Rand v.s. All
clear
clc
fprintf('TEST 08: Rand v.s. All\r');
im = imtest('phan',175);
angles = 0:3:180;
figure('name','Original Image: Phantom')
imshow(im)
siz = size(im);
eq_num_ratio = 0.5;
load weightmat.mat
[W,p] = buildWeightMatrixArea(im,angles);
%180 Angles, W4 p4
n_eq_all = size(W,1);
n_eq_rand = ceil(n_eq_all*eq_num_ratio);
ix = randperm(n_eq_all,n_eq_rand);
fprintf('60 angles, iteration 3500\r')
im_rec = tomo_recon_lsqr(W(1:n_eq_all,:), p(1:n_eq_all,:), siz, 1e-3, 3500);
im_rec = uint8(imscale(im_rec));
figure('name','60 angles, use all equations')
imshow(im_rec);
im_rec = tomo_recon_lsqr(W(ix,:), p(ix), siz, 1e-3, 3500);
im_rec = uint8(imscale(im_rec));
figure('name','60 angles, use random equations')
imshow(im_rec);
