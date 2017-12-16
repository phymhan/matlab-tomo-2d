%% Structure and Arrangement of this toolbox
% There are two main methods for tomographic reconstruction: one is based
% on Radon transform and its inverse transform such as filtered
% back-projection (FBP), another is based on solving linear algebra
% equations such as algebraic reconstruction technique (ART).
%
% This toolbox includes tools for creating projections and reconstructing
% the image from projections:
% tomo_projection_2d, computes projections of a given image;
% build_weight_matrix, builds weighting factor matrix used for algebraic
% methods;
% tomo_reconstruction_bp, reconstructs the image from its projections using
% BP method;
% tomo_reconstruction_fbp, reconstructs the image from its projections
% using FBP method;
% tomo_reconstruction_art, reconstructs the image from its projections
% using ART method;
% tomo_reconstruction_sart, reconstructs the image from its projections
% using SART method;
% tomo_reconstruction_lsqr, reconstructs the image from its projections,
% equations are solved by calling the build-in lsqr function;
%
%% References:
%  http://www.sv.vt.edu/xray_ct/parallel/Parallel_CT.html
%  http://www.dspguide.com/ch25/5.htm
%  http://www.owlnet.rice.edu/~elec539/Projects97/cult/node2.html
%  http://www.owlnet.rice.edu/~elec539/Projects97/cult/node3.html#SECTION00021000000000000000
%  http://www.owlnet.rice.edu/~elec539/Projects97/cult/node4.html#SECTION00022000000000000000
%  http://en.wikipedia.org/wiki/Radon_transform
%  http://en.wikipedia.org/wiki/Deconvolution
%  http://www.tnw.tudelft.nl/index.php?id=33826&binary=/doc/mvanginkel_radonandhough_tr2004.pdf
%
%% Demos

%% Show test images
fprintf(1,'Demo 00: Show test images.\r');
for k = 0:15
    figure
    imshow(imtest(k))
end

%% Test filters
load bone01;
n_it = 15;
tao = 0.1;
lambda = 30;
gn = 1;
fprintf('Demo 01: Compare Anisotropic Filter and Bilateral Filter\r');
fprintf('Parameters:\n\nAnisoDiff:\niteration %d\ntao       %f\nlambda    %f\r',n_it,tao,lambda);
fprintf('Bilateral:\ninteration %d\nmask size  %dx%d\nlambda_s   AUTO\nlambda_r   AUTO\r',n_it,2*gn+1,2*gn+1);
im0 = rgb2gray(im0);
im1 = anisodiff(im0,n_it,tao,lambda,1);
im2 = bilateral_filter(im0,n_it,gn);
figure('name','Original Image')
imshow(uint8(im0));
figure('name','Anisotropic Filter')
imshow(uint8(im1));
figure('name','Bilateral Filter')
imshow(uint8(im2));

%% Compute projections
fprintf(1,'Demo 02: Compute projections\r');
im = imtest('phantom');
ang = 0:1:179.99;
[projmat, ~] = tomo_projection_2d(im,ang);
improj = uint8(imscale(projmat));
figure('name','Projections');
imshow(improj);

%% Build weight matrix
fprintf(1,'Demo 03: Build weight matrix\r');
sz = [128,128];
ang = 0:2:179.99;
[W, ~] = build_weight_matrix(sz, ang, 0, 'area');
figure('name','Weight Matrix');
spy(W);

%% Filtered back-projection reconstruction
fprintf(1,'Demo 04: Filtered back-projection\r');
im = imtest('phantom');
ang = 0:1:179.99;
[projmat, ~] = tomo_projection_2d(im,ang);
tomo_reconstruction_fbp(projmat,ang,true);

%% Algebraic methods
fprintf(1,'Demo 05: Algebraic methods\r');
clear;
clc;
im = imtest('phantom');
sz = size(im);
ang = 0:3:179.99;
[W, p, ~, ~] = build_weight_matrix(im, ang, 1, 'area');
fprintf('LSQR solver\r')
im_rec = tomo_reconstruction_lsqr(W, p, sz, 0.5*1e-3, 200);
im_rec = uint8(imscale(im_rec));
figure('name','LSQR')
imshow(im_rec);

%% ART
fprintf(1,'Demo 06_01: ART methods\r');
clear;
clc;
im = imtest('phantom',[128,128]);
sz = size(im);
ang = 0:3:179.99;
[W, p, ~, ~] = build_weight_matrix(im, ang, 1, 'area');
fprintf('ART solver\r')
im_rec = tomo_reconstruction_art(W, p, sz);
im_rec = uint8(imscale(im_rec));
figure('name','ART')
imshow(im_rec);

fprintf(1,'Demo 06_02: ART methods\r');
clear;
clc;
im = imtest('phantom',[256,256]);
sz = size(im);
ang = 0:12:179.99;
[W, p, ~, ~] = build_weight_matrix(im, ang, 1, 'area');
fprintf('ART solver\r')
im_rec = tomo_reconstruction_art(W, p, sz);
im_rec = uint8(imscale(im_rec));
figure('name','ART')
imshow(im_rec);

%% Two ways of calling function 'build_weight_matrix'
clear;
clc;
im = imtest('phantom',[128,128]);
sz = size(im);
ang = 0:12:179.99;
% Method 1
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
% Method 2
fprintf('The second way to call ''build_weight_matrix''\r')
[W, p, ~, ~] = build_weight_matrix(im, ang, true, 'area');
im_rec = tomo_reconstruction_lsqr(W, p, sz, 0.5*1e-3, 200);
im_rec = uint8(imscale(im_rec));
figure('name','The second method')
imshow(im_rec);

%% Crop reconstructed images
im = imtest('lena',512);
angles = 0:1:179;
im = double(im);
% Compute projection matrix
[projmat,D] = tomo_projection_2d(im,angles);
% Reconstruct image
im_rec = tomo_reconstruction_fbp(projmat,angles);
% Crop reconstructed image
[m,n] = size(im);
m_pad = floor((D-m)/2);
n_pad = floor((D-n)/2);
im_crop = im_rec(m_pad+1:m_pad+m,n_pad+1:n_pad+n);
im_crop = uint8(imscale(im_crop));
figure
imshow(im_crop)
% Crop reconstructed image by calling function 'imcrop_tomo'
im_crop2 = imcrop_tomo(im_rec,[m,n]);
im_crop2 = uint8(imscale(im_crop2));
figure
imshow(im_crop2)

%% An animation illustrates how to compute weighting factor matrix
test_buildweight
