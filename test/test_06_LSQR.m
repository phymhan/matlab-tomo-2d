%TEST 06: LSQR Method
clear
clc
fprintf('TEST 06: LSQR Method\r');
im = imtest('phan',128);
figure('name','Original Image: Phantom')
imshow(im)
siz = size(im);
load weightmat.mat
%10 Angles
fprintf('10 angles:\r')
fprintf('iteration 1000\r')
im_rec1 = tomo_recon_lsqr(W1, p1, siz, 1e-6, 1000);
im_rec1 = uint8(imscale(im_rec1));
figure('name','10 angles, LSQR(iteration 1000)')
imshow(im_rec1);
%30 Angles
fprintf('30 angles:\r')
fprintf('iteration 1000\r')
im_rec2 = tomo_recon_lsqr(W2, p2, siz, 1e-6, 1000);
im_rec2 = uint8(imscale(im_rec2));
figure('name','30 angles, LSQR(iteration 1000)')
imshow(im_rec2);
%60 Angles
fprintf('60 angles:\r')
fprintf('iteration 1000\r')
im_rec3 = tomo_recon_lsqr(W3, p3, siz, 1e-6, 1000);
im_rec3 = uint8(imscale(im_rec3));
figure('name','60 angles, LSQR(iteration 1000)')
imshow(im_rec3);
%180 Angles
fprintf('180 angles:\r')
fprintf('iteration 1000\r')
im_rec4 = tomo_recon_lsqr(W4, p4, siz, 1e-6, 1000);
im_rec4 = uint8(imscale(im_rec4));
figure('name','180 angles, LSQR(iteration 1000)')
imshow(im_rec4);
