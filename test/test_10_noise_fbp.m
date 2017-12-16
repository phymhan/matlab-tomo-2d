%TEST 10: Noise FBP
clear
clc
fprintf('TEST 10: Noise FBP\r');
angles = 0:1:179.999;
%Phantom
im = imtest('phan',256);
im = im+uint8(abs(25.5*5*rand(256,256)+double(im)));
[M,N] = size(im);
D = ceil(sqrt(M^2+N^2));
M_pad = ceil((D-M)/2)+1;
N_pad = ceil((D-N)/2)+1;
figure('name','Original Image: Phantom')
imshow(im)
[projmat,~] = tomoproj2d(im,angles);
im_rec = tomo_recon_anisofbp(projmat,angles);
im_rec = im_rec(M_pad:D-M_pad,N_pad:D-N_pad);
im_rec = uint8(imscale(im_rec));
figure
imshow(im_rec)
im_rec = tomo_recon_fbp(projmat,angles);
im_rec = im_rec(M_pad:D-M_pad,N_pad:D-N_pad);
im_rec = uint8(imscale(im_rec));
figure
imshow(im_rec)