%TEST 11: Continuity of f(x,y)
clear
clc
fprintf('TEST 11: Continuity of f(x,y)\r');
angles = 0:1:179.999;
%Phantom
im = imtest('noise',256);
im = double(im);
[M,N] = size(im);
D = ceil(sqrt(M^2+N^2));
M_pad = ceil((D-M)/2)+1;
N_pad = ceil((D-N)/2)+1;
stp = 2;
% %
MM = ceil(M/stp);
NN = ceil(N/stp);
idx = false(M,N);idx(1:stp:M,1:stp:N) = true;
figure('name','origin','color',[1 1 1])
surf(reshape(im(idx),MM,NN))
box on
camproj perspective
[projmat,~] = tomoproj2d(im,angles);
im_rec = tomo_recon_fbp(projmat,angles);
im_rec = im_rec(M_pad:D-M_pad,N_pad:D-N_pad);
im_rec = imscale(im_rec);
M = M-2;
N = N-2;
MM = ceil(M/stp);
NN = ceil(N/stp);
idx = false(M,N);idx(1:stp:M,1:stp:N) = true;
figure('name','recon','color',[1 1 1])
surf(reshape(im_rec(idx),MM,NN))
box on
camproj perspective
im = im(2:M+1,2:N+1);
err = im-im_rec;
figure('name','err','color',[1 1 1])
surf(reshape(err(idx),MM,NN))
box on
camproj perspective
im_std = std(err(:));
fprintf('standard deviation: %f\r',im_std);
