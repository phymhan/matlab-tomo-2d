clear all;
close all;

I=imread('flame1.jpg');
figure;
imshow(I);

I1=rgb2gray(I);
im = I1;

%Creating Sinogram.
fprintf(1,'Demo 01: Compute projections\r');
% im=imtest('I1');
ang=0:1:179.99;
[projmat,~]=tomo_projection_2d(im,ang);
improj=uint8(imscale(projmat));
figure('name','Projections');
imshow(improj);

%Build weight matrix.
fprintf(1,'Demo 02: Build weight matrix\r');
sz=[128,128];
ang=0:2:179.99;
[W,~]=build_weight_matrix(sz,ang,0,'area');
figure('name','Weight Matrix');
spy(W);

%ART reconstructions.
fprintf(1,'Demo 03:ART reconstruction\r');
clear;
clc;
im=imtest('I1',[256,256]);
sz=size(im);
ang=0:3:179.99;
[W,p,~,~]=build_weight_matrix(im,ang,1,'area');
fprintf('ART solver\r')
im_rec=tomo_reconstruction_art(W,p,sz);
im_rec=uint8(imscale(im_rec));
figure('name','ART')
imshow(im_rec);

fprintf(1,'Demo 04: ART reconstruction\r');
clear;
clc;
im=imtest('I1',[256,256]);
sz=size(im);
ang=0:12:179.99;
[W,p,~,~]=build_weight_matrix(im,ang,1,'area');
fprintf('ART solver\r')
im_rec=tomo_reconstruction_art(W,p,sz);
im_rec=uint8(imscale(im_rec));
figure('name','ART')
imshow(im_rec);