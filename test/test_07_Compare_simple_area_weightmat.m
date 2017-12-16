%TEST 07: Simple v.s. Area (two methods of building weight matrix)
clear
clc
pic = 'phan';
fprintf('TEST 07: Simple v.s. Area (two methods of building weight matrix)\r');
im = imtest(pic,128);
figure('name',['Original Image: ' pic])
imshow(im)
siz = size(im);
% angles = 0:2.5:179;
% [W1,p1] = buildWeightMatrixSimple(im,angles);
% [W2,p2] = buildWeightMatrixArea(im,angles);
% %Simple
% fprintf('Reconstruct from simple-W:\r')
% im_rec1 = tomo_recon_lsqr(W1, p1, siz, 1e-3, 3500);
% %im_rec1 = tomo_recon_myart(W1, p1, siz, 500);
% im_rec1 = uint8(imscale(im_rec1));
% figure('name','Simple')
% imshow(im_rec1);
% %Area
% fprintf('Reconstruct from area-W:\r')
% im_rec2 = tomo_recon_lsqr(W2, p2, siz, 1e-3, 3500);
% %im_rec2 = tomo_recon_myart(W2, p2, siz, 500);
% im_rec2 = uint8(imscale(im_rec2));
% figure('name','Area')
% imshow(im_rec2);

% angles = 0:6:179;
% [W,p] = buildWeightMatrixArea(im,angles);
% im_rec = tomo_recon_lsqr(W, p, siz, 1e-3, 3500);
% im_rec = uint8(imscale(im_rec));
% figure('name','180')
% imshow(im_rec);
% angles = 0:4.5:134;
% [W,p] = buildWeightMatrixArea(im,angles);
% im_rec = tomo_recon_lsqr(W, p, siz, 1e-3, 3500);
% im_rec = uint8(imscale(im_rec));
% figure('name','135')
% imshow(im_rec);
% angles = 0:3:79;
% [W,p] = buildWeightMatrixArea(im,angles);
% im_rec = tomo_recon_lsqr(W, p, siz, 1e-3, 3500);
% im_rec = uint8(imscale(im_rec));
% figure('name','90')
% imshow(im_rec);
% angles = 0:1.5:44;
% [W,p] = buildWeightMatrixArea(im,angles);
% im_rec = tomo_recon_lsqr(W, p, siz, 1e-3, 3500);
% im_rec = uint8(imscale(im_rec));
% figure('name','45')
% imshow(im_rec);

% [M,N] = size(im);
% D = ceil(sqrt(M^2+N^2));
% M_pad = ceil((D-M)/2)+1;
% N_pad = ceil((D-N)/2)+1;
% angles = 0:2:179;
% [W,p] = buildWeightMatrixArea(im,angles);
% im_rec = tomo_recon_lsqr(W, p, siz, 1e-3, 3500);
% im_rec = uint8(imscale(im_rec));
% figure('name','lsqr180')
% imshow(im_rec);
% [pm,~] = tomoproj2d(im,angles);
% im_rec = tomo_recon_fbp(pm,angles);
% im_rec = im_rec(M_pad:D-M_pad,N_pad:D-N_pad);
% im_rec = uint8(imscale(im_rec));
% figure('name','fbp180')
% imshow(im_rec);
% angles = 0:1:89;
% [W,p] = buildWeightMatrixArea(im,angles);
% im_rec = tomo_recon_lsqr(W, p, siz, 1e-3, 3500);
% im_rec = uint8(imscale(im_rec));
% figure('name','lsqr90')
% imshow(im_rec);
% [pm,~] = tomoproj2d(im,angles);
% im_rec = tomo_recon_fbp(pm,angles);
% im_rec = im_rec(M_pad:D-M_pad,N_pad:D-N_pad);
% im_rec = uint8(imscale(im_rec));
% figure('name','fbp90')
% imshow(im_rec);

% angles = 0:3:179;
% [W,p] = buildWeightMatrixArea(im,angles);
% im_rec = tomo_recon_lsqr(W, p, siz, 1e-3, 2000);
% im_rec = uint8(imscale(im_rec));
% figure('name','lsqr')
% imshow(im_rec);
% im_rec = tomo_recon_myart(W, p, siz, 2000);
% im_rec = uint8(imscale(im_rec));
% figure('name','myart')
% imshow(im_rec);

% angles = 0:6:89;
% [W,p] = buildWeightMatrixArea(im,angles);
% im_rec = tomo_recon_lsqr(W, p, siz, 1e-3, 3500);
% im_rec = uint8(imscale(im_rec));
% figure('name','0:6:90')
% imshow(im_rec);
% angles = 0:4:89;
% [W,p] = buildWeightMatrixArea(im,angles);
% im_rec = tomo_recon_lsqr(W, p, siz, 1e-3, 3500);
% im_rec = uint8(imscale(im_rec));
% figure('name','0:4:90')
% imshow(im_rec);
% angles = 0:2:89;
% [W,p] = buildWeightMatrixArea(im,angles);
% im_rec = tomo_recon_lsqr(W, p, siz, 1e-3, 3500);
% im_rec = uint8(imscale(im_rec));
% figure('name','0:2:90')
% imshow(im_rec);
% angles = 0:1:89;
% [W,p] = buildWeightMatrixArea(im,angles);
% im_rec = tomo_recon_lsqr(W, p, siz, 1e-3, 3500);
% im_rec = uint8(imscale(im_rec));
% figure('name','0:1:90')
% imshow(im_rec);
% angles = 0:0.5:89;
% [W,p] = buildWeightMatrixArea(im,angles);
% im_rec = tomo_recon_lsqr(W, p, siz, 1e-3, 3500);
% im_rec = uint8(imscale(im_rec));
% figure('name','0:0.5:90')
% imshow(im_rec);
% angles = 0:0.25:89;
% [W,p] = buildWeightMatrixArea(im,angles);
% im_rec = tomo_recon_lsqr(W, p, siz, 1e-3, 3500);
% im_rec = uint8(imscale(im_rec));
% figure('name','0:0.25:90')
% imshow(im_rec);
% angles = 0:0.1:89;
% [W,p] = buildWeightMatrixArea(im,angles);
% im_rec = tomo_recon_lsqr(W, p, siz, 1e-3, 3500);
% im_rec = uint8(imscale(im_rec));
% figure('name','0:0.1:90')
% imshow(im_rec);
% save('phan_0..0.1..90.mat','W','p','angles')
