% im = imtest('head',128);
% a = 0:6:179;
% sz = size(im);
% [projmat,~] = tomoproj2d(im,a);
% [W,ix] = buildWeightMatrixArea(sz,a,0);
% [W,p] = combWeightProj(W,ix,projmat);
% im_rec = tomo_recon_lsqr(W,p,sz,1e-3,3500);
% im_rec = uint8(imscale(im_rec));
% imshow(im_rec);

fuck
sz = [192 192];
a0 = 0;
da = 2;
a1 = 179.99;
a = a0:da:a1;
[W,ix] = buildWeightMatrixArea(sz,a,0);
filename = [num2str(sz(1)) 'x' num2str(sz(2)) '_' ...
    num2str(a0) '..' num2str(da) '..' num2str(a1) '.mat'];
save(filename,'W','ix')
figure('color',[0 1 0])
