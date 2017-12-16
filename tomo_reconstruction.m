function im_rec = tomo_reconstruction(im,angles,method,n_it)
% TOMO_RECON   FBP, BP, ART
%
%   09-Aug-2013 10:19:03
%   hanligong@gmail.com
%

switch nargin
    case 3
        n_it = 100;
    case 2
        method = 'fbp';
end
if strcmpi(method,'fbp')
    [projmat,~] = tomo_projection_2d(im,angles);
    im_rec = tomo_reconstruction_fbp(projmat,angles);
elseif strcmpi(method,'bp')
    [projmat,~] = tomoproj2d(im,angles);
    im_rec = tomo_reconstruction_bp(projmat,angles);
elseif strcmpi(method,'art')
    im_rec = tomo_reconstruction_art_im(im,angles,n_it);
elseif strcmpi(method,'sart')
    im_rec = tomo_reconstruction_sart_im(im,angles,n_it);
elseif strcmpi(method,'lsqr')
    im_rec = tomo_reconstruction_lsqr_im(im,angles,1e-6,n_it);
else
    fprintf('Comming soon...\r')
    im_rec = [];
    return
end
