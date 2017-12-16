im = imtest('cat',512);
angles = 0:1:179;
im = double(im);
%% Compute projection matrix
[projmat,D] = tomo_projection_2d(im,angles);
%% Reconstruct image
im_rec = tomo_reconstruction_fbp(projmat,angles);
%% Crop reconstructed image
[m,n] = size(im);
m_pad = floor((D-m)/2);
n_pad = floor((D-n)/2);
im_crop = im_rec(m_pad+1:m_pad+m,n_pad+1:n_pad+n);
im_crop = uint8(imscale(im_crop));
figure
imshow(im_crop)
%% Crop reconstructed image by calling function 'imcrop_tomo'
im_crop2 = imcrop_tomo(im_rec,[m,n]);
im_crop2 = uint8(imscale(im_crop2));
figure
imshow(im_crop2)
