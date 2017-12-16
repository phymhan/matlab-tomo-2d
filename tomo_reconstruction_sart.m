function im_rec = tomo_reconstruction_sart(W, p, siz, n_it)
%  Tomographic reconstruction using SART method.
%
%   08-Aug-2013 20:44:52
%   hanligong@gmail.com

if nargin < 4
    n_it = 100;
end
%siz = size(im);
%[W, p, ~, ~] = buildWeightMatrix(im,angles);
f = solver_sart(W,p,n_it);
im_rec = reshape(f,siz);
end
