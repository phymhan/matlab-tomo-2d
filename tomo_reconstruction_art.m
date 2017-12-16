function im_rec = tomo_reconstruction_art(W, p, sz, n_it)
%Tomographic reconstruction using ART method.
%
%   08-Aug-2013 22:39:35
%   hanligong@gmail.com

if nargin < 4
    n_it = size(W,1);
end

f = solver_art(W,p,n_it);
im_rec = reshape(f,sz);

end
