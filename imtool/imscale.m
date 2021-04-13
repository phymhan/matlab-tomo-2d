function [im, scl] = imscale(im, scl, otyp)
% function im = imscale(im, scl)
%   im, image
%   scl, can by a 1-by-4 vector that defines the low, high, and range value for
%    scaling (scl = [im_lo im_hi sc_lo sc_hi]). scl can also be 'auto' for auto 
%    scaling, or a function handle.
% 11/14/2015

currtyp = class(im);
if ~exist('otyp', 'var')
    otyp = currtyp; % or 'uint8' for images
elseif strcmp(otyp, 'same') || isempty(otyp)
    otyp = currtyp;
end
if ~exist('scl', 'var')
    scl = 'auto';
elseif isempty(scl)
    scl = 'same';
end
% scl: [im_lo im_hi sc_lo sc_hi]
if isa(scl, 'function_handle')
    fun = scl;
elseif strcmp(scl, 'same')
    % fun = @(M) M;
    im_casting(otyp);
    return
else
    if strcmp(scl, 'auto')
        im_lo = min(im(:));
        im_hi = max(im(:));
        sc_lo = 0;
        sc_hi = 255;
    elseif length(scl) == 2
        im_lo = min(im(:));
        im_hi = max(im(:));
        sc_lo = scl(1);
        sc_hi = scl(2);
    else
        % scl: [im_lo im_hi sc_lo sc_hi]
        im_lo = scl(1);
        im_hi = scl(2);
        sc_lo = scl(3);
        sc_hi = scl(4);
    end
    a = im_lo-sc_lo;
    if im_hi == im_lo
        l = 1;
    else
        l = (sc_hi-sc_lo)/(im_hi-im_lo);
    end
    fun = @(M) (M-a).*l;
    scl = [im_lo im_hi sc_lo sc_hi];
end
% Scaling
im = fun(im);
% Casting
im_cast(otyp);

    function im_cast(outType)
        if ~isa(currtyp, outType)
            im = cast(im, outType);
        end
    end

end

% %% old version here: 
% % Scale matrix im to a given range (default: 0~255)
% 
% if nargin < 3
%     cut = [-inf inf];
%     if nargin < 2
%         scl = 255;
%     end
% end
% im((im < cut(1)) | (im > cut(2))) = nan;
% m_max = max(im(:));
% m_min = min(im(:));
% im = (im-m_min)/(m_max-m_min)*scl;
% end
