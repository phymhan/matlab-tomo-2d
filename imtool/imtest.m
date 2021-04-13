function im = imtest(im_type,m,n)
%IMTEST   Generate test image\
%  Testing Images:
%
%
%   08-Aug-2013 19:18:05
%   hanligong@gmail.com

switch nargin
    case 0
        im_type = 'phan';
        m = 256;
        n = 256;
    case 1
        if ~ischar(im_type)
            switch im_type
                case 0
                    im_type = 'c';
                case 1
                    im_type = 'c1';
                case 2
                    im_type = 'c2';
                case 3
                    im_type = 'c3';
                case 4
                    im_type = 'ring';
                case 5
                    im_type = 's';
                case 6
                    im_type = 's2';
                case 7
                    im_type = 's4';
                case 8
                    im_type = 'sq';
                case 9
                    im_type = 'phantom';
                case 10
                    im_type = 'lena';
                case 12
                    im_type = 'cat';
                case 13
                    im_type = 'head';
                case 14
                    im_type = 'gauss';
                case 15
                    im_type = 'noise';
                otherwise
                    im_type = 'non';
            end
        end
        m = 256;
        n = 256;
    case 2
        if length(m) ~= 2
            n = m;
        else
            n = m(2);
            m = m(1);
        end
end
M = 256;
N = 256;
[x,y] = meshgrid(1:M,1:N);
x = x';
y = y';
if strcmpi(im_type,'c')
    im = zeros(M,N);
    im((x-N/2).^2+(y-M/2).^2 <= 32^2) = 255;
elseif strcmpi(im_type,'c1')
    im = zeros(M,N);
    im((x-N/4).^2+(y-M/4).^2 <= 32^2) = 255;
elseif strcmpi(im_type,'c2')
    im = zeros(M,N);
    im((x-  N/4).^2+(y-M/2).^2 <= 48^2) = 255;
    im((x-3*N/4).^2+(y-M/2).^2 <= 16^2) = 255;
elseif strcmpi(im_type,'c3')
    im = zeros(M,N);
    im((x-  N/2).^2+(y-3*M/4).^2 <= 16^2) = 255;
    im((x-3*N/4).^2+(y-  M/4).^2 <= 16^2) = 255;
    im((x-  N/4).^2+(y-  M/4).^2 <= 16^2) = 255;
elseif strcmpi(im_type,'ring')
    im = zeros(M,N);
    im((x-N/2).^2+(y-M/2).^2<=102^2 & (x-N/2).^2+(y-M/2).^2>=90^2) = 255;
elseif strcmpi(im_type,'s')
    im = zeros(M,N);
    l = 32;
    for k1 = 0:M/l-1
        im(k1*l+1:k1*l+l/2,:) = 255;
    end
elseif strcmpi(im_type,'s2')
    im = zeros(M,N);
    l = 32;
    for k1 = 0:M/l-1
        im(k1*l+1:k1*l+l/2,1:N/2) = 255;
    end
    for k2 = 0:(N/2)/l-1
        im(:,N/2+(k2*l+1:k2*l+l/2)) = 255;
    end
elseif strcmpi(im_type,'s4')
    im = zeros(M,N);
    l = 32;
    for k1 = 0:(M/2)/l-1
        im((k1*l+1:k1*l+l/2),1:N/2) = 255;
    end
    for k2 = 0:(N/2)/l-1
        im(M/2+1:M,(k2*l+1:k2*l+l/2)) = 255;
    end
    for k1 = 0:(M/2)/l-1
        im(M/2+(k1*l+1:k1*l+l/2),N/2+1:N) = 255;
    end
    for k2 = 0:(N/2)/l-1
        im(1:M/2,N/2+(k2*l+1:k2*l+l/2)) = 255;
    end
elseif strcmpi(im_type,'sq')
    im = zeros(M,N);
    l = 64;
    im((x-(M/2-l/2)).*(x-(M/2+l/2))<=0 & ...
        (y-(N/2-sqrt(3)*l/2)).*(y-(N/2+sqrt(3)*l/2))<=0) = 255;
elseif strcmpi(im_type,'phantom')
    im = phantom(256);
    im = im*255;
elseif strcmpi(im_type,'lena')
    S = load('lena.mat');
    im = S.im;
    im = imresize(im,[M,N]);
elseif strcmpi(im_type,'cat')
    S = load('cat.mat');
    im = S.im;
    im = imresize(im,[M,N]);
elseif strcmpi(im_type,'head')
    S = load('head.mat');
    im = S.im;
    im = imresize(im,[M,N]);
elseif strcmpi(im_type,'gauss')
    im = 255*exp(-((x-M/2).^2+(y-N/2).^2)/64^2);
elseif strcmpi(im_type,'noise')
    im = 255*rand(M,N);
else
    im = zeros(M,N);
end
if m~=M || n~=N
    im = imresize(im,[m,n]);
end
im = uint8(im);
if nargout == 0
    imshow(im);
end
end