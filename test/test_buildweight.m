function test_buildweight
% This function is wrote for the purpose of debugging.

%Pad image
im = ones(4,5);
%disp(im)
angles = [23 67];
[im_pad,D] = impad(im,'diag',1); %avoid <=1 or >=D indices
%disp(im_pad)
sz = size(im);
m = sz(1); %rows
n = sz(2); %cols
n_p = length(angles);
N = m*n; %number of unknown variables
M = D*n_p; %number of equations
%Weighting Factor Matrix
W = zeros(M,N);
%proj_info = zeros(M,2); %[theta No.]
projmat = zeros(n_p,D); %projection matrix
p = zeros(M,1); %projection vector
ix = false(M,1);
%=========================%
rc = repmat([m;n]/2+0.5,1,N);
rcy = repmat(n/2+0.5,1,N);
%xypad = repmat(floor([(D-m)/2;(D-n)/2]),1,N);
ypad = repmat(floor((D-n)/2),1,N);
%=========================%
%cnt = 0;
fprintf('Building Weight Matrix...\r')
[x,y] = ind2sub(sz,1:N);
%xy = [x;y]; %x,y coordinate
%x,y coordinate with respect to rotation center
xy_c = [x;y]-rc;
%Parameter for accounting area ratio
d = 1; %side length of square
A0 = d^2;
H = @heaviside;
%% lines
% %
figure
axes('xlim',[0 D+1],'ylim',[0 D+1],'dataaspectratio',[1 1 1],...
    'xtick',1:D,'ytick',1:D,'box','on')
%padded image
for k = 0:D
    line([k+0.5 k+0.5],[0.5 D+0.5],'linestyle',':','color',[0 0 0]);
    line([0.5 D+0.5],[k+0.5 k+0.5],'linestyle',':','color',[0 0 0]);
end
%original image
xp = floor((D-m)/2);
yp = floor((D-n)/2);
lxh = zeros(1,m+1);
lyh = zeros(1,n+1);
%gx = xp+0.5:1:xp+m+0.5;
%gy = yp+0.5:1:yp+n+0.5;
gx = 0.5:1:m+0.5;
gy = 0.5:1:n+0.5;
for k = 1:m+1
    %lxh(k) = line([gx(k) gx(k)],[yp+0.5 yp+n+0.5],'color',[0 0 1]);
    lxh(k) = line([gx(k) gx(k)]+xp,[yp+0.5 yp+n+0.5],'color',[0 0 1]);
end
for k = 1:n+1
    %lyh(k) = line([xp+0.5 xp+m+0.5],[gy(k) gy(k)],'color',[0 0 1]);
    lyh(k) = line([xp+0.5 xp+m+0.5],[gy(k) gy(k)]+yp,'color',[0 0 1]);
end
txth = zeros(D,D);
for ii = 1:D
    for jj = 1:D
        txth(ii,jj) = text(ii,jj,num2str(im_pad(ii,jj)),...
            'HorizontalAlignment','center');
    end
end
pth = line(NaN,NaN,'marker','o','markeredgecolor',[1 0 0]);
lphc = line([0 D+1],[NaN NaN],'linestyle',':','color',[1 0 0]);
lphu = line([0 D+1],[NaN NaN],'color',[1 0 0]);
lphd = line([0 D+1],[NaN NaN],'color',[1 0 0]);
ly0h = line([0 D+1],[NaN NaN],'color',[0 1 0]);
% %
%% Compute
for kp = 1:n_p
    im_rot = imrotate(im_pad,-angles(kp),'bilinear','crop');
    pvec = sum(im_rot,1);
    projmat(kp,:) = pvec;
    t = -angles(kp)/180*pi;
    RR = [cos(t) -sin(t);sin(t) cos(t)];
    R = [sin(t) cos(t)];
    % %
    for k = 1:m+1
        %lxy1 = RR*([gx(k);yp+0.5]-[D/2+0.5;D/2+0.5])+[D/2+0.5;D/2+0.5];
        %lxy2 = RR*([gx(k);yp+n+0.5]-[D/2+0.5;D/2+0.5])+[D/2+0.5;D/2+0.5];
        lxy1 = RR*([gx(k);0.5]-[m/2+0.5;n/2+0.5])+[m/2+0.5;n/2+0.5]+[xp;yp];
        lxy2 = RR*([gx(k);n+0.5]-[m/2+0.5;n/2+0.5])+[m/2+0.5;n/2+0.5]+[xp;yp];
        set(lxh(k),'xdata',[lxy1(1) lxy2(1)],'ydata',[lxy1(2) lxy2(2)]);
    end
    for k = 1:n+1
        %lxy1 = RR*([xp+0.5;gy(k)]-[D/2+0.5;D/2+0.5])+[D/2+0.5;D/2+0.5];
        %lxy2 = RR*([xp+m+0.5;gy(k)]-[D/2+0.5;D/2+0.5])+[D/2+0.5;D/2+0.5];
        lxy1 = RR*([0.5;gy(k)]-[m/2+0.5;n/2+0.5])+[m/2+0.5;n/2+0.5]+[xp;yp];
        lxy2 = RR*([m+0.5;gy(k)]-[m/2+0.5;n/2+0.5])+[m/2+0.5;n/2+0.5]+[xp;yp];
        set(lyh(k),'xdata',[lxy1(1) lxy2(1)],'ydata',[lxy1(2) lxy2(2)]);
    end
    updateImpad
    % %
    fprintf('\nAngle %d(%d Degree)\n%.1f%% completed.\r',...
        kp,angles(kp),100*kp/n_p)
    %=========================%
    xy_rot = R*xy_c+rcy+ypad; %!!!
    %=========================%
    %idx = round(xy_rot(2,:));
    idx = round(xy_rot);
    %corresponding indice in W and p matrix
    ixM = D*(kp-1)+idx;
    %Calculate Area
    x1 = d*min(abs([sin(t) cos(t)]))+eps;
    x2 = d*max(abs([sin(t) cos(t)]))-eps;
    l = x1+x2;
    h = d^2/x2;
    %A0 = d^2;
    A1 = x1*h/2;
    G1 = @(x) H(x)-H(x-x1);
    G2 = @(x) H(x-x1)-H(x-x2);
    G3 = @(x) H(x-x2)-H(x-l);
    f_A = @(x) (x/x1).^2*A1.*G1(x)+...
        (A1+h*(x-x1)).*G2(x)+...
        (A0-((l-x)/x1).^2.*A1).*G3(x)+...
        A0*H(x-l);
    %x0 = xy_rot(2,:)-l/2;
    x0 = xy_rot-l/2;
    S1 = f_A((idx-0.5)-x0)/A0;
    S2 = f_A((idx+0.5)-x0)/A0;
    for kn = 1:N
        W(ixM(kn)-1,kn) = S1(kn);
        W(ixM(kn)  ,kn) = S2(kn)-S1(kn);
        W(ixM(kn)+1,kn) = 1-S2(kn);
        % %
        [xx,yy] = ind2sub(sz,kn);
        xxyy = RR*([xx;yy]-[m/2+0.5;n/2+0.5])+[m/2+0.5;n/2+0.5]+[xp;yp];
        set(pth,'xdata',xxyy(1),'ydata',xxyy(2));
        disp(xxyy(2))
        disp(idx(kn))
        set(lphc,'ydata',[idx(kn) idx(kn)]);
        set(lphu,'ydata',[idx(kn) idx(kn)]+0.5);
        set(lphd,'ydata',[idx(kn) idx(kn)]-0.5);
        set(ly0h,'ydata',[x0(kn) x0(kn)]);
        pause(0.5) %#
        fprintf('\narea1: %.2f\n',S1(kn));
        fprintf('\narea2: %.2f\n',S2(kn)-S1(kn));
        fprintf('\narea3: %.2f\n',1-S2(kn));
        % %
        %if ~ix(ixM(kn))
        %    p(ixM(kn)) = pvec(idx(kn));
        %    ix(ixM(kn)) = true;
        %end
        %=========================%
        %p(ixM(kn)) = projmat(kp,idx(kn));
        p(ixM(kn)) = pvec(idx(kn));
        % %
        disp(W((kp-1)*D+idx(kn),kn));
        disp(pvec(idx(kn)));
        % %
        %=========================%
        ix(ixM(kn)) = true;
    end
end
%% Delete all-zero rows in W
fprintf('\nDelete all-zero rows in W...\r')
%ix = sum(W,2)==0;
try
    W(~ix,:) = [];
catch expr
    fprintf([expr.message '\nConverting to sparse...\r'])
    W = sparse(W);
    W(~ix,:) = [];
end
p(~ix) = [];
n_eq = size(W,1);
fprintf('\nFinish building A. \n%d equations in total.\r',n_eq);
figure('name','Sparsity pattern of matrix W','numbertitle','off');
spy(W)
e = W*im(:)-p;
disp(e)

    function updateImpad
        for fii = 1:D
            for fjj = 1:D
                set(txth(fii,fjj),'string',num2str(im_rot(fii,fjj),'%.1f'),...
                    'HorizontalAlignment','center');
            end
        end
    end
end
