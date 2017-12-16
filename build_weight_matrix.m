function [W, varargout] = build_weight_matrix(im, angles, delz, opt)
%BUILD_WEIGHT_MATRIX  Build weighting factor matrix 'W', the projection
%   vector 'p', and projection matrix 'projmat'.
%   Reconstruct the image from equation W*f == p.
%   Weighting factor matrix W is a M-by-N matrix, where M represents the 
%   number of equations and N represents the number of unknowns (i.e. 
%   N = number_of_rows * number_of_columns).
%   If the input arguments 'delz' equals false, all-zero rows in W will not
%   be deleted.
%   
%   Syntax: 
%   [W, ix] = build_weight_matrix_area(im, angles, 0)
%   [W, p, D, projmat] = build_weight_matrix_area(im, angles, 1)
%   
%   If you just want to compute W and use the same weighting factor matrix 
%   to reconstruct different images of the same size, you can call the 
%   function in this way: 
%   [W, ix] = build_weight_matrix_area(im, angles, 0), 
%   where the output argument 'ix' is a vector indicates which row in W is 
%   all-zero. Denote the diameter of the given image as 'D'. Then each 
%   projection will be a 1-by-D vector, and W will be a 
%   (D*number_of_angles) by D matrix.
%   You can also delete all-zero rows at the same time when building W, and
%   get the corresponding value of p, D and projmat by calling the function
%   in this way: 
%   [W, p, D, projmat] = build_weight_matrix_area(im, angles, 1).
%
%   25-Aug-2013 12:44:31
%   hanligong@gmail.com

if nargin < 4
    opt = 'area';
    if nargin < 3
        delz = true;
    end
end
if ~ischar(opt)
    switch opt
        case 1
            opt = 'simple';
        case 2
            opt = 'area';
        otherwise
            fprintf('%d is not a property value of ''option''.\r',opt);
            return
    end
end
%Timer
try
    t1 = toc;
catch expr
    fprintf([expr.message '\nCalling TIC to start a stopwatch timer...\r'])
    tic
    t1 = toc;
end
sz = size(im);
if length(sz) ~= 2
    fprintf(['The input argument ''im'' must be an image or ' ...
        'a vector represents its size.\r'])
    return
end
if sz(1)*sz(2) == 2
    sz(1) = im(1);
    sz(2) = im(2);
end
if delz
    %Pad image
    if strcmpi(opt,'simple')
        [im_pad,D] = impad(im,'diag',0);
    else %area
        [im_pad,D] = impad(im,'diag',1); %avoid <=1 or >=D indices
    end
else
    D = ceil(sqrt(sz*sz'))+2;
end
m = sz(1); %rows
n = sz(2); %cols
n_p = length(angles);
N = m*n; %number of unknown variables
M = D*n_p; %number of equations
%Weighting Factor Matrix
try
    W = zeros(M,N);
    method = 'matrix';
catch expr
    fprintf([expr.message '\nGenerating a sparse...\r'])
    W = sparse(M,N);
    method = 'sparse';
end
%proj_info = zeros(M,2); %[theta No.]
if delz
    %projection matrix
    projmat = zeros(n_p,D);
    %projection vector
    p = zeros(M,1);
else
    projmat = [];
end
ix = false(M,1);
%rc = [m;n]/2;
%cnt = 0;
fprintf('Building Weight Matrix...\r')
%% SLOW
% for kp = 1:n_proj
%     im_rot = imrotate(im_pad,-angles(kp),'bilinear','crop');
%     projmat(kp,:) = sum(im_rot,1);
%     t = -angles(kp)/180*pi;
%     R = [cos(t) -sin(t);sin(t) cos(t)];
%     %L = m*abs(sin(t))+n*abs(cos(t));
%     %offset = ceil((D-L)/2);
%     fprintf('\nAngle No.%d(%d Degree)\r',kp,angles(kp))
%     for kn = 1:N
%         [x,y] = ind2sub(sz,kn);
%         xy_rot = R*([x;y]-rc)+D/2+0.5;
%         idx = round(xy_rot(2));%#
%         %corresponding indice in W and p matrix
%         ixM = D*(kp-1)+idx;
%         W(ixM,kn) = 1;
%         %W(D*(kp-1)+idx,kn) = 1;
%         %W(cnt+idx-offset,kn) = 1;
%         %p(cnt+idx-offset) = projmat(kp,idx);%#
%         if ~ix(ixM)
%             p(ixM) = projmat(kp,idx);
%             ix(ixM) = true;
%         end
%         %p(D*(kp-1)+idx) = projmat(kp,idx);
%         %proj_info(cnt+idx-offset,:) = [kp,idx];
%     end
%     %cnt = cnt+(D-2*offset);
%     %fprintf('%d equations built.\r',D-2*offset)
% end

%% FAST
[x,y] = ind2sub(sz,1:N);
%xy = [x;y]; %x,y coordinate
%x,y coordinate with respect to rotation center
xy_c = [x;y]-repmat([m;n]/2+0.5,1,N);
rc = (n/2+0.5)*ones(1,N);
%xypad = repmat(floor([(D-m)/2;(D-n)/2]),1,N);
xypad = repmat(floor((D-n)/2),1,N);
if strcmpi(opt,'simple')
    if strcmp(method,'sparse')
        try
            A = zeros(D,N);
        catch expr
            fprintf(['Building A...\n' ...
                expr.message '\nGenerating a sparse...\r'])
            A = sparse(D,N);
        end
        for kp = 1:n_p
            A(:) = 0; %#
            if delz
                im_rot = imrotate(im_pad,-angles(kp),'bilinear','crop');
                pvec = sum(im_rot,1);
                projmat(kp,:) = pvec;
            end
            t = -angles(kp)/180*pi;
            %R = [cos(t) -sin(t);sin(t) cos(t)];
            R = [sin(t) cos(t)]; %!!
            fprintf('\nAngle %d(%.1f Degree)\n%.1f%% completed.\r',...
                kp,angles(kp),100*kp/n_p)
            xy_rot = R*xy_c+rc+xypad; %!!
            %idx = round(xy_rot(2,:));
            idx = round(xy_rot); %!!
            %corresponding indice in W and p matrix
            ixM = D*(kp-1)+idx;
            for kn = 1:N
                A(idx(kn),kn) = 1;
                %if ~ix(ixM(kn))
                %    p(ixM(kn)) = pvec(idx(kn));
                %    ix(ixM(kn)) = true;
                %end
                %p(ixM(kn)) = pvec(idx(kn)); %!!
                %p(ixM(kn)) = projmat(kp,idx(kn));
                ix(ixM(kn)) = true;
            end
            if delz
                for kn = 1:N
                    p(ixM(kn)) = pvec(idx(kn)); %!!
                end
            end
            W((D*(kp-1)+(1:D)),:) = A; %#
        end
    else
        for kp = 1:n_p
            if delz
                im_rot = imrotate(im_pad,-angles(kp),'bilinear','crop');
                pvec = sum(im_rot,1);
                projmat(kp,:) = pvec;
            end
            t = -angles(kp)/180*pi;
            %R = [cos(t) -sin(t);sin(t) cos(t)];
            R = [sin(t) cos(t)]; %!!
            fprintf('\nAngle %d(%.1f Degree)\n%.1f%% completed.\r',...
                kp,angles(kp),100*kp/n_p)
            xy_rot = R*xy_c+rc+xypad; %!!
            %idx = round(xy_rot(2,:));
            idx = round(xy_rot); %!!
            %corresponding indice in W and p matrix
            ixM = D*(kp-1)+idx;
            for kn = 1:N
                W(ixM(kn),kn) = 1;
                %if ~ix(ixM(kn))
                %    p(ixM(kn)) = pvec(idx(kn));
                %    ix(ixM(kn)) = true;
                %end
                %p(ixM(kn)) = pvec(idx(kn)); %!!
                %p(ixM(kn)) = projmat(kp,idx(kn));
                ix(ixM(kn)) = true;
            end
            if delz
                for kn = 1:N
                    p(ixM(kn)) = pvec(idx(kn)); %!!
                end
            end
        end
    end
elseif strcmpi(opt,'area')
    %Parameter for accounting area ratio
    d = 1; %side length of square
    A0 = d^2;
    H = @heaviside;
    if strcmp(method,'sparse')
        try
            A = zeros(D,N);
        catch expr
            fprintf(['Building A...\n' ...
                expr.message '\nGenerating a sparse...\r'])
            A = sparse(D,N);
        end
        for kp = 1:n_p
            A(:) = 0; %#
            if delz
                im_rot = imrotate(im_pad,-angles(kp),'bilinear','crop');
                pvec = sum(im_rot,1);
                projmat(kp,:) = pvec;
            end
            t = -angles(kp)/180*pi;
            %R = [cos(t) -sin(t);sin(t) cos(t)];
            R = [sin(t) cos(t)]; %!!
            fprintf('\nAngle %d(%.1f Degree)\n%.1f%% completed.\r',...
                kp,angles(kp),100*kp/n_p)
            xy_rot = R*xy_c+rc+xypad; %!!
            %idx = round(xy_rot(2,:));
            idx = round(xy_rot); %!!
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
            x0 = xy_rot-l/2; %!!
            S1 = f_A((idx-0.5)-x0)/A0;
            S2 = f_A((idx+0.5)-x0)/A0;
            for kn = 1:N
                A(idx(kn)-1,kn) = S1(kn);
                A(idx(kn)  ,kn) = S2(kn)-S1(kn);
                A(idx(kn)+1,kn) = 1-S2(kn);
                %if ~ix(ixM(kn))
                %    p(ixM(kn)) = pvec(idx(kn));
                %    ix(ixM(kn)) = true;
                %end
                %p(ixM(kn)) = pvec(idx(kn)); %!!
                %p(ixM(kn)) = projmat(kp,idx(kn));
                ix(ixM(kn)) = true;
            end
            if delz
                for kn = 1:N
                    p(ixM(kn)) = pvec(idx(kn)); %!!
                end
            end
            W((D*(kp-1)+(1:D)),:) = A; %#
        end
    else
        for kp = 1:n_p
            if delz
                im_rot = imrotate(im_pad,-angles(kp),'bilinear','crop');
                pvec = sum(im_rot,1);
                projmat(kp,:) = pvec;
            end
            t = -angles(kp)/180*pi;
            %R = [cos(t) -sin(t);sin(t) cos(t)];
            R = [sin(t) cos(t)]; %!!
            fprintf('\nAngle %d(%.1f Degree)\n%.1f%% completed.\r',...
                kp,angles(kp),100*kp/n_p)
            xy_rot = R*xy_c+rc+xypad; %!!
            %idx = round(xy_rot(2,:));
            idx = round(xy_rot); %!!
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
            x0 = xy_rot-l/2; %!!
            S1 = f_A((idx-0.5)-x0);
            S2 = f_A((idx+0.5)-x0);
            for kn = 1:N
                W(ixM(kn)-1,kn) = S1(kn);
                W(ixM(kn)  ,kn) = S2(kn)-S1(kn);
                W(ixM(kn)+1,kn) = 1-S2(kn);
                %if ~ix(ixM(kn))
                %    p(ixM(kn)) = pvec(idx(kn));
                %    ix(ixM(kn)) = true;
                %end
                %p(ixM(kn)) = pvec(idx(kn)); %!!
                %p(ixM(kn)) = projmat(kp,idx(kn));
                ix(ixM(kn)) = true;
            end
            if delz
                for kn = 1:N
                    p(ixM(kn)) = pvec(idx(kn)); %!!
                end
            end
        end
    end
else
    fprintf('%s is not a property value of ''option''.\r',opt);
    return
end
t2 = toc;
fprintf('\nTotal elapsed time: %f seconds\n\r',t2-t1);
if delz
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
    % figure('name','Sparsity pattern of matrix W','numbertitle','off');
    % spy(W)
end
nout = nargout;
if nout > 1
    if delz
        varargout{1} = p;
    else
        varargout{1} = ix;
    end
    if nout > 2
        varargout{2} = D;
        if nout > 3
            varargout{3} = projmat;
        end
    end
end
end
