%TEST 04: 3D Reconstruction
clear
close all
clc
fprintf('TEST 05: 3D Reconstruction\nUsing a downloaded toolbox ''tomobox.''\r');
%% Dimensions to use in simulation
r1_max   = 17;
N        = 2*r1_max+1;      % 3D image side length
dims     = N*[1,1,1];
u_max    = 23;              % Projection side length
nuv      = [1,1];           % Number of subpixels, see buildSystemMatrix
vpRatio  = 1;               % Voxel-to-pixel ratio, see buildSystemMatrix
num_proj = 30;              % Number of projections to reconstruct from
rnl      = 0.01;            % Relative noise level

%% Set up 3D test image
X_true = phantom3d('Modified Shepp-Logan',N);
x_true = X_true(:);

%% Choose a number of random projection directions

% Choose x,y,z as Gaussian triple and normalize. Since Gaussians are
% rotation-symmetric, the directions obtained are samples from the uniform
% distribution over the unit sphere.
v_list     = randn(num_proj,3);
v_listnorm = sqrt(sum(v_list.^2,2));
v_list     = v_list./repmat(v_listnorm,1,3);

%% Set up the parallel beam system matrix
[A,p_all] = buildSystemMatrix(r1_max,u_max,v_list,nuv,vpRatio);

%% Compute projections
b_orig = A*x_true;

%% Add Gaussian noise
e = getNoise(rnl,b_orig);
b = b_orig + e;
B = reshape(b,2*u_max+1,2*u_max+1,num_proj);

%% Reconstruct
tol = 1e-6;
maxit = 200;
x_sol = lsqr(A,b,tol,maxit);
X_sol = reshape(x_sol,dims);
%% Display

figure('name','Original')
plotLayers(X_true)
suptitle('Original')

figure('name','Projections')
plotLayers(B)
suptitle('Projections')

figure('name','Reconstruction')
plotLayers(X_sol)
suptitle(['Reconstruction, iteration ',num2str(maxit)])

%%
%myslicer
figure('name','MySlicer')
X_true = myimadj(X_true);
%T = [1 0 0 0;0 1 0 0;0 0 2.5 0];
myslicer(X_true);
colormap gray;
axis off;
