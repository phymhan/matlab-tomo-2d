%TEST 12: Radon Transform and Hough Transform
clear
clc
fprintf('TEST 12: Radon Transform and Hough Transform\r');
% im = zeros(128);
% im(31:32,63:64) = 255;
% im(63:64,96:97) = 255;
% im(96:97,63:64) = 255;
% im(63:64,31:32) = 255;
% figure
% imshow(uint8(im));
% a = 0:1:179;
% [pm,~] = tomoproj2d(im,a);
% figure
% imshow(uint8(imscale(pm)));

% figure('color',[1 1 1])
% axes('xlim',[0 2*pi],'ylim',[0 90])
% box on
% xlabel('\phi')
% ylabel('r')
% line(pi,32,'marker','o','markeredgecolor',[1 0 0])
% line(pi/2,32,'marker','o','markeredgecolor',[1 0 0])
% line(0,32,'marker','o','markeredgecolor',[1 0 0])
% line(3*pi/2,32,'marker','o','markeredgecolor',[1 0 0])

im = zeros(128);
im(83:84,63:64) = 255;
% figure
% imshow(uint8(im));
im(83:84,31:32) = 255;
% figure
% imshow(uint8(im));
im(83:84,01:02) = 255;
% figure
% imshow(uint8(im));
im(83:84,96:97) = 255;
figure
imshow(uint8(im));
