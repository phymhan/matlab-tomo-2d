function fb = create_filter_bank(fsize, fnum)
sigma = linspace(0.05, (fsize-1)/4, fnum);
fb = zeros([fsize fsize 1 fnum]);
for i = 1:fnum
    fb(:,:,1,i) = fspecial('Gaussian', fsize, sigma(i));
end