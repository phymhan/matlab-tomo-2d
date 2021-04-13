function imp = impiece(im1,im2,ix)
% For each layer, im1(ix) = im2(ix). Be sure im1 and im2 have the same
% scale.

for kc = 1:size(im1,3)
    m1 = im1(:,:,kc);
    m2 = im2(:,:,kc);
    m1(ix) = m2(ix);
    im1(:,:,kc) = m1;
end
if nargout < 1
    caller = inputname(1);
    assignin('caller',caller,im1);
end
imp = im1;
end
