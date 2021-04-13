function [I, cnt] = flow_unwarp(J, flow, inpaint)
if ~exist('inpaint', 'var') || isempty(inpaint)
    inpaint = true;
end

sz = [size(J,1) size(J,2)];
if sz(1)~=size(flow,1) || sz(2)~=size(flow,2)
    flow = imresize(flow, sz, 'bicubic');
end
I = nan(size(J));
cnt = zeros(sz);

ind = 1:prod(sz);
[y, x] = ind2sub(sz, ind');
x_ = x + reshape(flow(:,:,1), [], 1);
y_ = y + reshape(flow(:,:,2), [], 1);
x1 = round(x_);
y1 = round(y_);
for i = ind
    if x1(i)>=1 && x1(i)<=sz(2) && y1(i)>=1 && y1(i)<=sz(1)
        if cnt(y1(i),x1(i)) == 0
            I(y1(i),x1(i),:) = J(y(i),x(i),:);
        else
            % average:
            %(cnt(y1(i),x1(i))*I(y1(i),x1(i))+J(y(i),x(i),:))/(cnt(y1(i),x1(i))+1);
            % random pick:
            if rand < 1/(cnt(y1(i),x1(i))+1)
                I(y1(i),x1(i),:) = J(y(i),x(i),:);
            end
        end
        cnt(y1(i),x1(i)) = cnt(y1(i),x1(i))+1;
    end
end

if inpaint
    se = strel('disk', 5);
    ix = cnt > 0;
    ix = imerode(imdilate(ix, se), se);
    I(~ix) = 0;
    I = inpaintn(I);
    I(~ix) = nan;
end
I = cast(I, 'like', J);
