function [W,p] = combWeightProj(W,ix,projmat)
%COMBWEIGHTPROJ   Combine weight matrix and projection matrix and generate
%   projection vector p, which is used for solving linear equations. The
%   output argument W will not contain all zero rows.
%
%   Phymhan
%   25-Aug-2013 13:40:36

%M = length(ix);
[n_p,D] = size(projmat);
if n_p*D < size(W,1)
    offset = -1;
    D = D+2;
else
    offset = 0;
end
n_eq = length(find(ix));
p = zeros(n_eq,1);
cnt = 0;
for kp = 1:n_p
    for k = 1:D
        if ix((kp-1)*D+k)
            cnt = cnt+1;
            p(cnt) = projmat(kp,k+offset);
        end
    end
end
try
    W(~ix,:) = [];
catch expr
    fprintf([expr.message '\nConverting to sparse...\r'])
    W = sparse(W);
    W(~ix,:) = [];
end
end
