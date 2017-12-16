function x = solver_art(A,b,n_it)
% N = size(A,2);
% x = zeros(N,1);
% ix = (A~=0)';
% Z = sum(ix,2)*N;
% lambda = 1;
% cnt = 0;
% for ki = 1:n_it
%     cnt = cnt+1;
%     if cnt >= n_it/10
%         fprintf('\nIteration %d\r',ki);
%         cnt = 0;
%     end
%     x = x-lambda*(ix*(A*x-b))./Z;
% end

[M,N] = size(A);
x = zeros(N,1);
A_sq = sum(A.*A,2);
lambda = 1;
cnt = 0;
for ki = 1:n_it
    cnt = cnt+1;
    if cnt >= n_it/10
        fprintf('\nIteration %d\r',ki);
        cnt = 0;
    end
    ii = rem(ki-1,M)+1;
    x = x-lambda*((A(ii,:)*x-b(ii))/A_sq(ii))*A(ii,:)';
end
end
