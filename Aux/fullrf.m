function [Q,M,rk] = fullrf(A)
%[Q,M,rk]=fullrf(A)
%Full rank factorization : A=Q.M
%with range(Q)=range(A) and ker(M)=ker(A),
%Q full column rank , M full row rank
% rk = rank(A) = #columns(Q) = #rows(M)

na1 = norm(A,1);
if na1 < 1e-10
    Q = [];
    M = [];
    rk = 0;
    return
end

[U,s,V] = svd(A);
s(s<sqrt(eps)) = 0;
sq = sqrt(s);
Q = U*sq;
M = sq*V';
rk = rank(s);
Q = Q(:,1:rk);
M = M(1:rk,:);
    
end