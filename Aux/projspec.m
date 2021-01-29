function [S,P,D,index]=projspec(A,tol)
%[S,P,D,index]=projspec(A)
%Spectral characteristics of A at 0
%S = reduced resolvent at 0 (S=-Drazin_inverse(A))
%P = spectral projection at 0
%D = Nilpotent operator at 0
%index = index of the 0 eigenvalue

if nargin < 2
    tol = 1e-6;
end

n = size(A,1);
if norm(A,1) < eps*n^2
    P = eye(n);
    D = A;
    S = 0*P;
    index=1;
    return %% mod
end

if rcond(A) > tol && ~isempty(A)
    S = inv(A);
    P = 0*eye(n);
    D = P;
    index = 0;
    return
end

index = 1;
[B,C,dim] = fullrf(A);

if dim == 0  
    P = eye(n);
    S = 0*P;
    D = A;
    return
end

Ck = C;
Bk = B;
tst = rcond(Ck*Bk);
if size(Ck,1)==1
    tst = norm(Ck*Bk,1);
end
if tst > tol
%     M = inv(C*B);
    M = C*B;
%     P = eye(n)-B*M*C;
    P = eye(n)-(B/M)*C;
%     S = B*M*M*C;
    S = ((B/M)/M)*C;
    D = 0*A;
    return
end

for k = 2:n
    [B,C,dim] = fullrf(C*B);
    if dim == 0
        P = eye(n);
        S = 0*P;
        D = A;
        index=k;
        return
    end
    Bk = Bk*B;
    Ck = C*Ck;
    index = k;
    if norm(C*B)*rcond(C*B) > tol
%         M = inv((C*B)^index);
        M = (C*B)^index;
%         P = eye(n) - Bk*M*Ck;
        P = eye(n) - (Bk/M)*Ck;
%         S = Bk*M*inv(C*B)*Ck;
        S = ((Bk/M)/(C*B))*Ck;
        D = 0.5*(A*P+P*A);
        return
    end
end
P = eye(n);
S = 0*P;
D = A;
index = k;
end