function [c, cn] = labudde(A,B)
% LABUDDE   Calculates the characteristc polynomial of a Matrix A or a matrix
%           pair (A,B) using La Budde's method.
%
%           Example: [c,cn] = labudde(A,B);
%                   
%           where:  A and B are (n x n) matrices.
%                   c is a vector of length (n+1) with the coefficients of
%
%                   p(x) = c(1)*x^n + c(2)*x(n-1) + ... + c(n)*x + c(n+1)
%     
%                   cn is the normalized polynomial coefficients: cn = c/c(1)
%
% Developed by Ataide Neto as an extension of the original method.
% Chemical Engineering Program - COPPE/UFRJ (May 2017)
%  
% Reference: Rehman, R. and Ipsen, I.C.F. "La Budde's method for computing
%           characteristic polynomials" 2011.

if nargin < 1
    error('Not enough input arguments. :(')
elseif nargin == 1
    [n,nn] = size(A);
    if n ~= nn
        error('Matrix must be square. :(')
    end
    T = eye(n);
    H = hess(A);
elseif nargin == 2
    [nA,mA] = size(A);
    [nB,mB] = size(B);
    if nA ~= mA || nB ~= mB
        error('Matrix must be square. :(')
    elseif nA ~= nB
        error('Matrix dimensions must agree. :(')
    end
    n = nA;
    [H,T] = hess(A,B);
else
    warning('Extra inputs will be neglected. :|')
end

% %Some definitions 
t = diag(T,0);
a = diag(H,0);
b = diag(H,-1);

% %Alocation and initialization
c = zeros(n+1,n+1);

c(1,1) = 1;
c(1,2) = t(1);
c(2,2) = -a(1);

% %Compute coefficients recursively
for i = 3:n+1
    
    c(1,i) = prod(t(1:i-1));
    for j = 2:i-1
        
        s1 = 0;
        for m = 1:j-2
            s1 = s1 + H(i-m-1,i-1)*prod(b(i-m-1:i-2))*c(j-m-1,i-m-1) ;           
        end
        
        s2 = 0;
        for m = 1:j-1
            s2 = s2 + T(i-m-1,i-1)*prod(b(i-m-1:i-2))*c(j-m,i-m-1);
        end
        
        c(j,i) = t(i-1)*c(j,i-1) - a(i-1)*c(j-1,i-1) - s1 + s2;       
    end
    
    s = 0;
    for m = 1:i-2
        s = s + H(i-m-1,i-1)*prod(b(i-m-1:i-2))*c(i-m-1,i-m-1);
    end
    
    c(i,i) = -a(i-1)*c(i-1,i-1) - s;
end

% %Extract the core result
c = c(:,n+1);

% %Trim null coefficients
while true
    if abs(c(1)) < n*sqrt(eps)
        c(1) = [];
    else
        break
    end    
end

% %Normalize polynomial
cn = c/(c(1));

end