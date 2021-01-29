function [xo,Ot,nS] = argeSwarm(S,x0,Lb,Ub,nger,npas,ObjectiveLimit,MaxFunEvals,debug,varargin)

n = length(x0);
problem = -1;
c1 = 1;
c2 = 1;

w = 0.9;
tw = (w-4e-3)/nger;
Ot = feval(S,x0,varargin{:})*problem;

p = zeros(n,npas);
v = zeros(n,npas);
x = zeros(n,npas);

del = Ub-Lb;
 
for j = 1:npas
    x(:,j) = Lb + rand(n,1).*del;
    v(:,j) = 0.1*rand(n,1).*del;
    p(:,j) = x(:,j);
end	

x(:,1) = x0;
p(:,1) = x0;
nS = 1;
ipg = 1;

if debug
    fprintf('\n F-Count               Best value\n');
end

y = zeros(npas,1);
y(1) = Ot;
for j = 2:npas
    y(j) = feval(S,x(:,j),varargin{:})*problem;
    nS = nS + 1;
    if y(j) > Ot
        Ot = y(j);
        ipg = j;
   end
end

if debug
    fprintf(' %7.0d %24.6e\n',[nS,abs(Ot)])
end

xo = p(:,ipg);

for ig = 1:nger
    rnd1 = ones(n,1)*rand(1,npas);
    rnd2 = ones(n,1)*rand(1,npas);
    v = w*v + c1*rnd1.*(p-x) + c2*rnd2.*( p(:,ipg)*ones(1,npas) - x );
    x = x + v;
    
    for i = 1:n
        j = find(x(i,:) > Ub(i));
        if ~isempty(j)
            x(i,j) = Ub(i);
            v(i,j) = 0;
        end

        j = find(x(i,:) < Lb(i));
        if ~isempty(j)
            x(i,j) = Lb(i);
            v(i,j) = 0;
        end
    end

    for j = 1:npas
        val = feval(S,x(:,j),varargin{:})*problem;
        nS = nS + 1;

        if val > y(j)
            y(j) = val;
            p(:,j) = x(:,j);
            if val > Ot
                Ot = val;
                ipg = j;
            end
        end
    end
  
    if debug
        fprintf(' %7.0d %24.6e\n',[nS,abs(Ot)])
    end
    
    xo = p(:,ipg);

    w = w - tw;        
    
    if Ot*problem <= ObjectiveLimit || nS > MaxFunEvals
        break
    end
end

Ot = Ot*problem;

end