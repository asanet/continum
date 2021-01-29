function [x,fval,eflag,out,jac] = hfsolve(fun,x0,opt,nP,debug,varargin)

n = length(x0);

if nargin < 2 || isempty(fun) || isempty(x0)
    error('hfsolve needs at least a function and an initial guess :)')
end

if nargin < 3 || isempty(opt)
    opt = optimoptions(@fsolve); 
end

if nargin < 4 || isempty(nP)
    nP = 25; 
end
    
if nargin < 5 || isempty(debug)
    debug = false;
end

if nargin < 6
    hmodel = @(x,alpha)  alpha*fun(x) + (1-alpha)*(x - ones(n,1));
    numF = length(fun(x0));
else
    hmodel = @(x,alpha,par)  alpha*fun(x,par) + (1-alpha)*(x - ones(n,1));
    numF = length(feval(fun,x0,varargin{:}));
    
end

if numF > n
    error('You have more equations than variables. System must be square :(')
elseif numF < n 
    error('You have less equations than variables. System must be square :(')
end

if isempty(opt)
    opt = optimoptions(@fsolve);
end

if nP < 1 || isempty(nP)
    nP = 25;
end

    x = x0(:);
    optD = opt.Display;
    opt.Display = 'none';
    for alpha = 0:1/nP:1-1/nP
        [x,fval] = fsolve(hmodel,x,opt,alpha,varargin{:});
        if debug
            fprintf('alpha = %2.4f\tNorm(fval) = %2.8e\n',alpha, norm(fval))
        end
    end
    opt.Display = optD;
    [x,fval,eflag,out,jac] = fsolve(hmodel,x,opt,1,varargin{:});

end
