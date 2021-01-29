function dydt = cstr_simple(~, y, par)

% Variables
u = y(1); v = y(2);

% Parameters
Da = par(1); B = par(2); beta = par(3);

% Equations
dydt(1,1) = 1 - u.*(1 + Da*exp(v));
dydt(2,1) = B*Da*u.*exp(v) -v*(1+beta);

end