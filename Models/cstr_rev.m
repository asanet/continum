function dy = cstr_rev(~,y,par)

% %Variables
x1 = y(1); x2 = y(2); x3 = y(3); x4 = y(4); x5 = y(5);

% %Parameters
Da  = par(1);    g  = par(2);    beta = par(3); 
c   = par(4);    H1 = par(5);    H2   = par(6);
Keq = par(7);  

% %Equations
dy(1,1) = 1 - x1 - x5;
dy(2,1) = -x2 - x4 + x5;
dy(3,1) = 1 - x3 + H1*x5 +H2*x4 - beta*(x3-c);
dy(4,1) = x4 - Da*exp(-g/x3)*x2;
dy(5,1) = x1*Keq - x2;


% %lb = [0,0,0,0.5,0,0,0]
% %ub = [10,100,10,10,10,100,100]
% %p = [9.9644, 30.4871,  3.2207,  4.9193,  6.7836, 59.4390, 54.7310]

end