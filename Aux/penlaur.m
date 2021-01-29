function [Si,Pi,Di,order] = penlaur(E,A)
% [Si,Pi,Di,order]=penlaur(E,A)
% First Laurent coefficients of (s*E-A)^-1;
% (s*E-A)^-1 = ... + Si/s - (Pi + s*Di + ... +s^order Ni) at s = infinity
% order = order of the singularity
% The matrix s*E-A should be invertible.
% Experimental version: troubles when bad conditioning of
% (so*E-A)...)

tests = [-0.7616491, 0.6755538, 1.4739762, 1.1443051, 0.8529776, 0.4529708, 0.7223316, 1.9273332, 0.6380837, -0.8498894];
conditions=0*tests;k=1;
for s0=tests, conditions(k)=cond(s0*E-A);k=k+1;end
[w,k1]=min(conditions);
if w>1e20
    error('penlaur');
end
s0=tests(k1);
J=s0*E-A;
[Se,Pe,~,i1]=projspec(J\E);
[Sa,~,~,~]=projspec(J\A);
order=i1-1;
Si=Se/J;
Pi=(Pe*Sa)/J;
Di=Pi*E*Pi;
if order==0 
    Di=0*Di;
end
    
end