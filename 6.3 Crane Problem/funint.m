function out = funint(dvar,m,optODE,x0_orig)
% myfunint2:
% integration function for piecewise constant control
%% Part1. Pars Settings
z0 = x0_orig;
u1vec = dvar(      1:  m.N);
u2vec = dvar(  m.N+1:2*m.N);
u3vec = dvar(2*m.N+1:3*m.N);
for ks = 1 : m.N % N
    [~,z] = ode45(@(t,x)dyneqn(t,x,u1vec(ks),u2vec(ks),u3vec(ks)),...
        [m.ts(ks) m.ts(ks+1)],z0,optODE);
    z0 = z(end,:)';
end
out = z0; 
end