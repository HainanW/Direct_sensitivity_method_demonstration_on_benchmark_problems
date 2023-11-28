function out = myfunint1(dvar,N,odefunvec,optODE,ts,x0)
%myfunint1:integration function for piecewise constant control
%% Part1. Pars Settings
z0 = x0;
%% Part2. Solve ODE with Sensitivity
for i = 1 : N 
    u1 = dvar(i);
    u2 = dvar(i+N);
    odefun = odefunvec{i};
    odefun = @(t,Y) odefun(t,Y,u1,u2);
    [~,z] = ode45(odefun,[ts(i) ts(i+1)],z0,optODE);
    z0= z(end,:)';
    
end
out = z0; 
end