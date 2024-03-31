function [c,ceq] = consfunTS(dvar, m, odefunvec, optODE, x0_orig)
% [c,ceq,cgd,ceqgd]=confun1
z0 = x0_orig;
for i = 1 : m.N % N
    u1 = dvar(i+0*m.N);
    u2 = dvar(i+1*m.N);
    u3 = dvar(i+2*m.N);
    odefun = odefunvec.statefun;
    odefun = @(t,Y) odefun(t,Y,u1,u2,u3);
    [~,z_red] = ode45(odefun,[m.ts(i) m.ts(i+1)],z0,optODE);
    z0 = z_red(end,:)';
end
z_last = z0;
%% Part2. Gradient grad by Sensitivity Equation
c = [];
ceq(1) = z_last(1) - 10;
ceq(2) = z_last(2) - 14;
ceq(3) = z_last(3) - 0.0;
ceq(4) = z_last(4) - 2.5; 
ceq(5) = z_last(5);
ceq(6) = z_last(6);
ceq(7) = sum(dvar(2*m.N+1:end)) - 9;
end 