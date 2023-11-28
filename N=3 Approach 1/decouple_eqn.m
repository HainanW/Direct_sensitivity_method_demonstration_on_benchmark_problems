clear all 


z0 = [0, 0];
optODE = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
% Random control in lb/ub
% u = lb + (ub-lb)*rand
u = 0.1 + (2-0.1)*rand(1);
[~,zout] = ode45(@(t,x)decoupledeqn(t,x,u),[0, 8],z0,optODE)


function  dyn = decoupledeqn(t,x,u)
e1  = 18000;
e2  = 30000;
k10 = 0.535*10^11;
k20 = 0.461*10^18;
alpha = e2/e1;
c = k20/(k10)^(e2/e1);
x(1) = -u*x(1);
x(2) = u*x(1) - c*u^alpha*x(2);

dyn = [x(1);x(2)];
end