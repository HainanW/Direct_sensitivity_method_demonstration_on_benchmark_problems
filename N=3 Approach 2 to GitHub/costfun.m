function [J, grad]= costfun(dvar,m, odefunvec, optODE, x0)
% Function Input: dvar, m, odefunvec, optODE,x0
% dvar(  1:  N): u1, u2, u3 
% dvar(N+1:2*N): theta1, theta2, theta3 
%% Part1. Objective J for piecewise constant control

tspan = zeros(1,m.N);
for i = 1 : m.N
    tspan(i+1) = sum(dvar(1+m.N:m.N+i));
end

theta1 = dvar(m.N+1);
theta2 = dvar(m.N+2);
theta3 = dvar(m.N+3);

z0 = x0;

for i = 1 : m.N 
    u1 = dvar(i);
    odefun = odefunvec{i};
    odefun = @(t,Y) odefun(t, Y, u1, theta1, theta2, theta3);
    [~,z] = ode45(odefun,[tspan(i) tspan(i+1)],z0,optODE);
    z0= z(end,:)';
end

J = -z0(2);
%% Part2. Gradient grad by Sensitivity Equation
if nargout > 1 % Gradient Required
    grad = - z0(4:2:end);
else 
    grad = [];
end

end % Func end