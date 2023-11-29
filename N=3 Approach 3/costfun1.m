function [J,grad, zvec]= costfun1(dvar,N,odefunvec,optODE,ts,x0)
%dvar,N,odefunvec,optODE,ts,x0
%% Part1. Objective J for piecewise constant control
z0 = x0;
%% Part2. Solve ODE with Sensitivity
zvec = [];
for i = 1 : N 
    u1 = dvar(i);
    u2 = dvar(i+N);
    odefun = odefunvec{i};
    odefun = @(t,Y) odefun(t,Y,u1,u2);
    [t,z] = ode45(odefun,[ts(i) ts(i+1)],z0,optODE);
    z0= z(end,:)';
    z_rec = [t z];
    zvec = [zvec; z_rec];
end

J = -z0(2);
% Part2. Gradient grad by Sensitivity Equation 
% if nargout > 1 % Gradient Required
grad = -z0(4:2:14);
% end

end % Func end 