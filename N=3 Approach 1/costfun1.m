function [J, grad, ResVec]= costfun1(dvar,N,odefunvec,optODE,x0)
%dvar,N,odefunvec,optODE,ts,x0
%% Part1. Objective J for piecewise constant control
nx = 2;
e1  = 18000;
e2  = 30000;
k10 = 0.535*10^11;
k20 = 0.461*10^18;
% alpha = e2/e1;
% c = k20/(k10)^alpha;
% variableName = sprintfc('x%d(t)',1:nx);
% syms(variableName{:})
% u = sym('u');
% alpha = sym('alpha');
% c = sym('c');
alpha = e2/e1;
c = k20/(k10)^(e2/e1);
% f(1,1) = -u*x1;
% f(2,1) =  u*x1-c*(u^alpha)*x2;

z0 = x0;
ts = zeros(1,N+1);
for i = 1 : N
    ts(i+1) = sum(dvar(1+N:N+i));
end
ResVec = [];
for i = 1 : N
    u1 = dvar(i);
    odefun = odefunvec{i};
    odefun = @(t,Y) odefun(t,Y,u1);
    [t,zout] = ode45(odefun,[ts(i) ts(i+1)],z0,optODE);
    if i ~= N
        z0 = zout(end,:);
        % z0(2*i+7)=z0(2*i+7)+...
        %   (subs(f(1),{u,x1,x2},{dvar(i),  z0(1),z0(2)})...
        %   -subs(f(1),{u,x1,x2},{dvar(i+1),z0(1),z0(2)}));
        z0(2*i+7)=z0(2*i+7)+...
            -  dvar(i)*z0(1)...
            -(-dvar(i+1)*z0(1));
        % f(2,1) =  u*x1-c*(u^alpha)*x2;
%         z0(2*i+8)=z0(2*i+8)+...
%             (subs(f(2),{u,x1,x2,alpha,c},{dvar(i),  z0(1),z0(2),e2/e1,k20/(k10)^(e2/e1)})...
%             -subs(f(2),{u,x1,x2,alpha,c},{dvar(i+1),z0(1),z0(2),e2/e1,k20/(k10)^(e2/e1)}));
        z0(2*i+8)=z0(2*i+8)+...
            dvar(i)*z0(1) - c*(dvar(i)^alpha)*z0(2) - ...
           (dvar(i+1)*z0(1) - c*(dvar(i+1)^alpha)*z0(2));
        ResVec = [ResVec;t zout];
        % t(end) z0];
    else % i == N 
        zfinal = zout(end,:);
        zfinal(13) = -dvar(3)*zfinal(1);
        zfinal(14) = dvar(3)*zfinal(1) - c*(dvar(3)^alpha)*zfinal(2);
        zout(end,:) = zfinal;
        ResVec = [ ResVec;
                   t zout];
    end
end
J = -zfinal(2);
%% Part2. Gradient grad by Sensitivity Equation
if nargout > 1 % Gradient Required
    % nargout: Number of function output arguments
    % https://www.mathworks.com/help/matlab/ref/nargout.html
    Jt1 = -zfinal(10);
    Jt2 = -zfinal(12);
    % f(2,1) =  u*x1-c*(u^alpha)*x2;
    Jt3 = - zfinal(14);
    % double(subs(-f(2),{u,x1,x2,alpha,c},{dvar(3),zfinal(1),zfinal(2),e2/e1,k20/(k10)^(e2/e1)}));

    Jdelt1 = 1*Jt1 + 1*Jt2 + 1*Jt3;
    Jdelt2 = 0*Jt1 + 1*Jt2 + 1*Jt3;
    Jdelt3 = 0*Jt1 + 0*Jt2 + 1*Jt3;

    grad(1:3) = -zfinal(4:2:8);
    grad(4:6) = [Jdelt1, Jdelt2, Jdelt3];
end

end % Func end