function  [J,grad] = costfunTS(dvar, m, odefunvec, optODE, x0_red, epsilon, rho)            
%% Part1. Objective J for piecewise constant control
% f = myfunint1(dvar,N,odefunvec,optODE,ts,x0_red,epsilon);
z0 = x0_red;
for i = 1 : m.N % N
    u1 = dvar(i+0*m.N);
    u2 = dvar(i+1*m.N);
    u3 = dvar(i+2*m.N);
    odefun = odefunvec.augfun{i};
    odefun = @(t,Y) odefun(t,Y,u1,u2,u3,epsilon);
    [~,z_red] = ode45(odefun,[m.ts(i) m.ts(i+1)],z0,optODE);
    z_last = z_red(end,:)';
    if i ~= m.N
        z0(1:m.nx) = z_last(1:m.nx); %#1
        z0(m.nx+1:m.nx+1+m.nx*i-1) = z_last(m.nx+1:m.nx+1+m.nx*i-1); %#2
        z0(m.nx+1+m.nx*i:m.nx+1+m.nx*(i+1)-1) = 0; %#3
        z0(m.nx+1+m.nx*(1*i+1):m.nx+1+m.nx*(2*i+1)-1) = z_last(m.nx+1+m.nx*i:m.nx+1+m.nx*(2*i)-1); %#4
        z0(m.nx+1+m.nx*(2*i+1):m.nx+1+m.nx*(2*i+2)-1) = 0; %#5
        z0(m.nx+1+m.nx*(2*i+2):m.nx+1+m.nx*(3*i+2)-1) = z_last(m.nx+1+m.nx*(2*i):m.nx+1+m.nx*(3*i)-1); %#6
        z0(m.nx+1+m.nx*(3*i+2):m.nx+1+m.nx*(3*i+3)-1) = 0; %#7
    else
        break 
    end
end

J = z_last(7)+rho*sum(z_last(8:9));
%% Part2. Gradient grad by Sensitivity Equation
grad = zeros(1,m.nu*m.N)';
if nargout > 1 % Gradient Required 
    for i = 1 : m.nu*m.N
        grad(i) = z_last(16+(i-1)*m.nx) + ... 
         rho*sum(z_last(16+(i-1)*m.nx+1: 16+(i-1)*m.nx+2));
    end 
else 
    grad = [];
end

end % Func end