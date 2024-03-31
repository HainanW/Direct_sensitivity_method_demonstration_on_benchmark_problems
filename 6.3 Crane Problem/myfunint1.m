function out = myfunint1(dvar,N,odefunvec,optODE,ts,x0,epsilon)
%myfunint1:integration function for piecewise constant control
%% Part1. Pars Settings
z0 = x0;
for i = 1 : N % N
    u1 = dvar(i);
    u2 = dvar(i+N);
    odefun = odefunvec{i};
    odefun = @(t,Y) odefun(t,Y,u1,u2,epsilon);
    [~,z_red] = ode45(odefun,[ts(i) ts(i+1)],z0,optODE);
    z_last = z_red(end,:)';
    if i~= N 
        z0(1:3) = z_last(1:3);
        z0(4:4+3*i-1) = z_last(4:4+3*i-1); 
        z0(4+3*i:4+3*(i+1)-1) = 0; 
        z0(4+3*(1*i+1):4+3*(2*i+1)-1) = z_last(4+3*i:4+3*(2*i)-1); 
        z0(4+3*(2*i+1):4+3*(2*i+2)-1) = 0; 
    else 
        break
    end 
end
out = z_last; 
end