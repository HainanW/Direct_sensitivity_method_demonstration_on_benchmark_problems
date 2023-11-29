clear all 
N = 3; % # of time stages 
nx = 2;
nu = 2;
nx_t = nx + nx*(nu*N); % 2+ 2*(2*6) = 26 
e1  = 18000;
e2  = 30000;
k10 = 0.535*10^11;
k20 = 0.461*10^18;
beta1 = 0.53;
beta2 = 0.43;
alpha = e2/e1;
c = k20/(k10)^alpha;
x0 = [beta1 beta2 zeros(1,nx_t-nx)]';

%% Part 3.Build sensitivity function
variableName = sprintfc('x%d(t)',1:nx_t);
syms(variableName{:})
syms u v
vars = str2sym(variableName);
f(1,1) = -v*(u*(x1));
f(2,1) =  v*(u*x1-c*(u^alpha)*x2);
pfpx = Jacob_fun(f,nx); % pfpx = jacobian(f,sym('x',[1, nx]))
pfpu = jacobian(f,[u v]);
F = sym(zeros(nx_t,1));
F(1) = f(1,1);
F(2) = f(2,1);
S = sym(zeros(nx,nu*N));
for i = 1 : nx*(nu*N) % 24
    S(i) = str2sym(sprintfc('x%d(t)', i+nx));
end
odefunvec = cell(1,N);
for i = 1 : N
    pupw = zeros(nu,nu*N); % 2x12
    pupw(1,i) = 1;
    pupw(2,i+N) = 1;
    res = pfpx*S + pfpu*pupw; % f(2) = subs(f(2),[a,b,c],[1,0,1]);
    for j = nx+1: nx_t % Be advised, size f change
        F(j,1) = res(j-nx);
    end
    odefunvec{i} = odeFunction(F,vars,u,v);
end
save('functionHandle.mat','odefunvec')