clear all 
%% 1 Equation and Parameter Definition 
N = 3; % # of time stages 
nx = 2;
nu = 1;
nx_t = nx + nx*(nu+1)*N; % 2+ 2*(2*3) = 14
e1  = 18000;
e2  = 30000;
k10 = 0.535*10^11;
k20 = 0.461*10^18;
beta1 = 0.53;
beta2 = 0.43;
alpha = e2/e1;
c = k20/(k10)^alpha;
x0 = [beta1 beta2 zeros(1,nx_t-nx)]';

%% Part 2.Build sensitivity function
x_Name = sprintfc('x%d(t)', 1:nx_t);
theta_Name = sprintfc('theta%d', 1:N);
syms(x_Name{:})
syms(theta_Name{:})
syms u
xvars = str2sym(x_Name);
thetavars = str2sym(theta_Name);

f_orig = sym(zeros(nx, 1)); 
f_orig(1,1) = -u*x1;
f_orig(2,1) =  u*x1 - c*(u^alpha)*x2;

% f_over_theta = sym(zeros(nx, N));
f_over_theta = f_orig./thetavars; 


pfpx = Jacob_fun(f_orig,nx); 
pfpu =  jacobian(f_orig, u);

Fvec = sym(zeros(nx_t,1));
Fvec(1:nx) = f_orig; 
S_u = sym(zeros(nx,N));
S_theta = sym(zeros(nx,N));
for i = 1 : nx*N % 12
    S_u(i) = str2sym(sprintfc('x%d(t)', i+nx));
    S_theta(i) = str2sym(sprintfc('x%d(t)', i + nx + nx*N));
end
odefunvec = cell(1,N);

pfpx_S_u = pfpx*S_u;
pfpx_S_theta = pfpx*S_theta;
for i = 1 : N
    pupw = zeros(1,N);
    pupw(1,i) = 1;
    res_u = pfpx_S_u + pfpu*pupw;

    f_over_theta_chi = sym(zeros(size(pfpx_S_theta)));
    f_over_theta_chi(:,i) = f_over_theta(:,i);
    res_theta = pfpx_S_theta + f_over_theta_chi; 
    for j = nx+1: nx+nx*N % 3-8
        Fvec(j,1) = res_u(j-nx);
    end
    for k = nx+nx*N+1:nx_t
        Fvec(k,1) = res_theta(k-nx-nx*N);
    end 
    odefunvec{i} = odeFunction(Fvec,xvars,u,theta1,theta2,theta3);
end
save('functionHandle.mat','odefunvec')