clear all 
N = 3; % # of time stages 
nx = 2;
nu = 2;
nx_t = nx + nx*(nu*N); % 2+ 2*(2*3) = 14
e1  = 18000;
e2  = 30000;
k10 = 0.535*10^11;
k20 = 0.461*10^18;
beta1 = 0.53;
beta2 = 0.43;
% alpha = e2/e1;
% c = k20/(k10)^alpha;
x0 = [beta1 beta2 zeros(1,nx_t-nx)]';

%% Part 3.Build sensitivity function
variableName = sprintfc('x%d(t)',1:nx_t);
syms(variableName{:})
syms u alpha c % theta3 
vars = str2sym(variableName);
f(1,1) = -u*x1;
f(2,1) =  u*x1-c*(u^alpha)*x2;
pfpx = Jacob_fun(f,nx); 
pfpu =  jacobian(f, u);
F = sym(zeros(nx_t,1));
F(1:2,:) = f; 
S = sym(zeros(nx,nu*N));
for i = 1 : nx*(nu*N) % 12
    S(i) = str2sym(sprintfc('x%d(t)', i+nx));
end
odefunvec = cell(1,N);
for i = 1 : N
    pupw = zeros(1,6);
    pupw(1,i) = 1;
    res = pfpx*S + pfpu*pupw;
    %     if i ~= N
    for j = nx+1: nx_t
        F(j,1) = res(j-nx);
    end
    %     else
    %         res(:,end) = res(:,end) + [f(1,1)/theta3;f(2,1)/theta3];
    %         for j = nx+1: nx_t
    %             F(j,1) = res(j-nx);
%         end 
%     end
    F = subs(F,[alpha,c],[e2/e1,k20/(k10)^(e2/e1)]);
    odefunvec{i} = odeFunction(F,vars,u);
end
save('functionHandle.mat','odefunvec')