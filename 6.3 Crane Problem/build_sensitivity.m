clear all

%% 1. Pars Settings
N = 10; % # of time stages
nx = 9;
nu = 3;
nx_t = nx+nx*(nu*N); % 9+ 9*(3*10) = 279
allvariableName = sprintfc('x%d(t)',1:nx_t); % x1~x279
syms(allvariableName{:})
syms u1 u2 u3 epsilon
allxvars = str2sym(allvariableName);

%% 2. Constraint Related
G1 = formula(x4^2-2.5^2);
G2 = formula(x5^2-1.0^2);

P1G1 = 0.5*(G1+(G1^2+4*epsilon^2)^(1/2));
P1G2 = 0.5*(G2+(G2^2+4*epsilon^2)^(1/2));

%% 3. Calculate Sensitivity eqn

% Faux : Auxilary system generated 
% Forig: Original Function Vector 
% Fun  : Function stored in Function handle in each stage

Forig = sym(zeros(nx,1));
Forig(1,1) = u3*x4;
Forig(2,1) = u3*x5;
Forig(3,1) = u3*x6;
Forig(4,1) = u3*(u1+17.2656*x3);
Forig(5,1) = u3*(u2);
Forig(6,1) = u3*(-(u1+27.0756*x3+2*x5*x6)/x2);
Forig(7,1) = u3*(x3^2 + x6^2);
% Forig(8-9,1) are constraints 
Forig(8,1) = u3*P1G1;
Forig(9,1) = u3*P1G2;


pfpx = Jacob_fun(Forig,nx); % pfpx = jacobian(f,sprintfc('x%d(t)',1:3))
pfpu = jacobian(Forig,[u1 u2 u3]);
S = sym(zeros(nx,nu*N)); % 9,30  
for i = 1 : nx*(nu*N) % 9*(3*10) = 270
    S(i) = str2sym(sprintfc('x%d(t)', i+nx));
end
pfpx_S = pfpx*S;
%% 4. Build Function Handle
Faux = sym(zeros(nx_t,1));
Faux(1:nx) = Forig(1:nx);
augfuncell = cell(1,N);
for i = 1 : N
    tic 
    Fun = sym(zeros(nx*nu*i+nx,1));
    pupw = zeros(nu,nu*N); % 3x30(90)
    pupw(1,i) = 1;
    pupw(2,i+N) = 1;
    pupw(3,i+2*N) = 1; 
    res = pfpx_S + pfpu*pupw;
    for j = nx+1 : nx_t % 10: 279
        Faux(j,1) = res(j-nx);
    end
    Faux = simplify(Faux);
    Fun(1:nx) = Forig;
    Fun(nx*0*i+1+nx:nx*1*i+nx) = Faux(nx*N*0+nx+1 : nx*N*0+nx+1+i*nx-1);
    Fun(nx*1*i+1+nx:nx*2*i+nx) = Faux(nx*N*1+nx+1 : nx*N*1+nx+1+i*nx-1);
    Fun(nx*2*i+1+nx:nx*3*i+nx) = Faux(nx*N*2+nx+1 : nx*N*2+nx+1+i*nx-1);
    vars = [allxvars(          1 : nx*N*0+nx+1+i*nx-1)...
            allxvars(nx*N*1+nx+1 : nx*N*1+nx+1+i*nx-1)...
            allxvars(nx*N*2+nx+1 : nx*N*2+nx+1+i*nx-1)];
    augfuncell{i} = odeFunction(Fun,vars,u1,u2,u3,epsilon);
    toc 
end
statefuncell = odeFunction(Forig(1:7),allxvars(1:7),u1,u2,u3);
odefunvec.augfun = augfuncell;
odefunvec.statefun = statefuncell; 
save('functionhandle.mat','odefunvec')