clear all
clc 
%% Part 1.Options for ODE & NLP Solvers
odeTol = 1e-10;
NLPTol = 1e-8;
optODE = odeset('RelTol', odeTol, 'AbsTol', odeTol);
optNLP = optimset('GradObj','on','GradConstr','off',...
    'DerivativeCheck','off', 'Display', 'iter', ...
    'TolX', NLPTol,'TolFun', NLPTol, 'TolCon', NLPTol, ...
    'MaxFunEvals',3000,'MaxIter',1000,'Algorithm','interior-point');
%% Part 2.Model Pars and Algorithm Pars
m.N = 10; % # of time stages
m.nx = 9;
m.nu = 3;
m.nx_t = m.nx + m.nx*(m.nu*m.N); % 9+ 9*(3*10) = 279
m.tf = 9; % 
m.ts = linspace(0,m.N,m.N+1);
x0_orig = [0 22 0 0 -1 0 0];
x0_const = [x0_orig 0 0];
x0_red  = [x0_const zeros(1,m.nx*m.nu)]; 

lb = [-2.83374*ones(1,m.N) -0.80865*ones(1,m.N)  zeros(1,m.N)];
ub = [+2.83374*ones(1,m.N) +0.71265*ones(1,m.N) 9*ones(1,m.N)];
sigma = 1e-6;
rho = 1;
q = 2;
epsilon = 2*sigma/((1+sqrt(5))*rho*q*m.N);
%% Part 3.load sensitivity function handle 
% Built sens. func in build_sensitivity.m
% Now load it from functionhandle.m
load('functionhandle.mat');
load('OptRes032223_Round6.mat','points')

dvarO_Pyomo = [2.833596064491985, -0.0038133355140898676, 0.4174487666110814,0.6999549590642804, 0.9186759922513857, 1.073815774457458, 0.9302904224573103, ...
0.5775724258679925, -0.11150914836543367, 2.833728593932601, ...
0.05230522738102885, -0.005638956767162429, 0.0006594018879565948, 0.003373992470865873, 0.11085863354045986, -0.07241185587227153, 0.00043001732048645014, 0.7126324959937229, ...
0.7126440735216086, 0.7126418427879623, ... 
0.13090470552752673, 0.8578752237482009, 1.0782375580771617, 1.2597902963899896, 1.3633850365023688, 2.1795797057182953, 0.7278324107552094, 0.6545771050678583, 0.520398859396996, ...
0.2274190988163944]; %J0 = 0.01030695
dvarO_Matlab = [2.83341871110044	-0.00517151251504766	0.413207853501381	0.696202442880169	0.917050698436233	1.07364670681192	0.931352022434512	0.577862536225025	-0.110586018843131	2.83371421527009	5.97243592117872e-05	-5.71420813365664e-06	1.57830340561305e-06	4.91402279998671e-06	0.122206022649508	-0.0767819834875334	-2.22906159192196e-06	0.712610376104849	0.712636626397327	0.712631573221118	0.131023465916379	0.850544636678883	1.06818936948555	1.25993745367952	1.37253091898834	2.18463005015648	0.729879428371767	0.654779132616060	0.521181125042945	0.227304419064018];


[JO_Py,gradO_Py]= costfunTS(dvarO_Pyomo, m, odefunvec, optODE, x0_red, epsilon, rho);
[JO_Matlab,gradO_Matlab]= costfunTS(dvarO_Matlab, m, odefunvec, optODE, x0_red, epsilon, rho);
dvar0 = points(4,:); 
%% Part 5. 
dvarOvec = []; 
[JO,gradO]= costfunTS(dvar0, m, odefunvec, optODE, x0_red, epsilon, rho);
JOvec = [];
counter = 1;
JOvec(counter) = JO; 
tendvec = [];
while counter <=10 
    %%
    tstart = tic;
    problem = createOptimProblem('fmincon','objective',...
        @(dvar)costfunTS(dvar, m, odefunvec, optODE, x0_red, epsilon, rho),...
        'x0',dvar0,'options',optNLP,'lb',lb,'ub',ub,'nonlcon',...
        @(dvar)consfunTS(dvar, m, odefunvec, optODE, x0_orig));
    [dvarO,JO,exitflag,output,lambda,gradOpt,hessian]= fmincon(problem);
    tend = toc(tstart);
    tendvec = [tendvec tend];
    fprintf('tf= %9.8f\n',JO)
    if abs(JO - JOvec(counter)) > sigma
        epsilon = epsilon*0.1; 
        rho = rho*10;
        JOvec(counter+1) = JO;
        dvar0 = dvarO;
        dvarOvec = [dvarOvec; dvarO];
    else 
        JOvec(counter+1) = JO;
        dvarOvec = [dvarOvec; dvarO];
        break 
    end 
  
    %%
    counter = counter +1; 
end
% [JOTS,gradTS]= costfun1(dvarO,nu,nx,nx_t,N,odefunvec,optODE,rho,ts,x0_red,epsilon);
% [JOFD,gradFD]= costfunFD(dvarO,N,optODE,rho,x0_FD,epsilon);

ShowMulti_Control(dvarO,m)
% save('OptRes012923.mat')
return
%% Multiple Start
% ms = MultiStart('StartPointsToRun','bounds');
% stpoints = RandomStartPointSet('NumStartPoints',10, ...
%     'ArtificialBound',10);
% [dvarO,JO,flag,outpt,allmins] = run(ms,problem,stpoints);
%% Examine Result
trec = [];
trec_FD = [];
zrec = [];
zrec_FD = [];
z0 = [0.02 0]';
z0_FD = [0.02 0 0];
u1 = dvarO(1:N);
u2 = dvarO(N+1:2*N);
for ks = 1 : N
    [t,z] = ode15s(@(t,x)dyneqn1(t,x,u1,u2,ks,N),...
        [ts(ks),ts(ks+1)],z0,optODE);
    trec = [trec;t];
    zrec = [zrec;z];
    z0 = z(end,:)';
end

ts_FD = zeros(1,N+1);
for i = 1: N
    ts_FD(i+1) = sum(dvarO(N+1:N+i));
end

for ks = 1: N
    [t_FD,z_FD] = ode15s(@(t,x)dyneqnFD(t,x,u1,ks,epsilon),...
        [ts_FD(ks),ts_FD(ks+1)],z0_FD,optODE);
    trec_FD = [trec_FD;t_FD];
    zrec_FD = [zrec_FD;z_FD];
    z0_FD = z_FD(end,:)';
end

%% Visualize results
figure(1)
yyaxis left 
a = axis; 
plot(trec_FD,zrec_FD(:,1),[0 a(2)],[0.026 0.026],'b--')
ylabel('x_1')
yyaxis right 
plot(trec_FD,zrec_FD(:,2))
ylabel('x_2')

figure(2)

showResult_var(dvarO(1:N),dvarO(N+1:2*N),N);
%% Get the value of x3(tau=1) for each iteration 
zvec_Rec = [];
sigma = 2e-4;
rho = 1;
q = 1;
tf = 1;
epsilon_Ini = 2*sigma/((1+sqrt(5))*rho*q*tf);
for i = 1: 6
    dvarin = dvarOvec(i,:);
    epsilon = epsilon_Ini/(10^(i-1));
    zvec = myfunint1(dvarin,N,odefunvec,optODE,ts,x0_sen,epsilon);
    zvec_Rec = [zvec_Rec zvec];
end 
