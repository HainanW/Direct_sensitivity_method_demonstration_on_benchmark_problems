%% 
clear all 
%% Part 1.Options for ODE & NLP Solvers
odetol = 1e-8;
NLPtol = 1e-6;
optODE = odeset('RelTol', odetol, 'AbsTol', odetol);
optNLP = optimset('GradObj','on', 'Display', 'iter', ...
    'TolX', NLPtol, 'TolFun', NLPtol, 'TolCon', NLPtol,...
    'MaxFunEvals',10000,'Algorithm','SQP');
%% Part 2.Model Pars and Algorithm Pars  
m.N = 3; % # of time stages 
m.nx = 2;
m.nu = 2;
m.nx_t = m.nx + m.nx*(m.nu*m.N); % 2+ 2*(2*3) = 14  
m.e1  = 18000;
m.e2  = 30000;
m.k10 = 0.535*10^11;
m.k20 = 0.461*10^18;
m.beta1 = 0.53;
m.beta2 = 0.43;
m.alpha = m.e2/m.e1;
m.c = m.k20/(m.k10)^m.alpha;

x0 = [m.beta1 m.beta2 zeros(1,m.nx_t-m.nx)]';
lb = [0.1*ones(1,m.N)  zeros(1,m.N)]; 
ub = [2.0*ones(1,m.N) 8*ones(1,m.N)];
%% Part 3.Run optimization with fmincon 
% dvar0 = [0.166771842272956	0.144829433223544	0.129227256151918 ...
%          1.63129998575985	2.49250008431194	3.87619992992821]; %JO_scaled = -0.67936690; 
dvar0 = [1*ones(1,m.N) 8/3*ones(1,m.N)];

load('functionHandle.mat')
[J0, grad0] = costfun(dvar0,m,odefunvec,optODE,x0);

problem = createOptimProblem('fmincon','objective',...
    @(dvar)costfun(dvar,m,odefunvec,optODE,x0),...
    'x0',dvar0,'options',optNLP,'lb',lb,'ub',ub,...
    'nonlcon',  @(dvar)NLPconst(dvar,m));
tstart = tic; 
[dvarO,JO,exitflag,output,lambda,grad,hessian]  = fmincon(problem);
tend = toc(tstart);
%% Part 4.This is the states + Sensitivity Information from 3 Approaches 

load('ResVecO_Approach1.mat','ResVecO') % Approach1 ! Sensitivity wrt time
load('zvec_scaled.mat', 'zvec_scaled')% Approach3 with CPET
dvar0 = [0.166771842272956	0.144829433223544	0.129227256151918 ...
         1.63129998575985	2.49250008431194	3.87619992992821]; %JO_scaled = -0.67936690; 

[J0, grad0, zvec0] = costfun_withSensitivityOut(dvar0,m,odefunvec,optODE,x0);
ts_scaled = zvec_scaled(:,1); 
topt_vec = zeros(1,m.N); 
for i = 1 : m.N
    topt_vec(i) = sum(dvar0(4:4+i-1));
end 
ts_convert = zeros(length(ts_scaled),1);

% Approach 3: Convert the time scaled time to original time 
for i = 1 : length(ts_scaled)
    if ts_scaled(i) <= 1 
       ts_convert(i) = dvar0(4)*ts_scaled(i);
    elseif ts_scaled(i)  > 1 && ts_scaled(i) <= 2
       ts_convert(i) = dvar0(4) + dvar0(5)*(ts_scaled(i)-1);
    else 
       ts_convert(i) = sum(dvar0(4:5)) + dvar0(6)*(ts_scaled(i)-2);
    end 
end 

% Approach 1: Convert 
t1 = ResVecO(:,1);

s11_end = ResVecO(end,10) + ResVecO(end,12) + ResVecO(end,14);
s12_end = ResVecO(end,12) + ResVecO(end,14);
s13_end = ResVecO(end,14);

s21_end = ResVecO(end,11) + ResVecO(end,13) + ResVecO(end,15);
s22_end = ResVecO(end,13) + ResVecO(end,15);
s23_end = ResVecO(end,15);
% p1: Approach1 t as decision variables 
% p2: Approach2 delt as decision variables
% p3: Approach3 CPET 
clf 
figure(1)
p1 = plot(t1(end), s11_end,'k^', t1(end), s21_end,'k^',...
          t1(end), s12_end,'k^', t1(end), s22_end,'k^',...
          t1(end), s13_end,'k^', t1(end), s23_end,'k^');

p1(1).MarkerSize = 8;
p1(2).MarkerSize = 8;
p1(3).MarkerSize = 8;
p1(4).MarkerSize = 8;
p1(5).MarkerSize = 8;
p1(6).MarkerSize = 8;
hold on 
p2 = plot(zvec0(1:2:end,1), zvec0(1:2:end,10:15), 'k-');
p3 = plot(ts_convert(1:2:end,1), zvec_scaled(1:2:end,10:15),'ko','Markersize',3);
hold off 
legend({'Approach 1 at \itt_{\itF}','','','','','',...
        'Approach 2','','','','','',...
        'Approach 3'})
xlim([0 9])
xticks(0:1:8)
xticklabels(0:1:8)
xlabel('Time (min)')
ylabel('Sensitivity w.r.t. \it{{\delta}_k} ')

print -dsvg Figure2-R.svg
% End 

