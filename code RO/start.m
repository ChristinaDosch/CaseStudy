% 1.) Initialization of parameters, objective function and constraints.
% 2.) Performing robust optimization.
% 3.) Plotting solutions and data
%% Initialize parameters
cost = 1;
epsilon = 0.1;
penalty = @(x) x*2; % linear penalty
N = 24*1; % one hour schedule
t = linspace(0,24,N);
%% Initialize e_l, e_r
pd = makedist('Normal');
e = 0.2 + pdf(pd,(t-12)/sqrt(12)); % gaussian radiance distribution during the day
mu = 0.2; % unsertanty interval width parameter
e_l = e*(1-mu);
e_u = e*(1+mu);
%% Initialize optimization model
x0 = e; % Starting guess for pattern search
objfct = @(x) big_objective(x,e_l,e_u,cost,penalty,epsilon); % Objective function

% constraints, see "First approach with RO" for explanation
x_max = 0.55; % tba
delta = 0.025; % tba
B = [-eye(N-1) zeros(N-1,1)] + [zeros(N-1,1) eye(N-1)];
A = [B; -B];
b = ones(2*(N-1),1)*delta;
%% Performing optimization
% [x_opt, obj_opt] = ga(objfct,N,A,b,[],[],0,x_max); % genetic algorithm
tic
[x_opt, obj_opt] = patternsearch(objfct,x0,A,b,[],[],zeros(1,N),x_max*ones(1,N)); % pattern search
toc
%% Plot the solutions and data
x1 = (penalty(e_u) + cost*e_l)/...                                          % Is an intersection point of two lines:
                    (cost*(1-epsilon) + penalty(1+epsilon));                % penalty*((1-epsilon)x-e_l) and cost*(e_u-(1+epsilon)x)
figure, hold on
plot(t,x1,'ob',... % solution for the uncoupled (unconstrained) problem
     t,x_opt,'*r',... % solution computed by ga or patternsearch
     [t(1) t(end)], [x_max x_max], 'k--',... % x_max
     [t; t], [e_l; e_u], 'k+-.',... % uncertainty intervals
     t,e,'ko') % centers of uncertainty intervals
legend('opt. sol. for unconstrained problem',...
       'calculated opt. sol.',...
       'x_{max}',...
       'uncertainty intervals')
xlabel('time'), ylabel('energy, kWh'), hold off