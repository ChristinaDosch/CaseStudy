% This script contains the entire simulation using SO, conducting the following steps: 
% 1.) Initialization of common parameters and constraints (independent of ansatz)
% 2.) Initialization of parameters and obj. function corresponding to SO
% 3.) Performing optimization using SO

% Important Notes:
% this implements the closed-form formula (see "Summary" in "Ansatz 1:
% Stochastic Optimization")
%% 1.) Initialize parameters
cost = 1;
epsilon = 0.05;
penalty = @(x) x*2; % linear penalty
N = 1; % one hour schedule
t = linspace(0,24,N);

% Constraints
[x_min, x_max, delta, A, b] = constraints(N);

%% 2.) Initialize optimization model for SO
H1=@(y) normcdf((1+epsilon)*y,5,30); %normal distribution function, mu and sigma might be added
H2=@(y) normcdf((1-epsilon)*y,5,30);
objfct = @(x) obj_SO_closed_form(x,H1,H2,cost,penalty,epsilon); %objective function using the closed-form version of SO

%% 3.) Performing optimization using SO
% [x_opt, obj_opt] = ga(objfct,N,A,b,[],[],0,x_max); % genetic algorithm
tic
[x_opt, obj_opt] = patternsearch(objfct,x0,A,b,[],[],x_min*ones(1,N),x_max*ones(1,N)); % pattern search
toc
