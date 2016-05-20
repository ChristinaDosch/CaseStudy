% This script contains the entire simulation using SO, conducting the following steps: 
% 1.) Initialization of common parameters and constraints (independent of ansatz)
% 2.) Initialization of parameters and obj. function corresponding to SO
% 3.) Performing optimization using SO

% Important Notes:
% this implements the closed-form formula (see "Summary" in "Ansatz 1:
% Stochastic Optimization")
%% 1.) Initialize parameters
T = 3; % one hour schedule

cost = 1*ones(1,T); % cost vector
epsilon = 0.05;
penalty = @(x) 2*ones(1,T)*x; % linear penalty function
P=3.8; %nominal power of PV-element

t = linspace(0,24,T);

% Constraints
[x_min, x_max, delta, A, b] = constraints(T);

%% 2.) Initialize optimization model for SO

mu = 5; sigma = 30; % mu and sigma for normal distribution, still random
H = @(y) normcdf(y,mu,sigma); %distribution function for E, here: normal distribution function
exp = 5*ones(T,1); %expected value E[E]

objfct = @(x) obj_SO_closed_form(x,H,exp,cost,penalty,epsilon,P); %objective function using the closed-form version of SO

%% 3.) Performing optimization using SO
% [x_opt, obj_opt] = ga(objfct,N,A,b,[],[],0,x_max); % genetic algorithm
x0 = zeros(T,1);
tic
[x_opt, obj_opt] = patternsearch(objfct,x0,A,b,[],[],x_min*ones(1,T),x_max*ones(1,T)); % pattern search
toc
