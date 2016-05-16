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
P=3.8; %nominal power of PV-element
T = 1; % one hour schedule
t = linspace(0,24,T);

% Constraints
[x_min, x_max, delta, A, b] = constraints(T);

%% 2.) Initialize optimization model for SO

H=@(y) normcdf(y,5,30); %distribution function for E, here: normal distribution function, mu and sigma are still random
exp=5;                  %expected value E[E]
objfct = @(x) obj_SO_closed_form(x,H,exp,cost,penalty,epsilon,P); %objective function using the closed-form version of SO

%% 3.) Performing optimization using SO
% [x_opt, obj_opt] = ga(objfct,N,A,b,[],[],0,x_max); % genetic algorithm
x0=0;
tic
[x_opt, obj_opt] = patternsearch(objfct,x0,A,b,[],[],x_min*ones(1,T),x_max*ones(1,T)); % pattern search
toc
