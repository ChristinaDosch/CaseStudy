% This script contains the entire simulation using SO, conducting the following steps: 
% 1.) Initialization of common parameters and constraints (independent of ansatz)
% 2.) Initialization of parameters and obj. function corresponding to SO
% 3.) Performing optimization using SO
% 4.) Plotting solutions and data

% Important Note:
% this implements the closed-form formula (see "Summary" in "Ansatz 1:
% Stochastic Optimization") for general T and multivariate distribution
%% 1.) Initialize parameters
T = 3; % one hour schedule

cost = 2 * ones(1,T); % cost vector
epsilon = 0.05;
penalty = @(x) 2*ones(1,T) * x; % linear penalty function, output: scalar
P=3.8; % nominal power of PV-element

t = linspace(0,24,T);

% Constraints
[x_min, x_max, delta, A, b] = constraints(T);

%% 2.) Initialize optimization model for SO

% multivariate normal distribution for E
mu = 5*ones(T,1); sigma = 30*ones(T,1); % mu and sigma for normal distributions, still random
H = cell(T); % cell array of function handles, H(i) contains function handle on cdf for E_i
for i=1:T
    H{i} = @(y) normcdf(y,mu(i),sigma(i));
end    

exp = mu; % vector of expected values of E

% all E_i independently distributed
%mu = 5; sigma = 30; % mu and sigma for normal distribution, still random
%H = @(y) normcdf(y,mu,sigma); %distribution function for E, here: normal distribution function
%exp = 5*ones(T,1); %expected value E[E]

objfct = @(x) obj_SO_closed_form(x,H,exp,cost,penalty,epsilon,P); %objective function using the closed-form version of SO

%% 3.) Performing optimization using SO
% [x_opt, obj_opt] = ga(objfct,T,A,b,[],[],0,x_max); % genetic algorithm
x0 = zeros(T,1);
tic
[x_opt, obj_opt] = patternsearch(objfct,x0,A,b,[],[],x_min*ones(1,T),x_max*ones(1,T)); % pattern search
toc

%% 4.) Plot the solutions and data
figure, hold on
plot(t,x_opt,'*r',... % solution computed by ga or patternsearch
     t,x_opt+epsilon*P,'^r',...
     t,x_opt-epsilon*P,'vr',...
     [t(1) t(end)], [x_max x_max], 'k--',... % x_max
     [t(1) t(end)], [x_min x_min], 'k--')    % x_min
legend('calculated opt. sol.',...
       'upper no-penalty bound',...
       'lower no-penalty bound',...
       'x_{max}, x_{min}')
xlabel('time'), ylabel('energy, kWh')

% size
xlim([0 24])
v = axis;
ylim([0 v(4)]);

% info-box at the top left corner
text(0.05*v(2),0.98*v(4),['optimal value = ', num2str(obj_opt)])
%text(0.05*v(2),0.93*v(4),['cost = ', num2str(cost)])
%text(0.05*v(2),0.90*v(4),['penalty = ', func2str(penalty)])
%text(0.05*v(2),0.87*v(4),['[x_{min} x_{max}] = ', '[', num2str(x_min), ' ', num2str(x_max), ']'])
%text(0.05*v(2),0.84*v(4),['\Delta = ', num2str(delta)])
hold off
