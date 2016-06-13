function [T, P, cost, penalty, penalty_grad, epsilon, C, SOC_0, t, mu, sigma, lambda] = init_parameters

%T = 1 + 24*1;       % one hour schedule (works for RO and start_SO_closed, for the latter only with long running time)
%T = 3;              % for start_SO_closed with an acceptable runtime
%T = 1440;           % every-minute schedule (required for start_SO_discretization and start_SO_discretization_battery when using PVdata2)
T = 96;              % 15min-schedule (required for start_SO_discretization and start_SO_discretization_battery when using sample_normal_independent)
         
P = 3.8;             % nominal power of the PV element
epsilon = 0.05;
cost = ones(1,T)*1;   % cost(j) is a price for an energy unit during j's hour
penalty = @(x) x.^2;  % quadratic penalty
penalty_grad = @(x) 2*x; % derivative of the penalty function
%penalty = @(x) 1*x; % linear penalty
%penalty_grad = @(x) 1;
C = 2.6;              % battery capacity
SOC_0 = 0.25;         % state of charge (SOC) at day break
t = linspace(0,24,T);

pd = makedist('Normal');
%mu = P*(0.2 + pdf(pd,(t-12)/sqrt(12)));
mu = P*(pdf(pd,(t-12)/sqrt(12)));
sigma = 0.2*ones(1,T);
lambda = 1.96;        % lambda*sigma is the width of the confidence intervals in RO
