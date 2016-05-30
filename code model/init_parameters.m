function [T, P, cost, penalty, epsilon, t, mu, sigma, lambda] = init_parameters

%T = 1 + 24*1;        % one hour schedule
T = 1440;            % every-minute schedule
P = 3.8;             % nominal power of the PV element
epsilon = 0.05;
cost = 1;
penalty = @(x) x.^2;  % quadratic penalty
t = linspace(0,24,T);

pd = makedist('Normal');
mu = P*(0.2 + pdf(pd,(t-12)/sqrt(12)));
sigma = 0.2*ones(1,T);
lambda = 1.96;
