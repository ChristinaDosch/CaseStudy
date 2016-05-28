function [T, P, cost, penalty, epsilon, t, mu, sigma, lambda] = init_parameters

T = 1 + 24*1;        % one hour schedule
P = 3.8;
epsilon = 0.02;
cost = 1;
penalty = @(x) x*2;  % linear penalty
t = linspace(0,24,T);

pd = makedist('Normal');
mu = P*(0.2 + pdf(pd,(t-12)/sqrt(12)));
sigma = 0.2*ones(1,T);
lambda = 1.96;
