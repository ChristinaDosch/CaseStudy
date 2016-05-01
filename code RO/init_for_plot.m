% Initialization of parameters.
% Was used for preparation before plotting with plot_objective.

% penalty functions
penalty = @(x) x*2; % linear penalty
% penalty = @(x) 0.1*x.^2; % quadratic penalty

cost = 1;
epsilon = 0.1;
e = 50;
mu = 0.05;
e_l = e*(1-mu);
e_u = e*(1+mu);

x = 20:0.1:80;
