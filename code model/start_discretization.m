%% DISCRETIZATION APPROACH
% This script solves our well known optimization problem using the
% discretization approach and including
%       * ramping constraints and x_min, x_max
%       * no battery constraints
%       * general penalty function (also quadratic is possible)

%% Initialize parameters
cost = 1;
epsilon = 0.05;
penalty = @(x) x.^2; % quadratic penalty
P = 3.8; % nominal power
T = 1440; % every minute schedule for a day

%% Constraints
x_min = 0;  
x_max = 0.7*P; 
delta = 0.03*P; % maximum deviation allowed
B = [-eye(T-1) zeros(T-1,1)] + [zeros(T-1,1) eye(T-1)];
A = [B; -B];
b = ones(2*(T-1),1)*delta;

%% Example scenarios

K = 372; % number of realizations
E = reshape(PVdata2,372,1440); % array of K realizations (one per row) with data per minute of one day, respectively
F = cell(K,1);
for k = 1:K % determine revenue function F(x,e^k) for every e^1,...,e^K
e = E(k,:);
F(k) = { @(x) revenue(x,e,cost,penalty,epsilon,P)};
end

objfct = @(x) 1/K * sum(cellfun(@(f)f(x),F)); % weighted (all weights=1/K) sum of F(x,e^k)

%% Performing optimization
tic
[x_opt, obj_opt] = fmincon(objfct,x0,A,b,[],[],x_min*ones(1,T),x_max*ones(1,T)); % pattern search
toc
