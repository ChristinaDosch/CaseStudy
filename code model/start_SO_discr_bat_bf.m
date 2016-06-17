function [x_opt, obj_opt, runningTime] = start_SO_discr_bat_bf(ToPlotOrNotToPlot)
%% DISCRETIZATION APPROACH with BATTERY BRUTE FORCE
% This script solves our well known optimization problem using the
% discretization approach and including
%       * ramping constraints and x_min, x_max
%       * battery using the BRUTE FORCE approach
%       * general penalty function (also quadratic is possible)

%% Initialize parameters
if nargin == 0, ToPlotOrNotToPlot = true; end
[T, P, cost, penalty, penalty_grad, epsilon, C, SOC_0, t] = init_parameters;

%% Constraints
[x_min,x_max,delta, ~, ~,A,b] = init_constraints(T,P,C,SOC_0);

%% Example scenarios

% using SAMPLE_NORMAL_INDEPENDENT.CSV:
% for this example T = 96 is required in init_parameters!!!
E = load('sample_normal_independent.csv');
E = 1/1000 * E; % since we need kWh (and in the samples it's in Wh)
%K = size(E,1); % number of realizations
K = 1;

% using PVDATA2.MAT
% for this example T = 1440 is required in init_parameters!!! 
%E = load('PVdata2');
%E = 1/1000 * E;
%K = 372; % number of realizations
%E = reshape(PVdata2,372,1440); % array of K realizations (one per row) with data per minute of one day, respectively
%E = reshape(PVdata2(:,1),31,1440); % array of 31 realizations with minute values from January

F = cell(K,1);

% determine revenue function F(x,\tilde{x}^k) for every k=1,...,K:
for k = 1:K 
x_tilde = @(x) battery(E(k,:),x,C,SOC_0); % compute \tilde{x}^k
F(k) = { @(x) obj_SO_discr(x,x_tilde(x),cost,penalty,penalty_grad,epsilon,P)};
end

objfct = @(x) 1/K * sum(cellfun(@(f)f(x),F)); % weighted (all weights=1/K) sum of F(x,e^k)

%% Performing optimization
x0 = zeros(1,T);
tic
[x_opt, obj_opt] = patternsearch(objfct,x0,A,b,[],[],x_min*ones(1,T),x_max*ones(1,T));
%[x_opt, obj_opt] = fmincon(objfct,x0,A,b,[],[],x_min*ones(1,T),x_max*ones(1,T));
runningTime = toc
%% 4.) Plot the solutions and data
if ToPlotOrNotToPlot
    figure, hold on,
    plot(t,x_opt,'*r',... % solution computed by ga or patternsearch
         t,x_opt + epsilon*P,'^r',...
         t,x_opt - epsilon*P,'vr',...
         [t(1) t(end)], [x_max x_max], 'k--',... % x_max
         [t(1) t(end)], [x_min x_min], 'k--') % x_min 
    legend('calculated opt. sol.',...
           'upper no-penalty bound',...
           'lower no-penalty bound',...
           'x_{max}, x_{min}')
    xlabel('time'), ylabel('energy, kWh')
    % size
    xlim([0 24]);
    v = axis;
    ylim([0 v(4)]);
    % info-box at the top left corner
    text(0.05*v(2),0.98*v(4),['optimal value = ', num2str(obj_opt)])
    text(0.05*v(2),0.93*v(4),['cost = ', num2str(cost)])
    text(0.05*v(2),0.90*v(4),['penalty = ', func2str(penalty)])
    text(0.05*v(2),0.87*v(4),['[x_{min} x_{max}] = ', '[', num2str(x_min), ' ', num2str(x_max), ']'])
%    text(0.05*v(2),0.84*v(4),['\Delta = ', num2str(delta)])
    title('SO discretization brute force')
    hold off
end
end