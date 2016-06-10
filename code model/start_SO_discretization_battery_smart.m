function [x_opt, obj_opt, runningTime] = start_SO_discretization_battery_smart(ToPlotOrNotToPlot)
%% DISCRETIZATION APPROACH with BATTERY SMART APPROACH
% This script solves our well known optimization problem using the
% discretization approach and including
%       * ramping constraints and x_min, x_max
%       * battery using the SMART approach
%       * general penalty function (also quadratic is possible)

%% Initialize parameters
if nargin == 0, ToPlotOrNotToPlot = true; end
[T, P, cost, penalty, epsilon, C, SOC_0, t] = init_parameters;

%% Constraints
[x_min, x_max, SOC_min, SOC_max, delta,~,~,~,~,A_smart,b_smart] = init_constraints(T,P,C,SOC_0);

%% Example scenarios

% using SAMPLE_NORMAL_INDEPENDENT.CSV:
% for this example T = 96 is required in init_parameters!!! (15min interv)
E = load('sample_normal_independent.csv');
E = 1/1000 * E; % since we need kWh (and in the samples it's in Wh)
%K = size(E,1); % number of realizations
K = 300;

% using PVDATA2.MAT
% for this example T = 1440 is required in init_parameters!!! (1min interv)
%E = load('PVdata2');
%E = 1/1000 * E;
%K = 372; % number of realizations
%E = reshape(PVdata2,372,1440); % array of K realizations (one per row) with data per minute of one day, respectively
%E = reshape(PVdata2(:,1),31,1440); % array of 31 realizations with minute values from January

F = cell(K,1);

% determine revenue function F(x,\tilde{x}^k) for every k=1,...,K:
for k = 1:K 
x_tilde = @(x) E(k,:)+0.95*x((3*T+1):(4*T))-x((2*T+1):(3*T)); % compute \tilde{x}^k
F(k) = { @(x) obj_SO_discr(x(1:T),x_tilde(x),cost,penalty,epsilon,P)};
end

objfct = @(x) 1/K * sum(cellfun(@(f)f(x),F)); % weighted (all weights=1/K) sum of F(x,e^k)

%% Performing optimization
x0 = zeros(1,4*T);
tic
%[x_opt, obj_opt] = patternsearch(objfct,x0,A_smart,b_smart,[],[],...
 %   [x_min*ones(1,T), SOC_min*ones(1,T),0*ones(1,(2*T))],[x_max*ones(1,T), SOC_max*ones(1,T),2*P*ones(1,(2*T))]);
[x_opt, obj_opt] = fmincon(objfct, x0, A_smart, b_smart,[],[],...
    [x_min*ones(1,T), SOC_min*ones(1,T),0*ones(1,(2*T))],[x_max*ones(1,T), SOC_max*ones(1,T),2*P*ones(1,(2*T))]);
runningTime = toc
%% 4.) Plot the solutions and data
if ToPlotOrNotToPlot
    figure, hold on,
    plot(t,x_opt(1:T),'*r',... % solution computed by ga or patternsearch
         t,x_opt(1:T) + epsilon*P,'^r',...
         t,x_opt(1:T) - epsilon*P,'vr',...
         t,x_opt(T+1:2*T)*C, 'bo',...
         t,x_opt(2*T+1:3*T), '-b',...
         t,x_opt(3*T+1:4*T), '-g',...
         [t(1) t(end)], [x_max x_max], 'k--',... % x_max
         [t(1) t(end)], [x_min x_min], 'k--') % x_min 
    legend('calculated opt. sol.',...
           'upper no-penalty bound',...
           'lower no-penalty bound',...
           'SOC',...
           'b in',...
           'b out',...
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
    text(0.05*v(2),0.84*v(4),['\Delta = ', num2str(delta)])
    title('SO discretization brute force')
    hold off
end
end