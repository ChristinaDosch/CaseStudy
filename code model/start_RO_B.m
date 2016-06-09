function [xb_opt, obj_opt, runningTime] = start_RO_B(ToPlotOrNotToPlot)
% This script contains the entire simulation using RO WITH BATTERY, conducting the following steps: 
% 1.) Initialization of common parameters and constraints (independent of ansatz)
% 2.) Initialization of parameters and obj. function corresponding to RO
% 3.) Performing optimization using RO
% 4.) Plotting solutions and data
%
%It follows the documentation in "Ansatz 2: Robust Optimization"

%% 1.) Initialize parameters and constraints
if nargin == 0, ToPlotOrNotToPlot = true; end
[T, P, cost, penalty, epsilon, C, t, mu, sigma, lambda] = init_parameters;
[x_min, x_max, delta, A, b, A_, b_] = init_constraints(T,P,C);

%% 2.) Initialize e_l, e_r for RO
e_l = mu - lambda*sigma; %
e_u = mu + lambda*sigma; %

% Initialize optimization model for RO
xb0 = [mu zeros(1,T)]; % Starting guess for pattern search
objfct = @(xb) obj_RO_B(xb,e_l,e_u,cost,penalty,epsilon,P); % Objective function for RO

%% 3.) Performing optimization using RO
A = [A zeros(size(A,1),size(A_,2)); zeros(size(A_,1),size(A,2)) A_];
tic
[xb_opt, obj_opt] = patternsearch(objfct,xb0,...
    A,[b; b_],[],[],[x_min*ones(1,T) -100*ones(1,T)],[x_max*ones(1,T) e_l]); % pattern search
runningTime = toc;
x_opt = xb_opt(1:T);
b_opt = xb_opt(T+1:end);
%% 4.) Plot the solutions and data
if ToPlotOrNotToPlot
    figure, hold on,
    plot(t,x_opt,'*r',... % solution computed by ga or patternsearch
         t,x_opt + epsilon*P,'^r',...
         t,x_opt - epsilon*P,'vr',...
         [t(1) t(end)], [x_max x_max], 'k--',... % x_max
         t, e_l, 'k+-.',... % uncertainty intervals
         t, cumsum(b_opt), 'bo',... % battary load
         [t; t], [zeros(1,T); b_opt], '-b',... % battery usage
         [t(1) t(end)], [x_min x_min], 'k--',... % x_min
         t,mu,'ko',... % centers of uncertainty intervals
         t, e_u, 'k+-.',... % uncertainty intervals
         [t(1) t(end)], [C C], 'b--',... % capacity
         t, e_l - b_opt, 'b+-.',... % new uncertainty intervals
         t, e_u - b_opt, 'b+-.') % new uncertainty intervals
        
    legend('calculated opt. sol.',...
           'upper no-penalty bound',...
           'lower no-penalty bound',...
           'x_{max}, x_{min}',...
           'uncertainty intervals',...
           'state of charge',...
           'battery usage')
    xlabel('time'), ylabel('energy, kWh')
    % size
    xlim([0 24])
    v = axis;
    ylim([-1 v(4)]);
    % info-box at the top left corner
    text(0.05*v(2),0.98*v(4),['optimal value = ', num2str(obj_opt)])
    text(0.05*v(2),0.93*v(4),['cost = ', num2str(cost)])
    text(0.05*v(2),0.90*v(4),['penalty = ', func2str(penalty)])
    text(0.05*v(2),0.87*v(4),['[x_{min} x_{max}] = ', '[', num2str(x_min), ' ', num2str(x_max), ']'])
    text(0.05*v(2),0.84*v(4),['\Delta = ', num2str(delta)])
    text(0.05*v(2),0.81*v(4),['C = ', num2str(C)])
    title('RO')
    hold off
end