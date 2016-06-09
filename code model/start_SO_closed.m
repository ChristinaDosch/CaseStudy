function [x_opt, obj_opt, runningTime] = start_SO_closed(ToPlotOrNotToPlot)
% This function contains the entire simulation using SO, conducting the following steps: 
% 1.) Initialization of common parameters and constraints (independent of ansatz)
% 2.) Initialization of parameters and obj. function corresponding to SO
% 3.) Performing optimization using SO
% 4.) Plotting solutions and data
%
% Important Note:
% this implements the closed-form formula (see "Summary" in "Ansatz 1:
% Stochastic Optimization") for general T and multivariate distribution

%% 1.) Initialize parameters and 
if nargin == 0, ToPlotOrNotToPlot = true; end
[T, P, cost, penalty, epsilon, C, SOC_0, t, mu, sigma] = init_parameters;
[x_min, x_max, delta, A, b] = init_constraints(T,P,C,SOC_0);

%% 2.) Initialize optimization model for SO
% multivariate normal distribution for E
H = cell(T); % cell array of function handles, H(i) contains function handle on cdf for E_i
for i=1:T
    H{i} = @(y) normcdf(y,mu(i),sigma(i));
end

objfct = @(x) obj_SO_closed_form(x,H,mu,cost,penalty,epsilon,P); % objective function using the closed-form version of SO

%% 3.) Performing optimization using SO
tic
[x_opt, obj_opt] = fmincon(objfct,mu,A,b,[],[],x_min*ones(1,T),x_max*ones(1,T));
runningTime = toc;

%% 4.) Plot the solutions and data
if ToPlotOrNotToPlot
    figure, hold on
    plot(t,x_opt,'*r',... % solution computed by ga or patternsearch
         t,x_opt + epsilon*P,'^r',...
         t,x_opt - epsilon*P,'vr',...
         [t(1) t(end)], [x_max x_max], 'k--',... % x_max
         t, mu - 1.960*sigma, 'k+-.', t, mu + 1.960*sigma, 'k+-.',... % uncertainty intervals
         [t(1) t(end)], [x_min x_min], 'k--',... % x_min
         t,mu,'ko') % centers of uncertainty intervals
    legend('calculated opt. sol.',...
           'upper no-penalty bound',...
           'lower no-penalty bound',...
           'x_{max}, x_{min}',...
           'uncertainty intervals (95%)')
    xlabel('time'), ylabel('energy, kWh')
    % size
    xlim([0 24])
    v = axis;
    ylim([0 v(4)]);
    % info-box at the top left corner
    text(0.05*v(2),0.98*v(4),['optimal value = ', num2str(obj_opt)])
    text(0.05*v(2),0.93*v(4),['cost = ', num2str(cost)])
    text(0.05*v(2),0.90*v(4),['penalty = ', func2str(penalty)])
    text(0.05*v(2),0.87*v(4),['[x_{min} x_{max}] = ', '[', num2str(x_min), ' ', num2str(x_max), ']'])
    text(0.05*v(2),0.84*v(4),['\Delta = ', num2str(delta)])
    title('SO, closed form')
    hold off
end
end