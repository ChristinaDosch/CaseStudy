function [x_opt, obj_opt, runningTime, obj_true] = start_SO_closed(ToPlotOrNotToPlot)
% This function contains the entire simulation using SO, conducting the following steps: 
% 1.) Initialization of common parameters and constraints (independent of ansatz)
% 2.) Initialization of parameters and obj. function corresponding to SO
% 3.) Performing optimization using SO
% 4.) Plotting solutions and data
%
% Important Note:
% this implements the closed-form formula (see "Summary" in "Ansatz 1:
% Stochastic Optimization") for general T and multidimensional distribution

%% 1.) Initialize parameters and constraints
display('1.) Initialize parameters and constraints')
if nargin == 0, ToPlotOrNotToPlot = true; end
[T, P, cost, penalty, ~, epsilon, C, SOC_0, t, mu, sigma] = init_parameters;
[x_min, x_max, delta, ~, ~, A, b] = init_constraints(T,P,C,SOC_0);

%% 2.) Initialize optimization model for SO
display('2.) Initialize optimization model for SO')
% multivariate normal distribution for E
H = cell(T); % cell array of function handles, H(i) contains function handle on cdf for E_i
for i=1:T
    H{i} = @(y) normcdf(y,mu(i),sigma(i));
end

objfct = @(x) obj_SO_closed_form(x,H,mu,cost,penalty,epsilon,P); % objective function using the closed-form version of SO, also containing the gradient in the second argument

%% 3.) Performing optimization
display('3.) Performing optimization')
tic
options = optimoptions('fmincon','GradObj','on');
[x_opt, obj_opt] = fmincon(objfct,mu,A,b,[],[],x_min*ones(1,T),x_max*ones(1,T),[],options);
runningTime = toc;

%% 4.) Computing true objective
display('4.) Computing true objective')
N = 10000;
Cost = ones(N,1) * cost;
E = mvnrnd(mu,sigma,N); E = max(E,zeros(size(E)));

X = ones(N,1)*x_opt;
    obj_true = sum(-Cost.*E + penalty(max(zeros(N,T), (X - P*epsilon) - E)) + max(zeros(N,T), E - (X + P*epsilon)).*Cost, 2);
    
%% 5.) Plot the solutions and data
if ToPlotOrNotToPlot
    display('5.) Plot the solutions and data')
    % Active ramping constraints (up to 0.001)
    arc = abs(abs(x_opt(2:end) - x_opt(1:end-1)) - delta) < 0.001;
    arc_= [false arc];
    arc = [arc false];
    
    figure, hold on
    plot(t,x_opt,'*r',... % solution computed by ga or patternsearch         
         t,mu,'ko',... % expected values
         [t(1) t(end)], [x_max x_max], 'k--',... % x_max
         [t(1) t(end)], [x_min x_min], 'k--',...% x_min
         [t(arc); t(arc_)], [x_opt(arc); x_opt(arc_)],'r') % active ramping constraints
          
    legend('calculated opt. sol.',...
           'expected solar radiation',...
           'x_{max}, x_{min}')
    xlabel('time'), ylabel('energy, kWh')
    % size
    xlim([0 24])
    v = axis;
    ylim([0 v(4)]);
    % info-box at the top left corner
    text(0.05*v(2),0.98*v(4),['true mean obj. value = ', num2str(mean(obj_true))])
    text(0.05*v(2),0.94*v(4),['true max obj. value = ', num2str(max(obj_true))])
    text(0.05*v(2),0.90*v(4),['optimal value = ', num2str(obj_opt)])
    text(0.05*v(2),0.86*v(4),['cost = ', num2str(cost)])
    text(0.05*v(2),0.82*v(4),['penalty = ', func2str(penalty)])
    text(0.05*v(2),0.78*v(4),['[x_{min} x_{max}] = ', '[', num2str(x_min), ' ', num2str(x_max), ']'])
    text(0.05*v(2),0.74*v(4),['\Delta = ', num2str(delta)])
    title('SO, closed form')
    hold off
end
end