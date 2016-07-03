function [x_opt, obj_opt, runningTime, obj_true] = start_RO(ToPlotOrNotToPlot, SmoothOrNonSmooth)
% This function contains the entire simulation using RO, conducting the following steps: 
% 1.) Initialization of common parameters and constraints (independent of ansatz)
% 2.) Initialization of parameters and obj. function corresponding to RO
% 3.) Performing optimization using RO
% 4.) Plotting solutions and data
%
% It follows the documentation in "Ansatz 2: Robust Optimization"

%% 1.) Initialize parameters and constraints
display('1.) Initialize parameters and constraints')
switch nargin
    case 0, ToPlotOrNotToPlot = true; SmoothOrNonSmooth = 'nonsmooth';
    case 1, SmoothOrNonSmooth = 'nonsmooth';
end

[T, P, cost, penalty, penalty_grad, epsilon, C, SOC_0, t, mu, sigma, lambda, ~] = init_parameters_inout;
[x_min, x_max, delta, ~, ~, A, b] = init_constraints_inout(T, P, C, SOC_0);
A = A(:,1:T);
%% 2.) Setting up optimization problem for RO
display('2.) Setting up optimization problem for RO')
e_l = max(mu - lambda*sigma,0); %
e_u = mu + lambda*sigma; %

% Initialize optimization model for RO
x0 = mu; % Starting guess for pattern search
switch SmoothOrNonSmooth
    case 'nonsmooth',
        objfct = @(x) obj_RO(x,e_l,e_u,cost,penalty,penalty_grad,epsilon,P,'nonsmooth'); % Nonsmooth objective function for RO
    case 'smooth'
        objfct = @(x) obj_RO(x,e_l,e_u,cost,penalty,penalty_grad,epsilon,P,'smooth'); % Smooth objective function for RO
    otherwise, error('SmoothOrNonSmooth must be either "smooth" or "nonsmooth"')
end
%% 3.) Performing optimization using RO
display('3.) Performing optimization using RO')
switch SmoothOrNonSmooth
    case 'nonsmooth'
        options = psoptimset('MaxFunEvals', 5000*T);
        tic
        [x_opt, obj_opt] = patternsearch(objfct,...
            x0, A, b, [], [], x_min*ones(1,T), x_max*ones(1,T), options); % pattern search
        runningTime = toc;
    case 'smooth'
        options = optimoptions(@fmincon,'MaxFunEvals', 5000*T);
        tic
        [x_opt, obj_opt] = fmincon(objfct,...
            x0, A, b, [], [], x_min*ones(1,T), x_max*ones(1,T), [], options); % ???
        runningTime = toc;
end

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
    
    % Plot
    figure, hold on,
    plot(t,x_opt,'*b',... % solution computed by interior point (fmincon) or patternsearch
         t,mu,'k',... % centers of uncertainty intervals
         [t(1) t(end)], [x_max x_max], 'k--',... % x_max
         t, e_l, 'k-.',... % uncertainty intervals
         [t(1) t(end)], [x_min x_min], 'k--',... % x_min
         [t(arc); t(arc_)], [x_opt(arc); x_opt(arc_)],'b') % active ramping constraints
    legend('calculated opt. sol.',...
           'expected solar radiation',...
           'x_{max}, x_{min}',...
           'lower bound of uncertainty intervals')
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
    title(['RO, ' SmoothOrNonSmooth])
    hold off
end