function [xb_opt, obj_opt, runningTime, obj_true] = start_RO_B_inout(ToPlotOrNotToPlot, SmoothOrNonSmooth, option)
% This script contains the entire simulation using RO WITH BATTERY, conducting the following steps: 
% 1.) Initialization of common parameters and constraints (independent of ansatz)
% 2.) Initialization of parameters and obj. function corresponding to RO
% 3.) Performing optimization using RO
% 4.) Plotting solutions and data
%
% It follows the documentation in "Ansatz 2: Robust Optimization"

%% 1.) Initialize parameters and constraints
display('1.) Initialize parameters and constraints')
switch nargin
    case 0, ToPlotOrNotToPlot = true; SmoothOrNonSmooth = 'nonsmooth'; option = 2;
    case 1, SmoothOrNonSmooth = 'nonsmooth'; option = 2;
    case 2, option = 2;
end
[T, P, cost, penalty, penalty_grad, epsilon, C, SOC_0, t, mu, sigma, lambda] = init_parameters_inout;
[x_min, x_max, delta, SOC_min, SOC_max, A, b, A_, b_, A__, b__] = init_constraints_inout(T, P, C, SOC_0);
A = [A; A_; A__];
b = [b; b_; b__];
%% 2.) Setting up optimization problem for RO
display('2.) Setting up optimization problem for RO')
e_l = max(mu - lambda*sigma,0); %
e_u = mu + lambda*sigma; %

% Initialize 
b_in0 = min(random('unif',0,0.1*C,1,T), e_l);
b_out0 = min(random('unif',0,0.1*C,1,T), e_l);
xb0 = [zeros(1,T) b_in0 b_out0]; % Starting guess for pattern search
k = 0;

while any(A*xb0' > b) && (k<1000) 
    b_in0 = min(random('unif',0,0.1*C,1,T), e_l);
    b_out0 = min(random('unif',0,0.1*C,1,T), e_l);
    xb0 = [zeros(1,T) b_in0 b_out0];
    k = k+1;
end
if k == 1000, error('failed to find an admissible starting point'), end
display([num2str(k), ' attempts were needed to find a feaseble starting point']);
objfct = @(xb) obj_RO_B_inout(xb,e_l,e_u,cost,penalty,penalty_grad,epsilon,P,SmoothOrNonSmooth,option); % Objective function for RO

%% 3.) Performing optimization using RO
display('3.) Performing optimization using RO')
switch SmoothOrNonSmooth
    case 'nonsmooth',
        options = psoptimset('MaxFunEvals', 5000*T, 'TolMesh', 1e-8);
        tic
        [xb_opt, obj_opt] = patternsearch(objfct,xb0,...
            A, b, [], [], [x_min*ones(1,T) zeros(1,2*T)], [x_max*ones(1,T) e_l C*ones(1,T)], options); % pattern search
        runningTime = toc;
    case 'smooth',
        options = optimoptions(@fmincon,'MaxFunEvals', 5000*T,'MaxIter',2000);
        tic
        [xb_opt, obj_opt] = fmincon(objfct,xb0,...
            A, b, [], [], [x_min*ones(1,T) zeros(1,2*T)], [x_max*ones(1,T) e_l C*ones(1,T)], [], options); % ???
        runningTime = toc;
end
x_opt = xb_opt(1:T);
b_in_opt = xb_opt(T+1:2*T);
b_out_opt = xb_opt(2*T+1:3*T);

% Complementrarity check for b_in_opt and b_out_opt
cl = sum(b_in_opt.*b_out_opt);
if cl < 1e-6,
    display('Complementrarity check passed up to a 1e-6 tolerance')
else
    display('Complementrarity check for b_in_opt and b_out_opt failed')
    display(['Complementrarity lack = ', num2str(cl)])
end

%% 4.) Computing true objective
display('4.) Computing true objective')
N = 10000;
Cost = ones(N,1) * cost;
B_in = ones(N,1) * b_in_opt;
B_out = ones(N,1) * b_out_opt;

E = mvnrnd(mu,sigma,N); E = max(E,zeros(size(E)));
E = E + 0.95*B_out - B_in;

X = ones(N,1)*x_opt;
    obj_true = sum(-Cost.*E + penalty(max(zeros(N,T), (X - P*epsilon) - E)) + max(zeros(N,T), E - (X + P*epsilon)).*Cost, 2);

%% 5.) Plot the solutions and data
display('5.) Plot the solutions and data')
if ToPlotOrNotToPlot
    % Active ramping constraints (up to 0.001)
    arc = abs(abs(x_opt(2:end) - x_opt(1:end-1)) - delta) < 0.001;
    arc_= [false arc];
    arc = [arc false];
    
    figure, hold on,
    plot(t,x_opt,'*b',... % solution computed by fmincon or patternsearch
         t,mu,'k') % centers of uncertainty intervals (expected radiation)
%     plot(t,x_opt + epsilon*P,'^r',...
%          t,x_opt - epsilon*P,'vr')
    plot([t(1) t(end)], [x_max x_max], 'k--',... % x_max
         t, e_l, 'k-.') % uncertainty intervals
    plot(t, cumsum(0.95*b_in_opt - b_out_opt) + SOC_0*C, 'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'b') % battary load
    plot([t; t], [zeros(1,T); b_in_opt], '-b',... % battery usage b_in_opt (charge)
         [t; t], [zeros(1,T); -b_out_opt], '-b',... % battery usage b_out_opt (discharge)
         [t(1) t(end)], [x_min x_min], 'k--',... % x_min
         [t(1) t(end)], [C*SOC_max C*SOC_max], 'b--',... % 95% of capacity
         [t(1) t(end)], [C*SOC_min C*SOC_min], 'b--') % 10% of capacity
    plot(t, e_l - (b_in_opt - 0.95*b_out_opt), 'b+-.',... % new uncertainty intervals
         t, e_u - (b_in_opt - 0.95*b_out_opt), 'b+-.') % new uncertainty intervals
    plot([t(arc); t(arc_)], [x_opt(arc); x_opt(arc_)],'b') % active ramping constraints
        
    legend('calculated opt. sol.',...
           'expected solar radiation',...
           'x_{max}, x_{min}',...
           'lower bound of uncertainty intervals',...
           'state of charge',...
           'battery usage')
    xlabel('time'), ylabel('energy, kWh')
    % size
    xlim([0 24])
    v = axis;
%     ylim([-1 v(4)]);
    % info-box at the top left corner
    text(0.05*v(2),0.98*v(4),['true mean obj. value = ', num2str(-mean(obj_true))])
    text(0.05*v(2),0.94*v(4),['true worst obj. value = ', num2str(-max(obj_true))])
    text(0.05*v(2),0.90*v(4),['optimal obj. value = ', num2str(-obj_opt)])
%     text(0.05*v(2),0.86*v(4),['cost = ', num2str(cost)])
%     text(0.05*v(2),0.82*v(4),['penalty = ', func2str(penalty)])
%     text(0.05*v(2),0.78*v(4),['[x_{min} x_{max}] = ', '[', num2str(x_min), ' ', num2str(x_max), ']'])
    text(0.05*v(2),0.85*v(4),['\Delta = ', num2str(delta)])
%     text(0.05*v(2),0.70*v(4),['C = ', num2str(C)])
    title(['RO inout, battery, ', SmoothOrNonSmooth, ', option ', num2str(option)])
    hold off
end