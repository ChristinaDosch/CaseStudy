function [x_opt, obj_opt, runningTime] = start_SO_discr_bat_smart_constr(ToPlotOrNotToPlot)
%% DISCRETIZATION APPROACH with BATTERY SMART APPROACH
% This script solves our well known optimization problem using the
% discretization approach and including
%       * ramping constraints and x_min, x_max
%       * battery using the SMART approach (using O(T*K) constraints)
%       * general penalty function (also quadratic is possible)

%% Initialize parameters
if nargin == 0, ToPlotOrNotToPlot = true; end
[T, P, cost, penalty, penalty_grad, epsilon, C, SOC_0, t, mu, ~, ~, penalty_hess] = init_parameters;

%% Example scenarios

% using SAMPLE_NORMAL_INDEPENDENT.CSV:
% for this example T = 96 is required in init_parameters!!! (15min intervalls)
E = load('sample_normal_independent.csv');
E = 1/1000 * E; % since we need kWh (and in the samples it's in Wh)
K = 1; % number of realizations to use
%E(1,:) = mu;

%% Initialize Constraints
[x_min, x_max, delta, SOC_min, SOC_max,~,~,~,~, A_smart, b_smart] = init_constraints(T,P,C,SOC_0,K,E(1:K,:));

%% Determine objective function
% x \in R^(T+3*K*T) is going to be the optimization variable in this function. 
% x = [x,SOC^1,...,SOC^K,b^{in,1},...b^{in,K},b^{out,1},...b^{out,K}]
objfct = @(x) obj_SO_discr_weighted_sum(x,E(1:K,:),K,cost,penalty,penalty_grad,penalty_hess,epsilon,P); % contains gradient as second argument
%hess = @(x, lambda) hess_lagr_SO_discr_smart( transpose(x), E,K,cost,penalty,penalty_grad,penalty_hess,epsilon,P );
%% Performing optimization
% create starting guess x0
x0 = zeros(1,T+3*K*T);
x0(T+1:2*T) = SOC_0;
for i=10:floor(T/2)
   x0(i) = min(x0(i-1)+delta,x_max); 
end 
for i=floor(T/2)+1:T-10
   x0(i) = max(0,x0(i-1)-delta);
end
options = optimoptions('fmincon','Algorithm','sqp','SpecifyObjectiveGradient',true,'Diagnostics','on');%,...
    %  'StepTolerance',1e-1000,'MaxFunEvals', 3000, 'MaxIterations', 1000);
%options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'Hessian','user-supplied','HessFcn',@hessianfcn,'MaxFunEvals',30000,'MaxIter',10000);%,'MaxFunEvals', 30000);

% possible options:
% * 'ScaleProblem','obj-and-constr': causes the algorithm to normalize all constraints and the objective function
%                                    didn't help at all


% start the solver
tic
[x_opt, obj_opt] = fmincon(objfct, x0, A_smart(1:2*(T-1)+K*T,:), b_smart(1:2*(T-1)+K*T),... % inequality constraints
    A_smart(2*(T-1)+K*T+1:2*(T-1)+2*(T*K),:),b_smart(2*(T-1)+K*T+1:2*(T-1)+2*(T*K)),...     % equality constraints
    [x_min*ones(1,T), SOC_min*ones(1,K*T),0*ones(1,(2*K*T))],...                    % lower bounds
    [x_max*ones(1,T), SOC_max*ones(1,K*T),2*P*ones(1,(2*K*T))],[],options);         % upper bounds
runningTime = toc
%% Plot the solutions and data
if ToPlotOrNotToPlot
    figure, hold on,
    plot(t,x_opt(1:T),'*r',... % solution computed by fmincon
         t,x_opt(1:T) + epsilon*P,'^r',...
         t,x_opt(1:T) - epsilon*P,'vr',...
         t,x_opt(T+1:2*T)*C, 'bo',... % \tilde{SOC^1}
         t,x_opt(T+K*T+1:T+K*T+T), '-b',... % \tilde{b^in,1}
         t,x_opt(T+2*K*T+1:T+2*K*T+T), '-g',... % \tilde{b^out,1}
         t,E(1,:)+0.95*x_opt((T+2*K*T+1):(T+2*K*T+T))-x_opt((T+K*T+1):(T+K*T+T)),'*b',...% \tilde{x^1}
         [t(1) t(end)], [x_max x_max], 'k--',... % x_max
         [t(1) t(end)], [x_min x_min], 'k--') % x_min 
    legend('calculated opt. sol.',...
           'upper no-penalty bound',...
           'lower no-penalty bound',...
           'SOC,1',...
           'b in,1',...
           'b out,1',...
           'x tilde 1',...
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
    title('SO discretization smart')
    hold off
end
end
