function [x_opt, obj_opt, runningTime, exitflag, output] = start_SO_discr_bat_smart_constr(ToPlotOrNotToPlot)
%% DISCRETIZATION APPROACH with BATTERY SMART APPROACH
% This script solves our well known optimization problem using the
% discretization approach and including
%       * ramping constraints and x_min, x_max
%       * battery using the SMART approach (using O(T*K) constraints)
%       * general penalty function (also quadratic is possible)

%% Initialize parameters
if nargin == 0, ToPlotOrNotToPlot = true; end
[T, P, cost, penalty, penalty_grad, epsilon, C, SOC_0, t, mu, ~, ~, penalty_hess] = init_parameters;

%% Load example scenarios SAMPLE_NORMAL_INDEPENDENT.CSV:
% for this example T = 96 is required in init_parameters!!! (15min intervalls)

E = load('sample_normal_independent.csv');
%E = load('sample_normal_sum.csv');
E = 1/1000 * max(E, 0); % since we need kWh (and in the samples it's in Wh)
% select "good" samples: 
%old good samples that turned out to be bad: E_good = zeros(11,96); E_good(1:3,:) = E(2:4,:); E_good(4,:) = E(6,:); E_good(5:9,:) = E(8:12,:); E_good(10:11,:) = E(14:15,:);
E_good = zeros(4,96); E_good(1,:) = E(2,:); E_good(2,:) = E(6,:); E_good(3,:) = E(9,:); E_good(4,:) = E(15,:);
E = E_good;

%Martina's good samples: 
% E_good = zeros(11,96); 
% E_good(1,:) = E(2,:);% E(5,:); 
% E_good(2,:) = E(10,:); 
% E_good(3,:) = E(16,:); 
% E_good(4,:) = E(21,:);
% E_good(5,:) = E(28,:);
% E_good(6,:) = E(29,:);% not nec. good
% E_good(7,:) = E(30,:);% not nec. good
% E = E_good;

% E_good for only T = 25 realizations
%E_short = zeros(11,25);
%E_short(1) = 0;
%for i = 1:24
%E_short(:,i+1) = 0.25*(E(:,(i-1)*4+1)+E(:,(i-1)*4+2)+E(:,(i-1)*4+3)+E(:,(i-1)*4+4));
%end
%E = E_short;

K = 1; % number of realizations to use

%% Initialize Constraints
[x_min, x_max, delta, SOC_min, SOC_max,~,~,~,~, A_smart, b_smart] = init_constraints(T,P,C,SOC_0,K);
b2 = reshape(E(1:K,:)',1,K*T); % upper bound for b^in (constructed here such that E does not have to be passed to init_constraints)

%% Determine objective function
% x \in R^(T+3*K*T) is going to be the optimization variable in this function. x = [x,SOC^1,...,SOC^K,b^{in,1},...b^{in,K},b^{out,1},...b^{out,K}]
objfct = @(x) obj_SO_discr_weighted_sum(x,E(1:K,:),K,cost,penalty,penalty_grad,penalty_hess,epsilon,P); % contains gradient as second argument
%hess = @(x, lambda) hess_lagr_SO_discr_smart( transpose(x), E,K,cost,penalty,penalty_grad,penalty_hess,epsilon,P );

%% Performing optimization
% create starting guess x0
x0 = zeros(1,(1+3*K)*T);
x0(T+1:2*T) = SOC_0;
for i=10:floor(T/2)
   x0(i) = min(x0(i-1)+delta,x_max); 
end 
for i=floor(T/2)+1:T-10
   x0(i) = max(0,x0(i-1)-delta);
end


options = optimoptions('fmincon','Algorithm','sqp','GradObj','on','Diagnostics','on');%,...
     % 'StepTolerance',1e-100,'MaxFunEvals', 3000, 'MaxIterations', 3000);

%options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'Hessian','user-supplied','HessFcn',@hessianfcn,'MaxFunEvals',30000,'MaxIter',10000);%,'MaxFunEvals', 30000);

% start the solver
tic
[x_opt, obj_opt, exitflag, output] = fmincon(objfct, x0, A_smart(1:2*(T-1),:), b_smart(1:2*(T-1)),... % inequality constraints 
    A_smart(2*(T-1)+1:2*(T-1)+(T*K),:),b_smart(2*(T-1)+1:2*(T-1)+(T*K)),...         % equality constraints
     [x_min*ones(1,T), SOC_min*ones(1,K*T),0*ones(1,(2*K*T))],...                    % lower bounds
     [x_max*ones(1,T), SOC_max*ones(1,K*T),b2,2*P*ones(1,K*T)],[],options);        % upper bounds
 %    [x_min*ones(1,T), SOC_0*ones(1,K*T),0*ones(1,(2*K*T))],...  % lower bounds without battery
 %    [x_max*ones(1,T), SOC_0*ones(1,K*T),0*ones(1,2*K*T)],[],options); % upper bounds without battery

runningTime = toc
%% Plot the solutions and data
if ToPlotOrNotToPlot
    
    % Active ramping constraints (up to 0.001)
    x_x_opt = x_opt(1:T);
    arc = abs(abs(x_x_opt(2:end) - x_x_opt(1:end-1)) - delta) < 0.001;
    arc_= [false arc];
    arc = [arc false];    
    
    figure, hold on,
    plot(t,x_opt(1:T),'*r',... % solution computed by fmincon
         t, E(1,:), '-k',...%t,mu,'-k',... % (estimated) expected value
         [t(1) t(end)], [x_max x_max], 'k--',... % x_max
         [t; t], [zeros(1,T); x_opt(T+K*T+1:T+K*T+T) - x_opt(T+2*K*T+1:T+2*K*T+T)], '-b',... % b^in - b^out, sodass Strich nach oben, wenn Batterie beladen wird und Strich nach unten, wenn Batterie entladen wird
         t,E(1,:)+0.95*x_opt((T+2*K*T+1):(T+2*K*T+T))-x_opt((T+K*T+1):(T+K*T+T)),'*b',...% \tilde{x^1}
         [t(1) t(end)], [x_min x_min], 'k--',...) % x_min  
         [t(arc); t(arc_)], [x_x_opt(arc); x_x_opt(arc_)],'r', 'markers',8); % active ramping constraints
     %t,x_opt(1:T) + epsilon*P,'^r',...
     %t,x_opt(1:T) - epsilon*P,'vr',...    
     %t,x_opt(T+K*T+1:T+K*T+T), '-ob',... % \tilde{b^in,1}         
     %t,x_opt(T+2*K*T+1:T+2*K*T+T), '-og',... % \tilde{b^out,1}        
     %t,E(1,:)+0.95*x_opt((T+2*K*T+1):(T+2*K*T+T))-x_opt((T+K*T+1):(T+K*T+T)),'*b',...% \tilde{x^1}
    %set(h,{'markers'},{12;5;9})
    fs = 16.5;
    set(gca,'FontSize',fs);
    legend1= legend('calculated opt. sol.',...
           'expected value',... 
           'x_{max}, x_{min}',...
           'battery usage');
    set(legend1,'Fontsize',fs)
    xlabel('time', 'Fontsize',fs), ylabel('energy, kWh', 'Fontsize',fs)
    %'upper no-penalty bound',...
    %'lower no-penalty bound',...
    
    % size
    xlim([0 24]);
    v = axis;
     ylim([-1 v(4)+0.5]);
    % info-box at the top left corner
    %text(0.03*v(2),0.98*(v(4)+0.5),['cost = ', num2str(cost(1))], 'Fontsize',fs)
    %text(0.03*v(2),0.935*(v(4)+0.5),['penalty = x^2'], 'Fontsize',fs)
    %text(0.03*v(2),0.87*(v(4)+0.5),['\Delta = ', num2str(delta)], 'Fontsize',fs)
    text(0.03*v(2),0.91*(v(4)+0.5),['optimal value = ', num2str(-obj_opt)], 'Fontsize',fs)
   % text(0.05*v(2),0.87*v(4),['[x_{min}, x_{max}] = ', '[', num2str(x_min), ', ', num2str(x_max), ']'])
    %title('SO discretization smart','Fontsize',18)
    hold off
end
end
