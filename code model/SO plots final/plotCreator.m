%% Before you run this file:
% * load corresponding results you wanna plot 
% * decide whether you want to include the first sample in your plot (and
%      check if the E that is loaded here is indeed the E you used for your
%      calculations)
% * decide whether you want to make the 10th schedule value to 0 artificially
% * change the penalty function in the info box
% * decide whether you wanna have x_min and x_max displayed in the info box
% * should cost=5 be changed to cost=5x to be consistend with the penalty term?
% * should the optimal obj. values be positive (since we speak of a maximization problem during the talk)

%% load corresponding x_opts and obj_opts
load('T=96_battery_K=3_Penalty7.mat')
x_opt_bat = x_opt;
obj_opt_bat = obj_opt;

load('T=96_noBattery_K=3_Penalty7_STEPSIZE.mat')
x_opt_noBat = x_opt;
obj_opt_noBat = obj_opt;

x_opt_bat(10) = 0;
x_opt_noBat(10) = 0;

%% Initialization
K = 3; 
[T, P, cost, penalty, penalty_grad, epsilon, C, SOC_0, t, mu, sigma, lambda, penalty_hess] = init_parameters;
[x_min, x_max, delta, SOC_min, SOC_max, A, b, A_, b_, A_smart, b_smart] = init_constraints(T, P, C, SOC_0,K);
%E = load('sample_normal_independent.csv');
%E = 1/1000 * max(E, 0); E_good(1,:) = E(2,:); E = E_good;

%% active ramping constraints
% for marking active ramping constraints in battery plan
x_x_opt_bat = x_opt_bat(1:T);
    arc = abs(abs(x_x_opt_bat(2:end) - x_x_opt_bat(1:end-1)) - delta) < 0.001;
    arc_= [false arc];
    arc = [arc false];    
% for marking active ramping constraints in the plan without battery
x_x_opt_noBat = x_opt_noBat(1:T);
    arc = abs(abs(x_x_opt_noBat(2:end) - x_x_opt_noBat(1:end-1)) - delta) < 0.001;
    arc_= [false arc];
    arc = [arc false];  
    
%% plot    
    figure, hold on,
    plot(t,x_opt_bat(1:T),'*b',... % solution computed by sqp with battery
         t,x_opt_noBat(1:T),'*r',... % solution without battery % t,E(1,:),'xk',...% sample 1
         t,mu,'-k',... % (estimated) expected value
         [t(1) t(end)], [x_max x_max], 'k--',... % x_max
         [t; t], [zeros(1,T); x_opt_bat(T+K*T+1:T+K*T+T) - x_opt_bat(T+2*K*T+1:T+2*K*T+T)], '-b',... % b^in - b^out, sodass Strich nach oben, wenn Batterie beladen wird und Strich nach unten, wenn Batterie entladen wird
         [t(1) t(end)], [x_min x_min], 'k--',...) % x_min 
         [t(arc); t(arc_)], [x_x_opt_bat(arc); x_x_opt_bat(arc_)],'b',... % active ramping constraints for battery schedule
         [t(arc); t(arc_)], [x_x_opt_noBat(arc); x_x_opt_noBat(arc_)],'r');
     %t,E(1,:)+0.95*x_opt((T+2*K*T+1):(T+2*K*T+T))-x_opt((T+K*T+1):(T+K*T+T)),'*b',...% \tilde{x^1}
     %t,x_opt(1:T) + epsilon*P,'^r',...
     %t,x_opt(1:T) - epsilon*P,'vr',...    
     %t,x_opt(T+K*T+1:T+K*T+T), '-ob',... % \tilde{b^in,1}         
     %t,x_opt(T+2*K*T+1:T+2*K*T+T), '-og',... % \tilde{b^out,1}        
         legend('optimal sol. with battery',...
           'optimal sol. without battery',...% 'first sample',...
           'expected value',... 
           'x_{max}, x_{min}',...
           'battery usage')
    xlabel('time'), ylabel('energy in kWh')
    %'upper no-penalty bound',...
    %'lower no-penalty bound',...
    
    % size
    xlim([0 24]);
    v = axis;
    ylim([-1 v(4)+0.5]);
    % info-box at the top left corner
    text(0.05*v(2),0.98*(v(4)+0.5),['cost = ', num2str(cost(1))])
    text(0.05*v(2),0.95*(v(4)+0.5),'penalty = 7x^2')
    %text(0.05*v(2),0.87*v(4),['[x_{min}, x_{max}] = ', '[', num2str(x_min), ', ', num2str(x_max), ']'])
    text(0.05*v(2),0.905*(v(4)+0.5),['\Delta = ', num2str(delta)])
    text(0.05*v(2),0.85*(v(4)+0.5),'optimal obj. value with battery = ')
    text(0.243*v(2),0.85*(v(4)+0.5), num2str(abs(obj_opt_bat)),'color','b')
    text(0.05*v(2),0.81*(v(4)+0.5),'optimal obj. value without battery = ')
    text(0.243*v(2),0.81*(v(4)+0.5),num2str(abs(obj_opt_noBat)),'color','r')
    title('SO discretization approach with (blue) and without (red) battery','Fontsize',18)