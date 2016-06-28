function [x_opt1, b_opt1, obj_opt1, obj1, x_opt2, b_opt2, obj_opt2, obj2] = start_RO_comparison
%% Parameter initialization
[T, P, cost, penalty, ~, epsilon, C, SOC_0, t, mu, sigma, lambda, ~] = init_parameters;
[x_min, x_max, delta, SOC_min, SOC_max, ~, ~, ~, ~, ~, ~] = init_constraints(T, P, C, SOC_0);

e_l = max(mu - lambda*sigma, 0);
e_u = mu + lambda*sigma;
%% Solving optimization plroblems
[xb_opt1, obj_opt1, runningTime1] = start_RO(false,'nonsmooth');
display(['first run finished in ', num2str(runningTime1), ' sec'])
x_opt1 = xb_opt1(1:T);
b_opt1 = xb_opt1(T+1:end);

[xb_opt2, obj_opt2, runningTime2] = start_RO_B(false,'smooth',2);
display(['second run finished in ', num2str(runningTime2), ' sec'])
x_opt2 = xb_opt2(1:T);
b_opt2 = xb_opt2(T+1:end);
%% Computing true objective
N = 10000;
cost = ones(N,1)*cost;
E = mvnrnd(mu,sigma,N);
X = ones(N,1)*x_opt1;
    obj1 = sum(-cost.*E + penalty(max(zeros(N,T), (X - P*epsilon) - E)) + max(zeros(N,T), E - (X + P*epsilon)).*cost, 2);
X = ones(N,1)*x_opt2;
    obj2 = sum(-cost.*E + penalty(max(zeros(N,T), (X - P*epsilon) - E)) + max(zeros(N,T), E - (X + P*epsilon)).*cost, 2);
%% Active ramping constraints (up to 0.001)
arc1 = abs(abs(x_opt1(2:end) - x_opt1(1:end-1)) - delta) < 0.001;
arc_1= [false arc1];
arc1 = [arc1 false];

arc2 = abs(abs(x_opt2(2:end) - x_opt2(1:end-1)) - delta) < 0.001;
arc_2= [false arc2];
arc2 = [arc2 false];
%% Plot
figure, hold on,
    plot(t,x_opt1,'*r', t,x_opt2,'*b',... % solution computed by fmincon or patternsearch
         t,x_opt1 + epsilon*P,'^r', t,x_opt2 + epsilon*P,'^b',...
         t,x_opt1 - epsilon*P,'vr', t,x_opt2 - epsilon*P,'vb',...
         [t(1) t(end)], [x_max x_max], 'k--',... % x_max
         t, e_l, 'k+-.', t, e_u, 'k+-.') % uncertainty intervals
    plot(t, cumsum(b_opt2) + SOC_0*C, 'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'b') % battery load
    plot([t; t], [zeros(1,T); b_opt2], '-b',... % battery usage
         [t(1) t(end)], [x_min x_min], 'k--',... % x_min
         t,mu,'ko',... % centers of uncertainty intervals
         [t(1) t(end)], [C*SOC_max C*SOC_max], 'b--',... % max load
         [t(1) t(end)], [C*SOC_min C*SOC_min], 'b--',... % min load
         t, e_l - b_opt2, 'b+-.',... % new uncertainty intervals
         t, e_u - b_opt2, 'b+-.',... % new uncertainty intervals
         [t(arc1); t(arc_1)], [x_opt1(arc1); x_opt1(arc_1)],'r', [t(arc2); t(arc_2)], [x_opt2(arc2); x_opt2(arc_2)],'b')
