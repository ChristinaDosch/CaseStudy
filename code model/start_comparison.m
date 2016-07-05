function [obj1, obj2] = start_comparison(...
    x1, x2,...
    ToPlotOrNotToPlot, my_title,...
    b_in1, b_out1, b_in2, b_out2)
%% 1~red, 2~blue

%% Parameter initialization
[T, P, cost, penalty, ~, epsilon, C, SOC_0, t, mu, sigma, lambda, ~] = init_parameters_inout;
[x_min, x_max, delta, SOC_min, SOC_max, ~, ~, ~, ~, ~, ~] = init_constraints_inout(T, P, C, SOC_0);

switch nargin
    case 2,...
        ToPlotOrNotToPlot = false; my_title = 'comparison';
        b_in1 = zeros(1,T); b_out1 = zeros(1,T); b_in2 = zeros(1,T); b_out2 = zeros(1,T);
        battery1 = false; battery2 = false;
    case 3,
        my_title = 'comparison';
        b_in1 = zeros(1,T); b_out1 = zeros(1,T); b_in2 = zeros(1,T); b_out2 = zeros(1,T);
        battery1 = false; battery2 = false;
    case 4,
        b_in1 = zeros(1,T); b_out1 = zeros(1,T); b_in2 = zeros(1,T); b_out2 = zeros(1,T);
        battery1 = false; battery2 = false;
    case 5,
        error('b_in_opt1 provided but b_out_opt1 is missing'),
    case 6,
        b_in2 = zeros(1,T); b_out2 = zeros(1,T);
        battery1 = true; battery2 = false;
    case 7,
        error('b_in_opt2 provided but b_out_opt2 is missing'),
    case 8,
        battery1 = true; battery2 = true;
end

e_l = max(mu - lambda*sigma, 0);
% e_u = mu + lambda*sigma;

%% Solving optimization plroblems
% [xb_opt1, obj_opt1, runningTime1] = start_RO(false,'nonsmooth');
% display(['first run finished in ', num2str(runningTime1), ' sec'])
% x_opt1 = xb_opt1(1:T);
% b_opt1 = xb_opt1(T+1:end);
% 
% [xb_opt2, obj_opt2, runningTime2] = start_RO_B(false,'nonsmooth',1);
% display(['second run finished in ', num2str(runningTime2), ' sec'])
% x_opt2 = xb_opt2(1:T);
% b_opt2 = xb_opt2(T+1:end);
%% Computing true objective
N = 10000;
Cost = ones(N,1) * cost;
B_in1 = ones(N,1) * b_in1;
B_out1 = ones(N,1) * b_out1;
B_in2 = ones(N,1) * b_in2;
B_out2 = ones(N,1) * b_out2;

E = mvnrnd(mu,sigma,N); E = max(E,zeros(size(E)));
E1 = E + 0.95*B_out1 - B_in1;
E2 = E + 0.95*B_out2 - B_in2;

X = ones(N,1)*x1;
    obj1 = sum(-Cost.*E1 + penalty(max(zeros(N,T), (X - P*epsilon) - E1)) + max(zeros(N,T), E1 - (X + P*epsilon)).*Cost, 2);
X = ones(N,1)*x2;
    obj2 = sum(-Cost.*E2 + penalty(max(zeros(N,T), (X - P*epsilon) - E2)) + max(zeros(N,T), E2 - (X + P*epsilon)).*Cost, 2);
%% Active ramping constraints (up to 0.001)
arc1 = abs(abs(x1(2:end) - x1(1:end-1)) - delta) < 0.001;
arc_1= [false arc1];
arc1 = [arc1 false];

arc2 = abs(abs(x2(2:end) - x2(1:end-1)) - delta) < 0.001;
arc_2= [false arc2];
arc2 = [arc2 false];
%% Plot
if ToPlotOrNotToPlot,
    % main plot with the optimal solutions
    figure, hold on
        plot(t, x2, '*r', t, x1, '*b',... % provided solutions
             [t(1) t(end)], [x_max x_max], 'k--',... % x_max
             [t(1) t(end)], [x_min x_min], 'k--',... % x_min
             t, e_l, 'k-.') % uncertainty intervals
        if battery1,
            plot(t, cumsum(0.95*b_in1 - b_out1) + SOC_0*C, 'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'b')
            plot([t; t], [zeros(1,T); b_in1], '-b',... % battery usage b_in1 (charge)
                 [t; t], [zeros(1,T); -b_out1], '-b') % battery usage b_out1 (discharge)) % battery load
        end
        if battery2,
            plot(t, cumsum(0.95*b_in2 - b_out2) + SOC_0*C, 'bo', 'MarkerSize', 4, 'MarkerFaceColor', 'r',...
                 [t; t], [zeros(1,T); b_in2], '-r',... % battery usage b_in2 (charge)
                 [t; t], [zeros(1,T); -b_out2], '-r') % battery usage b_out2 (discharge)) % battery load
        end
        if battery1 || battery2,
            plot([t(1) t(end)], [C*SOC_max C*SOC_max], 'b--',... % max load
                 [t(1) t(end)], [C*SOC_min C*SOC_min], 'b--') % min load
        end
        plot(t,mu,'k',... % centers of uncertainty intervals
             [t(arc1); t(arc_1)], [x1(arc1); x1(arc_1)],'b', [t(arc2); t(arc_2)], [x2(arc2); x2(arc_2)],'r')
        plot(t,cost/4.5,'g')
    % info-box at the top left corner
    xlim([0 24])
    v = axis;
    text(0.05*v(2),0.98*v(4),['true mean obj. value = ', num2str(-mean(obj1)), ' (blue), ', num2str(-mean(obj2)), ' (red)'])
    text(0.05*v(2),0.94*v(4),['true worst obj. value = ', num2str(-max(obj1)), ' (blue), ', num2str(-max(obj2)), ' (red)'])
%     if all(cost == cost(1)), cost_label = ['constant ', num2str(cost(1))];
%     else cost_label = 'variable'; end
%     text(0.05*v(2),0.85*v(4),['cost ', cost_label])
%     text(0.05*v(2),0.81*v(4),['penalty = ', func2str(penalty)])
%     text(0.05*v(2),0.77*v(4),['[x_{min} x_{max}] = ', '[', num2str(x_min), ' ', num2str(x_max), ']'])
    text(0.05*v(2),0.90*v(4),['\Delta = ', num2str(delta)])
%     text(0.05*v(2),0.71*v(4),['C = ', num2str(C)])
    title(my_title)
    hold off
%     % plot with realizations
%     figure, hold on
%     % plot with the best and worst realizations
%         i_good1  = find(obj1 == min(obj1), 1, 'first');
%         i_bad1 = find(obj1 == max(obj1), 1, 'first');
%         i_good2  = find(obj2 == min(obj2), 1, 'first');
%         i_bad2 = find(obj2 == max(obj2), 1, 'last');
%         plot(t, E(i_good1,:), 'r', t, E(i_bad1,:), 'r--',...
%              t, E(i_good2,:), 'b', t, E(i_bad2,:), 'b--',...
%              t, x1, '*r', t, x2, '*b')
%         plot(t, E(1:40:end,:), 'ok', 'MarkerSize', 1)
%         plot(t, e_l, 'k+-.', t, e_u, 'k+-.')
%         legend('best realizaion for the red solution',...
%                'worst realizaion for the red solution',...
%                'best realizaion for the blue solution',...
%                'worst realizaion for the blue solution',...
%                'first solution',...
%                'second solution',...
%                'realization of E')
%         title('Worst and best realizations for both solutions')
%     hold off
end