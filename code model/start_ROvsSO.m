function [x_opt_RO, obj_opt_RO, obj_RO,...
          x_opt_SO, obj_opt_SO, obj_SO] = start_ROvsSO(ToPlotOrNotToPlot)

[x_opt_RO, obj_opt_RO, ~] = start_RO(false,'smooth');
display('RO smooth finished')
[x_opt_SO, obj_opt_SO, ~] = start_SO_closed(false);
display('SO closed finished')

N = 10000;
[T, P, cost, penalty, ~, epsilon, C, SOC_0, t, mu, sigma, lambda, ~] = init_parameters;
[x_min, x_max, delta, ~, ~, ~, ~, ~, ~, ~, ~] = init_constraints(T, P, C, SOC_0);
variance = sigma.^2;
E = mvnrnd(mu,variance,N);

cost = ones(N,1)*cost;
X = ones(N,1)*x_opt_RO;
    obj_RO = sum(-cost.*E + penalty(max(zeros(N,T), (X - P*epsilon) - E)) + max(zeros(N,T), E - (X + P*epsilon)).*cost, 2);
X = ones(N,1)*x_opt_SO;
    obj_SO = sum(-cost.*E + penalty(max(zeros(N,T), (X - P*epsilon) - E)) + max(zeros(N,T), E - (X + P*epsilon)).*cost, 2);
if nargin == 0, ToPlotOrNotToPlot = true; end
if ToPlotOrNotToPlot
    figure, hold on,
    plot(t,x_opt_RO,'*r', t,x_opt_SO,'*b',... % solution computed by ga or patternsearch
         t,x_opt_RO + epsilon*P,'^r', t,x_opt_SO + epsilon*P,'^b',...
         t,x_opt_RO - epsilon*P,'vr', t,x_opt_SO - epsilon*P,'vb',...
         [t(1) t(end)], [x_max x_max], 'k--',... % x_max
         t, mu - lambda*sigma, 'k+-.', t, mu + lambda*sigma, 'k+-.',... % uncertainty intervals
         [t(1) t(end)], [x_min x_min], 'k--',... % x_min
         t,mu,'ko') % centers of uncertainty intervals
    xlabel('time'), ylabel('energy, kWh')
    % size
    xlim([0 24])
    v = axis;
    ylim([0 v(4)]);
    % info-box at the top left corner
    text(0.05*v(2),0.93*v(4),['cost = ', num2str(cost(1,:))])
    text(0.05*v(2),0.89*v(4),['penalty = ', func2str(penalty)])
    text(0.05*v(2),0.85*v(4),['[x_{min} x_{max}] = ', '[', num2str(x_min), ' ', num2str(x_max), ']'])
    text(0.05*v(2),0.81*v(4),['\Delta = ', num2str(delta)])
    title('RO(red) vs SO(blue)')
    hold off
end