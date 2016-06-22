function [obj,grad, hessian] = obj_SO_discr_weighted_sum(x,E,K,cost,penalty,penalty_grad,penalty_hess,epsilon,P)
%% OBJ_SO_DISCR_WEIGHTED_SUM computes 1/K * \sum_{i=1}^K F(x,\tilde{x}^k),
% i.e. the objective function of our discretized optimization problem as in
% 1.8 (smart approach)
%
%
% Input:
%       x - 1 by m+3*K*m vector          [x,SOC^1,...,SOC^K,b^{in,1},...,b^{in,K},b^{out,1},...,b^{out,K}] 
%                                        (corresponds to the optimization variable vector in outer opt problem
%                                        start_SO_discr_bat_smart_constr)
%       E - L by m vector                E(i,j) is the solar radiance during the
%                                        j^th hour in the i^th sample, L>=K
%       K - scalar                       sample size, i.e. we only use the
%                                        K first rows of E
%       cost - 1 by m vector             price per energy unit (time
%                                        dependent)
%       penalty - function handle        penalty function
%       penalty_grad function handle     derivative of the penalty function
%       epsilon - pos. scalar << 1       width of no-penalty interval
%       P - scalar                       nominal power of PV element
%
% Output:
%       obj - scalar                     weighted objective function value
%       grad - 1 by m+3*K*m vector       grad w.r.t. [x,SOC^1,...,SOC^K,b^{in,1},...,b^{in,K},b^{out,1},...,b^{out,K}]
%% Get number of time steps T
T = size(E,2);

%% Calculation
F = zeros(K,1); % F(k) = F(x,\tilde{x^k})
G = zeros(K,T+3*K*T); % G(k) = gradient w.r.t. [x,SOC^1,...,SOC^K,b^{in,1},...,b^{in,K},b^{out,1},...,b^{out,K}]
H = zeros(T+3*K*T,T+3*K*T);

for k = 1:K 
x_tilde = E(k,:)+0.95*x((T+2*K*T+(k-1)*T+1):(T+2*K*T+k*T))...
    -x((T+K*T+(k-1)*T+1):(T+K*T+k*T)); % compute \tilde{x}^k as e^k + 0.95 \tilde{b^out,k} - \tilde{b^in,k}
[F(k),grad_SOC_b,grad_x,~] = obj_SO_discr(x(1:T),x_tilde,cost,penalty,penalty_grad,epsilon,P,penalty_hess);

% In the "k^th objective" F(k) only b^{in,k} and b^{out,k} appear and thus
% only the corresponding entries in the gradient row G(k,:) are non-zero:
G(k,1:T) = grad_x; % gradient w.r.t. x
G(k,(T+(k-1)*T+1):(T+k*T)) = grad_SOC_b(1:T); % gradient w.r.t. SOC^k
G(k,((K+1)*T+(k-1)*T+1):((K+1)*T+k*T)) = grad_SOC_b(T+1:2*T); % gradient w.r.t. b^in,k
G(k,((2*K+1)*T+(k-1)*T+1):((2*K+1)*T+k*T)) = grad_SOC_b(2*T+1:3*T); % gradient w.r.t. b^out,k

% WENN MAN HESSE HABEN WILL, MUSS MAN OBEN BEI AUFRUF VON OBJ_SO_DISCR "~"
% DURCH HESSIAN ERSETZEN
%H(1:T,1:T) = H(1:T,1:T)+1/K*hessian(1:T,1:T);
%H(1:T,T+(k-1)*T+1:(T+k*T)) = H(1:T,T+(k-1)*T+1:(T+k*T))+1/K*hessian(1:T,T+1:2*T);
%H(1:T,(K+1)*T+(k-1)*T+1:((K+1)*T+k*T)) =H(1:T,(K+1)*T+(k-1)*T+1:((K+1)*T+k*T)) +1/K* hessian(1:T,2*T+1:3*T);
%H(1:T,(2*K+1)*T+(k-1)*T+1:((2*K+1)*T+k*T)) = H(1:T,(2*K+1)*T+(k-1)*T+1:((2*K+1)*T+k*T))+1/K*...
%    hessian(1:T,3*T+1:4*T);

%H(T+(k-1)*T+1:(T+k*T),1:(T)) = H(T+(k-1)*T+1:(T+k*T),1:(T))+1/K*hessian(T+1:2*T,1:T);
%H(T+(k-1)*T+1:(T+k*T),T+(k-1)*T+1:(T+k*T)) =                    1/K*hessian(T+1:2*T,T+1:2*T);
%H(T+(k-1)*T+1:(T+k*T),(K+1)*T+(k-1)*T+1:((K+1)*T+k*T)) =        1/K*hessian(T+1:2*T,2*T+1:3*T);
%H(T+(k-1)*T+1:(T+k*T),(2*K+1)*T+(k-1)*T+1:((2*K+1)*T+k*T)) =    1/K*hessian(T+1:2*T,3*T+1:4*T);

%H((K+1)*T+(k-1)*T+1:((K+1)*T+k*T),1:(T)) = H((K+1)*T+(k-1)*T+1:((K+1)*T+k*T),1:(T)) +1/K*hessian(2*T+1:3*T,1:T);
%H((K+1)*T+(k-1)*T+1:((K+1)*T+k*T),T+(k-1)*T+1:(T+k*T)) =                   1/K* hessian(2*T+1:3*T,T+1:2*T);
%H((K+1)*T+(k-1)*T+1:((K+1)*T+k*T),(K+1)*T+(k-1)*T+1:((K+1)*T+k*T)) =       1/K* hessian(2*T+1:3*T,2*T+1:3*T);
%H((K+1)*T+(k-1)*T+1:((K+1)*T+k*T),(2*K+1)*T+(k-1)*T+1:((2*K+1)*T+k*T)) =   1/K* hessian(2*T+1:3*T,3*T+1:4*T);

%H((2*K+1)*T+(k-1)*T+1:((2*K+1)*T+k*T),1:(T)) = H((2*K+1)*T+(k-1)*T+1:((2*K+1)*T+k*T),1:(T))+1/K*hessian(3*T+1:4*T,1:T);
%H((2*K+1)*T+(k-1)*T+1:((2*K+1)*T+k*T),T+(k-1)*T+1:(T+k*T)) =                    1/K*hessian(3*T+1:4*T,T+1:2*T);
%H((2*K+1)*T+(k-1)*T+1:((2*K+1)*T+k*T),(K+1)*T+(k-1)*T+1:((K+1)*T+k*T)) =        1/K* hessian(3*T+1:4*T,2*T+1:3*T);
%H((2*K+1)*T+(k-1)*T+1:((2*K+1)*T+k*T),(2*K+1)*T+(k-1)*T+1:((2*K+1)*T+k*T)) =    1/K*hessian(3*T+1:4*T,3*T+1:4*T);


end



obj = 1/K .* sum(F(:,1)); % weighted (all weights=1/K)x sum of F(x,\tilde{x}^k)
grad = 1/K .* sum(G,1); % weighted sum over all gradients
%hessian = H;



end

