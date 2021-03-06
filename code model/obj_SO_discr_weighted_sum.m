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
%H = zeros(T+3*K*T,T+3*K*T);

for k = 1:K 
x_tilde = E(k,:)+0.95*x(((2*K+k)*T+1):((1+2*K+k)*T))...
    -x(((K+k)*T+1):((1+K+k)*T)); % compute \tilde{x}^k as e^k + 0.95 \tilde{b^out,k} - \tilde{b^in,k}
[F(k),grad_SOC_b,grad_x,~] = obj_SO_discr(x(1:T),x_tilde,cost,penalty,penalty_grad,epsilon,P,penalty_hess);

% In the "k^th objective" F(k) only b^{in,k} and b^{out,k} appear and thus
% only the corresponding entries in the gradient row G(k,:) are non-zero:
G(k,1:T) = grad_x; % gradient w.r.t. x
G(k,(k*T+1):((1+k)*T)) = grad_SOC_b(1:T); % gradient w.r.t. SOC^k
G(k,((K+k)*T+1):((K+1+k)*T)) = grad_SOC_b(T+1:2*T); % gradient w.r.t. b^in,k
G(k,((2*K+1+k-1)*T+1):((2*K+1+k)*T)) = grad_SOC_b(2*T+1:3*T); % gradient w.r.t. b^out,k

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

obj = 1/K * sum(F(:,1)); % weighted (all weights=1/K)x sum of F(x,\tilde{x}^k)
grad = 1/K * sum(G,1); % weighted sum over all gradients
hessian = 0;

end

% function [obj,grad_SOC_b,grad_x,hessian] = obj_SO_discr_in(x,e,cost,penalty,penalty_grad,epsilon,P,penalty_hess)
% %function [obj,grad] = obj_SO_discr(x,e,cost,penalty,penalty_grad,P)
% % Calculates the objective function -F(x,E) as in 1.5 Objective function
% % Since F(x,E) = \sum_{i=1}^T F^{(i)}(x_i,E_i) (see also 1.5), obj_SO_discr
% % calculates the revenue values F^{(i)}(x_i,E_i) for all time steps 
% % and returns the sum of all values
% % 
% % Input:
% %       x - n by m matrix                x(i,j) corresponds to the j's hour of
% %                                        the i's day
% %       e - 1 by m vector                e(j) is the solar radiance during the
% %                                        i's hour
% %       cost - 1 by m vector             price for an energy unit
% %       penalty - function handle        penalty function, should be able to
% %                                        act on matrices
% %       penalty_grad function handle     derivative of the penalty function
% %       epsilon - pos. scalar << 1       defines the width of the no-penalty
% %                                        interval [x(1-epsilon),x(1+epsilon)]
% %       P - scalar                       nominal power of PV element
% %       var - text                       variables w.r.t. which the
% %                                        gradient is calculated
% % Output:
% %       obj - n by 1 vector              obj(i) = sum_j
% %                                               F(x(i,j),e(j))
% %       grad_SOC_b - n by 3m matrix      grad w.r.t. [SOC,b^{in},b^{out}]
% %       grad_x - n by m matrix           grad w.r.t. x
% %       hessian - 4m by 4m matrix        hessian of the objective w.r.t. [x,SOC,b^{in]},b^{out}]
% %% Check the size of x and e
% s = size(x); T = s(2); 
% if size(e,2) ~= s(2), error('Sizes of x and e do not match'); end
% if s(1) ~= 1, error('x schedules for than one day'); end
% %% Computation of the objective
% E = ones(s(1),1) * e;                        % E is an n by m matrix with e in each row
% Cost = ones(s(1),1) * cost;                  % C is an n by m matrix with cost in each row
% [yp,gradp,~,~] = smooth_ppart_in((x - epsilon*P) - E,1E-3,5);
% [yc,gradc,~,~] = smooth_ppart_in(E - (x + epsilon*P),1E-3,5);
% 
% obj = penalty(yp) + yc.*Cost - E.*Cost;  
% obj = sum(obj,2); % sum of all single revenue values
% 
% %% Computation of gradient w.r.t x
% grad_x = penalty_grad(yp).*gradp - gradc.*Cost; 
% 
% %% Computation of gradient w.r.t [SOC,b^in,b^out]
% Cost = ones(s(1),1) * [cost, cost, cost]; % C is an n by 3m matrix with three times cost in each row
% grad_x_tilde = [zeros(s(1),T), -1 * ones(s(1),T), 0.95 * ones(s(1),T)];
% grad_SOC_b = - penalty_grad([yp,yp,yp]).*[gradp, gradp, gradp].*grad_x_tilde + grad_x_tilde.*Cost.*[gradc, gradc, gradc] - grad_x_tilde.*Cost ;
% 
% %% Computation of the hessian w.r.t. [x,SOC,b^in,b^out]
% hessian = 0;
% % hessian(1:T,1:T) = ones(T,1)*(penalty_hess(yp).*gradp.*gradp + penalty_grad(yp).*hessp + cost.*hessp);
% % hessian(1:T,2*T+1:3*T) = ones(T,1)*(penalty_hess(yp).*gradp.*1.*gradp + penalty_grad(yp).*hessp + cost.*hessp);
% % hessian(1:T,3*T+1:4*T) = ones(T,1)*(penalty_hess(yp).*gradp.*(-0.95).*gradp + penalty_grad(yp).*hessp + cost.*hessp);
% % hessian(2*T+1:3*T,1:T) = hessian(1:T,2*T+1:3*T);
% % hessian(3*T+1:4*T,1:T) = hessian(1:T,3*T+1:4*T);
% % hessian(2*T+1:3*T,2*T+1:3*T) = ones(T,1)*(gradp.*penalty_hess(yp).*gradp);
% % hessian(2*T+1:3*T,3*T+1:4*T) = ones(T,1)*(- gradp.*penalty_hess(yp).*gradp.*0.95);
% % hessian(3*T+1:4*T,2*T+1:3*T) = hessian(T+1:2*T,2*T+1:3*T);
% % hessian(3*T+1:4*T,3*T+1:4*T) = ones(T,1)*(0.95.*0.95.*gradp.*penalty_hess(yp).*gradp);
% %% Alter Code:
% %obj = penalty(max(zeros(s), (x - epsilon*P) - E)) + ...
% %   max(zeros(s), E - (x + epsilon*P))*cost - E*cost;                     % just evaluating the formula
% end
% 
% function [y, grad, hess, p] = smooth_ppart_in(x,epsilon,N)
% % SMOOTH_PART(x,epsilon) berechnet den positiven Teil von x, falls
% % |x|>=epsilon, oder epsilon/4 + 0.5*x + 0.25/epsilon*x^2, falls
% % |x|<epsilon. Das entspricht der Gl?ttung der Funktion max(x,0) in dem
% % Intervall [-epsilon epsilon]. Here a quadratic polynomial interpolation
% % between the points -epsilon and epsilon is used.
% % 
% % Input:
% %       x - n by m matrix           
% %       epsilon - pos. scalar << 1            
% %       N - natural number > 5      can be used to produce twice
% %                                   differentiable approximation
% %
% % Output:
% %       y - n by m matrix           value of the smooth function at x
% %       grad - n by m matrix        value of the derivative at x
% %       hess - n by m matrix        value of the 2nd derivative at x 
% %   !!!! hess is only correctly implemented for N=2 !!!!
% %       p - 1 by 3 or
% %           1 by N+1 vector         vector of coefficients from the
% %                                   interpolating polynomial
% 
% s = size(x);
% 
% %% the first case is only relevant, if a twice differentiable function is needed
% if nargin == 3,
%     A = [epsilon.^(0:N);...
%         (-epsilon).^(0:N);...
%         (0:N).*[1 epsilon.^(0:N-1)];...
%         (0:N).*[1 (-epsilon).^(0:N-1)];...
%         (-1:N-1).*(0:N).*[1 1 epsilon.^(0:N-2)];...
%         (-1:N-1).*(0:N).*[1 1 (-epsilon).^(0:N-2)]
%         ];
%     A = A(:,1:end-1);
%     b = [epsilon; 0; 1; 0; 0; 0] - epsilon^N;
%     p = linsolve(A,b);
%     p = [1; p(end:-1:1)]';
% 
% elseif nargin == 2, N = 2; p = [0.25/epsilon 0.5 epsilon/4]; % quadratic interpolation with precomputed polynomial
% end
% 
% y = zeros(s);
% grad = zeros(s);
% hess = zeros(s);
% 
% %% |x|<eps
% I = (x > -epsilon) & (x < epsilon);
% y(I) = polyval(p,x(I));
% grad(I) = polyval((N:-1:1).*p(1:end-1), x(I));
% hess(I) = polyval((N:-1:2).*(N-1:-1:1).*p(1:end-2), x(I));
% %% x<-eps
% I = (x <= -epsilon);
% y(I) = 0;
% grad(I) = 0;
% hess(I) = 0;
% %% x>eps
% I = (x >= epsilon);
% y(I) = x(I);
% grad(I) = 1;
% hess(I) = 0;
% end
% 
