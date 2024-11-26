function [x_T,Tstop,rel_err,meas_err] = WirtingerFlowAlgorithms(A,y,x0,T,x,algorithm_version,rng_number,reg_par)

% Input:
%   A                   m times n measurement matrix 
%   y                   m times 1 vector containing the measurements
%   T                   maximum number of iterations 
%   x                   ground-truth object
%   algorithm_version   adressing the chosen version of gradient descent algorithm 
%   rng_number          seed for random experiment
%   reg_par             regularization parameter
%
% Output:
%   x_T         T-th iterate (reconstruction result)
%   Tstop       last iteration, either T or iteration before termination
%   rel_err     relative reconstruction error in all iterations
%   meas_err    measurement error in terms of the respective loss function


m = size(y,1);
n = size(x,1);
rel_err = zeros(T,1);
meas_err = zeros(T,1);

% Set seed
rng(rng_number)

% Initialization
xt = x0;

% Step size selection
[~,eig,~] = svd(A'*A);
mu = 1/max(max(eig)); mu = mu/2;
mu_AF = mu;


if algorithm_version == 0
    
    epsilon = 10^2;
    alpha = 2/3;
    y_max = max(y);
    normAr = zeros(1,m);
    for r = 1:m
        normAr(r) = norm(A(r,:));
    end 
    a_max = max(normAr);

    L0 = 2*y_max*a_max^2;
    L1 = 9/4*a_max^(4/3);
    K0 = L0*(2^(alpha^2/(1-alpha))+1);
    K1 = L1*2^(alpha^2/(1-alpha))*3^alpha;
    K2 = L1^(1/(1-alpha))*2^(alpha^2/(1-alpha))*3^alpha*(1-alpha)^(alpha/(1-alpha));
    gamma = epsilon^alpha/(12*(K0+K1+2*K2)+1);

elseif algorithm_version == 1

    reg_p = reg_par;
    mat_p = zeros(size(A));
    for i = 1:size(A,1)
     mat_p(i,:) = sqrt(1 + y(i)/(8*reg_p))*A(i,:); 
    end
    [~,eig,~] = svd(mat_p'*mat_p);
    mu_p = 1/max(max(eig));

end

% Run gradient descent in the chosen variant
for t = 1:T 

    y_xt = A*xt; 
    
    if algorithm_version == 0 % Wirtinger flow  
        grad  =    A'*(2*( abs(y_xt).^2- y) .* y_xt);  
        mu = gamma/(norm(grad)^(2/3));
    elseif algorithm_version == 1 % Poisson flow
        reg_p = reg_par;
        grad  = 2*A'*((1 - (y+0)./(abs(y_xt).^2+reg_p)).*(y_xt)); 
        mu = mu_p;
    elseif algorithm_version == 2 % Truncated Wirtinger flow
        alpha_z_lb = 0.1;
        alpha_z_ub = 5;
        alpha_h = 3;
        E1_1 = (alpha_z_lb <= sqrt(n) ./ vecnorm(A,2,2) .* abs(y_xt) / norm(xt));
        E1_2 = (sqrt(n) ./ vecnorm(A,2,2) .* abs(y_xt) / norm(xt) <= alpha_z_ub);
        Kt = 1/m * sum(abs(y - abs(y_xt).^2));
        E2 = (abs(y - abs(y_xt).^2) <= alpha_h * Kt * sqrt(n) ./ vecnorm(A,2,2) .* abs(y_xt) / norm(xt));     
        truncation = E1_1 .* E1_2 .* E2;
        grad =  - A'*((y - abs(y_xt).^2)./(conj(y_xt) ) .* truncation);
        mu = 2 * (0.2) /m;
    elseif algorithm_version == 3 % Flow with improved variance stabilizing transform without adaption for zero measurements
        c1 = 0.12;
        c2 = 0.27;
        grad_vec = (1/2*1/(2*0.25))*1/2*((sqrt(abs(y_xt).^2 + c1) - sqrt(y+c1)) + (sqrt(abs(y_xt).^2 + c2) - sqrt(y+c2))).* (1./sqrt(abs(y_xt).^2 + c1) + 1./sqrt(abs(y_xt).^2 + c2) ).*y_xt;
        grad = A'*grad_vec; 
        mu = mu_AF*2;
        mu = mu/(1/2*(3+sqrt(c2/c1)));
    elseif algorithm_version == 4 % Flow with improved variance stabilizing transform 
        c1 = 0.12;
        c2 = 0.27;
        grad_vec = 1/2*((sqrt(abs(y_xt).^2 +c1) - sqrt(y+c1)) + (sqrt(abs(y_xt).^2 +c2) - sqrt(y+c2))).* (1./sqrt(abs(y_xt).^2 + c1) + 1./sqrt(abs(y_xt).^2+c2) ).*y_xt;
        grad_vec_zeros =  y_xt;
        idx = y == 0;
        grad_vec(idx) = grad_vec_zeros(idx);
        grad = A'*grad_vec; 
        mu = mu_AF*2;        
        mu = mu/(1/2*(3+sqrt(c2/c1)));
    elseif algorithm_version == 5 % Amplitude flow
        if ~exist('reg_par','var')
            eps_af = 0; 
        else
            eps_af = reg_par;
        end
        grad =  2*A'*((1 - sqrt(y)./sqrt(abs(y_xt).^2 + eps_af)).*y_xt);       
    end
    
    % Gradient descent update
    xt = xt - mu * grad;
   
    % Compute relative reconstruction error
    theta_opt = fminbnd(@(theta) norm(exp(1.0i*theta)*xt-x),0,2*pi);
    rel_err(t) =norm(exp(1.0i*theta_opt)*xt-x)/norm(x);

    % Compute measurement error in terms of the respective loss function
    if algorithm_version == 0
        meas_err(t) = sum((y-abs(A*xt).^2).^2);
    elseif algorithm_version == 1
        meas_err(t) = sum(abs(A*xt).^2 - y.*log(abs(A*xt).^2 + reg_par));
    elseif algorithm_version == 2
        meas_err(t) = sum(abs(A*xt).^2 - y.*log(abs(A*xt).^2));  
    elseif algorithm_version == 3
        c1 = 0.12;
        c2 = 0.27; 
        approx_meas = abs(A*xt).^2;
        meas_err(t) = sum(1/2*(sqrt(y+c1) + sqrt(y+c2) - sqrt(approx_meas+c1) - sqrt(approx_meas+c2)).^2);    
    elseif algorithm_version == 4
        c1 = 0.12;
        c2 = 0.27; 
        y_0 = y == 0;
        approx_meas = abs(A*xt).^2;
        meas_err(t) = sum(1/2*(sqrt(y(~y_0)+c1) + sqrt(y(~y_0)+c2) - sqrt(approx_meas(~y_0)+c1) - sqrt(approx_meas(~y_0)+c2)).^2) + sum(approx_meas(y_0));    
    elseif algorithm_version == 5
        if ~exist('reg_par','var')
            eps_af = 0;  
        else
            eps_af = reg_par;
        end
        meas_err(t) = sum((sqrt(y)-sqrt(abs(A*xt).^2 + eps_af )).^2);
    end
     
    % Stopping criterion
    if t > 1 
       if abs(meas_err(t) - meas_err(t-1)) < 10^(-6) && meas_err(t) - meas_err(t-1) < 0      
        Tstop = t;
        break;
       end
    end
    
    Tstop = t;

end

theta_opt = fminbnd(@(theta) norm(exp(1.0i*theta)*x0-x),0,2*pi);
rel_err = [norm(exp(1.0i*theta_opt)*x0-x)/norm(x); rel_err];

% Reconstruction result
theta_opt_t = fminbnd(@(theta) norm(exp(1.0i*theta)*xt-x),0,2*pi);
x_T = exp(1.0i*theta_opt_t)*xt;

end