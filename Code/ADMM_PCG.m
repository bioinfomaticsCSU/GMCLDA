function [X, stat] = ADMM_PCG(m_train, m_val, m_test, prob_params)
    %% Problem parameters
    Ll = prob_params.Ll;
    Ld = prob_params.Ld;
    
    % the next are used instead of mask_train:
    AtA_op = prob_params.AtA_op;
    At_op = prob_params.At_op; 

    % the regularization parameters
    gamma_n = prob_params.gamma_n;  % nuclear norm
    gamma_l = prob_params.gamma_l;   % rows graph
    gamma_d = prob_params.gamma_d;  % columns graph

    size_X = prob_params.size_X;    % needed to go from m_train to M_train
    numel_X = prod(size_X);             % numel(X)

    %% PARAMETERS OF ADMM:
    TOL_ABS = 1e-2;
    TOL_REL = 1e-3; 
    maxit = 200;    % number of iteration
    rho_ADMM = 1;   % penalty parameter

    % for prox of nuclear norm
    param_nuclear.svds = false; 

    %% operators needed for STEP 2 of ADMM: (solve C(Y) = b)
    C_op = @(y) AtA_op(y) + vec(gamma_l * Ll * reshape(y, size_X) + gamma_d * reshape(y, size_X) * Ld + rho_ADMM * reshape(y, size_X));
    % to precondition Conjugate Gradients, we provide the inverse of the diagonal of the C operator matrix
    AtAM = reshape(At_op(m_train), size_X);

    %% ADMM Initialization
    Z = AtAM;
    if isa(m_train, 'single')   % ascertain type
        Y = single(zeros(size_X)); % is this a valid starting point?
    else
        Y = zeros(size_X);
    end
    Y_old = Y;
    
    res_pri  = nan(maxit,1);
    res_dual = nan(maxit,1);
    eps_pri  = nan(maxit,1);
    eps_dual = nan(maxit,1);
    
    pcg_tol = 1e-5;

    %% ADMM:
    for i = 1 : maxit  
        %  fprintf('iter %d:  ', i);    
        %% STEP 1: X <- argmin (gamma_n||X||_*  +  rho/2 ||X - (Y - Z)||_F^2)
        [X, info_nuclear] = prox_nuclearnorm(Y-Z, gamma_n/rho_ADMM, param_nuclear);
        param_nuclear.max_rank = round(max(5, info_nuclear.rank) * 1.5);
        
        %% STEP 2: Y <- argmin (    1/2 ||A(Y) - m||^2 
        %%                         + gamma_c/2 ||grad_Gc(Y)||_F^2
        %%                         + gamma_r/2 ||grad_Gr(Y)||_F^2 
        %%                         + rho/2    ||X - (Y - Z)||_F^2)
        [y, ~, ~, ~, resvec] = pcg(C_op, vec(AtAM + rho_ADMM * (X + Z)), pcg_tol, 30, [], [], vec(X));

        Y = reshape(y, size_X);
        
        Y(Y<0)=0;   % boundary
        Y(Y>1)=1;

        %% STEP 3: Z <- Z + X - Y
        Z = Z +1*( X - Y);    

        %% STOPPING CRITERION:
        res_pri(i) = norm(X - Y, 'fro');
        res_dual(i) = rho_ADMM * norm(Y_old - Y, 'fro');

        eps_pri(i)  = sqrt(numel_X) * TOL_ABS + TOL_REL * max(norm(X, 'fro'), norm(Y, 'fro'));
        eps_dual(i) = sqrt(numel_X) * TOL_ABS + TOL_REL * norm(Z, 'fro');

        if res_pri(i) < eps_pri(i) && res_dual(i) < eps_dual(i)
            break
        end
        
        Y_old = Y;

    end
    stat.last_iter = i;
end
