function [Q, W, H] = init_ntf_mu_kl_pow (X, K, init_rnd, init_iterations, after_iterations)

    %
    % [Q, W, H] = init_ntf_mu_kl_pow(X, K);
    %
    % Initialize NMF from STFTs of mixture using Kullback-Leibler (KL) divergence
    %
    %
    % input 
    % -----
    %
    % X                 : half STFTs of separated sources [F x N x nchann]
    % K                 : number of NMF components (NMF_CompPerSrcNum * nsrc)
    % init_rnd          : number of random initializations
    % init_iterations   : number of iterations for each random initialization
    % after_iterations  : number of iterations for best initialization after selection
    %
    % output
    % ------
    %
    % Q                 : estimated mixing matrix [nchann x K]
    % W                 : estimated matrix of bases [F x K]
    % H                 : estimated matrix of contributions [K x N]
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Modified by Igor Chame 
    % @ apr 2016
    % igorchame@poli.ufrj.br
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % some constants
    small_noise_var = 3e-11;       % small noise variance

    [F, N, nchann] = size(X);

    V = abs(X).^2 + randn(F, N, nchann).^2 * small_noise_var;

    % run multiplicative algorithm with init_iterations iterations and init_rnd random
    % initializations, and keep the best initialization
    for init_ind = 1:init_rnd
       
        Q_ini = abs(randn(nchann,K)) + ones(nchann,K);
        W_ini = abs(randn(F,K)) + ones(F,K);
        H_ini = abs(randn(K,N)) + ones(K,N);

        [Q,W,H,E] = ntf_mu_is(V, init_iterations, Q_ini, W_ini, H_ini);

        E_end = E(end);

        if init_ind > 1
            if E_end < E_end_best
                Q_best = Q;
                W_best = W;
                H_best = H;
                E_end_best = E_end;
            end;
        else
            Q_best = Q;
            W_best = W;
            H_best = H;        
            E_end_best = E_end;
        end;
    end;

    if after_iterations == 0,
        Q = Q_best; W = W_best; H = H_best;
    else
        % now run after_iterations iterations multiplicative algorithm with best initialization
        [Q, W, H, E] = ntf_mu_is(V, after_iterations, Q_best, W_best, H_best);
    end

    Q = max(Q, ones(nchann,K) .* small_noise_var * 100);
    W = max(W, ones(F,K) .* small_noise_var * 100);
    H = max(H, ones(K,N) .* small_noise_var * 10);
end