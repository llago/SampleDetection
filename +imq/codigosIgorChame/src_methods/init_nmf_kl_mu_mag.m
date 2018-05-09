function [W, H] = init_nmf_kl_mu (X, K, init_rnd, init_iterations, after_iterations)

    %
    % [W, H] = init_nmf_kl_mu(X, K);
    %
    % Initialize NMF from STFTs of mixture using Kullback-Leibler (KL) divergence
    %
    %
    % input 
    % -----
    %
    % X                 : half STFTs of separated sources [F x N x n_c]
    % K                 : number of NMF components (NMF_CompPerSrcNum * nsrc)
    % init_rnd          : number of random initializations
    % init_iterations   : number of iterations for each random initialization
    % after_iterations  : number of iterations for best initialization after selection
    %
    % output
    % ------
    %
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

    [F, N, n_c] = size(X);

    V = randn(F, N) * small_noise_var;
    for chann = 1:n_c, V = V + abs(X(:,:,chann)); end;
    mix_psd = (1/n_c) * (mean(V,2));

    % run multiplicative algorithm with init_iterations iterations and init_rnd random
    % initializations, and keep the best initialization
    for init_ind = 1:init_rnd
        % W is intialized so that its enegy follows mixture PSD
        W = 0.5 * (abs(randn(F,K)) + ones(F,K)); %.* (mix_psd * ones(1,K));
        H = 0.5 * (abs(randn(K,N)) + ones(K,N));

        [W,H,E] = nmf_kl_mu(V,W,H,init_iterations,small_noise_var);

        E_end = E(end);

        if init_ind > 1
            if E_end < E_end_best
                W_best = W;
                H_best = H;
                E_end_best = E_end;
            end;
        else
            W_best = W;
            H_best = H;        
            E_end_best = E_end;
        end;
    end;

    if after_iterations == 0,
        W = W_best; H = H_best;
    else
        % now run after_iterations iterations multiplicative algorithm with best initialization
        [W, H] = nmf_kl_mu(V,W_best,H_best,after_iterations,small_noise_var);
    end

    W = max(W, ones(F,K) .* small_noise_var * 100);
    H = max(H, ones(K,N) .* small_noise_var * 10);
end