function [part, L, D, W, H] = init_ntfclus_mu_kl_pow (X, K, nsrc, init_rnd, init_iterations, after_iterations)

    %
    % [L, D, W, H] = init_ntfclus_mu_kl_pow(X, K);
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

    % Make a labelling matrix
    L = zeros(nsrc,K);
    part = cell(1,nsrc);
    for j = 1:nsrc
        part{j} = (j-1)*(K/nsrc)+1: j*(K/nsrc);
        L(j,part{j}) = ones(1,length(part{j}));
    end

    % run multiplicative algorithm with init_iterations iterations and init_rnd random
    % initializations, and keep the best initialization
    for init_ind = 1:init_rnd
       
        Q_ini = abs(randn(nchann,K)) + ones(nchann,K);
        W_ini = abs(randn(F,K)) + ones(F,K);
        H_ini = abs(randn(K,N)) + ones(K,N);
        D_ini = abs(randn(nchann,nsrc)) + ones(nchann,nsrc);

        [D, W, H, E] = ntfclus_mu_is(V, init_iterations, L, D_ini, W_ini, H_ini);
        Q = D * L;

        E_end = E(end);

        if init_ind > 1
            if E_end < E_end_best
                Q_best = Q;
                W_best = W;
                H_best = H;
                D_best = D;
                E_end_best = E_end;
            end;
        else
            Q_best = Q;
            W_best = W;
            H_best = H;  
            D_best = D;      
            E_end_best = E_end;
        end;
    end;

    if after_iterations == 0,
        Q = Q_best; W = W_best; H = H_best; D = D_best;
    else
        % now run after_iterations iterations multiplicative algorithm with best initialization
        [D, W, H, E] = ntfclus_mu_is(V, after_iterations, L, D_best, W_best, H_best);
    end

    D = max(D, ones(nchann,nsrc) .* small_noise_var * 100);
    Q = max(Q, ones(nchann,K) .* small_noise_var * 100);
    W = max(W, ones(F,K) .* small_noise_var * 100);
    H = max(H, ones(K,N) .* small_noise_var * 10);
end