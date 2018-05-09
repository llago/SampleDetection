function [part, L, D, W, H, Sigma_b] = init_ntf_em_is_pow (X, K, nsrc, init_rnd, init_iterations, after_iterations)

    %
    % [L, D, W, H] = init_ntf_em_is_pow
    %
    % Initialize NMF from STFTs of mixture using EM IS
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

    V = abs(X).^2;
    mix_psd = (1/nchann) * (mean(sum(V,3),2));
    Sigma_b_ini = mix_psd / 100;
    %Sigma_b_ini = zeros(size(mix_psd));
    
    % Make a labelling matrix
    L = zeros(nsrc,K);
    part = cell(1,nsrc);
    for j = 1:nsrc,
        part{j} = (j-1)*(K/nsrc)+1: j*(K/nsrc);
        L(j,part{j}) = ones(1,length(part{j}));
    end;

    % run multiplicative algorithm with init_iterations iterations and init_rnd random
    % initializations, and keep the best initialization
    for init_ind = 1:init_rnd
       
        Q_ini = abs(randn(nchann,K)) + ones(nchann,K);
        W_ini = abs(randn(F,K)) + ones(F,K) .* (mix_psd * ones(1,K));
        H_ini = abs(randn(K,N)) + ones(K,N);
        D_ini = 0.5*(1.9*abs(randn(nchann,nsrc)) + 0.1*ones(nchann,nsrc));

        [W, H, D, Sigma_b, Se_EM, E] = multinmf_inst_em(X, W_ini, H_ini, D_ini, Sigma_b_ini, part, init_iterations);

        Q = D * L;

        E_end = E(end);

        if init_ind > 1
            if E_end > E_end_best
                Q_best = Q;
                W_best = W;
                H_best = H;
                D_best = D;
                Sigma_b_best = Sigma_b;
                E_end_best = E_end;
            end;
        else
            Q_best = Q;
            W_best = W;
            H_best = H;  
            D_best = D; 
            Sigma_b_best = Sigma_b;     
            E_end_best = E_end;
        end;
    end;

    if after_iterations == 0,
        Q = Q_best; W = W_best; H = H_best; D = D_best; Sigma_b = Sigma_b_best;
    else
        % now run after_iterations iterations multiplicative algorithm with best initialization
        [W, H, D, Sigma_b, Se_EM, E] = multinmf_inst_em(X, W_best, H_best, D_best, Sigma_b_best, part, after_iterations);
        Q = D * L;
    end

    D = max(D, ones(nchann,nsrc) .* small_noise_var * 100);
    Q = max(Q, ones(nchann,K) .* small_noise_var * 100);
    W = max(W, ones(F,K) .* small_noise_var * 100);
    H = max(H, ones(K,N) .* small_noise_var * 10);
end