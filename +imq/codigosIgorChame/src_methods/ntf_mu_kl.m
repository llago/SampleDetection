function [Q, W, H, cost] = ntf_mu_kl(V, n_iter, Q, W, H)

    % PARAFAC-NTF minimizing KL divergence through multiplicative updates
    %
    %
    % [Q, W, H, cost] = ntf_mu_kl(V, n_iter, Q_ini, W_ini, H_ini)
    %
    % Inputs :
    %
    % V : I x F x N nonnegative data tensor
    % n_iter : nb of algorithm iterations
    % Q_ini : positive matrix of dimensions I x K
    % W_ini : positive matrix of dimensions F x K
    % H_ini : positive matrix of dimensions K x N
    %
    % Outputs :
    %
    % Q : nonnegative mixing matrix
    % W : matrix of spectra
    % H : activation matrix
    % cost : cost function at every iteration
    %
    % See and cite paper
    %
    % C. Fevotte and A. Ozerov, "Notes on nonnegative tensor factorization of
    % the spectrogram for audio source separation~: statistical insights and
    % towards self-clustering of the spatial cues", Proc. 7th International
    % Symposium on Computer Music Modeling and Retrieval (CMMR'2010), 2010.
    
    % Author : C. Fevotte
    % Bugs & reports : fevotte -at- telecom-paristech.fr
    % June 2010

    % Modified by : Igor Chame
    % igorchame -at- poli.ufrj.br
    % May 2016

    V = permute (V, [3 1 2]);  % ntf_mu uses a diferent order of dimensions
    [I,F,N] = size(V);
    K = size(W,2);

    H = H';

    cost = zeros(1,n_iter);

    WoH = zeros(F,N,K);
    QoH = zeros(I,N,K);
    QoW = zeros(I,F,K);

    V_ap = zeros(I,F,N);

    % K-R product
    for k=1:K
        QoH(:,:,k) = Q(:,k)*H(:,k)';
    end

    % Tensor approximate
    for i=1:I
        V_ap(i,:,:) = W * squeeze(QoH(i,:,:))';
    end

    G_neg = V./V_ap;

    cost(1) = sum(V(:).*log(V(:)./V_ap(:)) - V(:) + V_ap(:));

    h = waitbar(0,'MU/NTF-KL');
    for iter = 2:n_iter, waitbar(iter/n_iter,h)

        %% Update Q %%

        % K-R product
        for k=1:K
            WoH(:,:,k) = W(:,k)*H(:,k)';
        end

        Q = Q .* (reshape(G_neg,I,F*N)*reshape(WoH,F*N,K))./(repmat(sum(reshape(WoH,F*N,K),1),I,1));

        % Tensor approximate
        for f=1:F
            V_ap(:,f,:) = Q * squeeze(WoH(f,:,:))';
        end

        G_neg = V./V_ap;


        %% Update W %%

        % K-R product
        for k=1:K
            QoH(:,:,k) = Q(:,k)*H(:,k)';
        end

        W = W .* (reshape(permute(G_neg,[2 1 3]),F,I*N)*reshape(QoH,I*N,K))./repmat(sum(reshape(QoH,I*N,K),1),F,1);

        % Tensor approximate
        for i=1:I
            V_ap(i,:,:) = W * squeeze(QoH(i,:,:))';
        end

        G_neg = V./V_ap;


        %% Update H %%

        % K-R product
        for k=1:K
            QoW(:,:,k) = Q(:,k)*W(:,k)';
        end

        H = H .* (reshape(G_neg,I*F,N)'*reshape(QoW,I*F,K))./repmat(sum(reshape(QoW,I*F,K),1),N,1);

        % Tensor approximate
        for i=1:I
            V_ap(i,:,:) = squeeze(QoW(i,:,:)) * H';
        end

        G_neg = V./V_ap;


        %% Compute cost %%
        cost(iter) = sum(V(:).*log(V(:)./V_ap(:)) - V(:) + V_ap(:));

        %% Solve scale ambiguities %%
        scale = sum(Q,1);
        Q = Q .* repmat(scale.^-1,I,1);
        W = W .* repmat(scale,F,1);

        scale = sum(W,1);
        W = W .* repmat(scale.^-1,F,1);
        H = H .* repmat(scale,N,1);

    end
    close(h);

    H = H';
end