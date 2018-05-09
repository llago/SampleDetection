function [W, H, cost] = is_nmf_ig(V, n_iter, W, H, alpha_H)

% Itakura-Saito NMF with inverse-Gamma Markov chain prior on H
%
% [W, H, cost] = is_nmf_ig(V, n_iter, W_ini, H_ini, alpha_H)
%
% Inputs:
%   - V: positive matrix data (F x N)
%   - n_iter: number of iterations
%   - W_ini: basis matrix (F x K)
%   - H_ini: gains matrix (K x N)
%   - alpha_H: IG Markov chains shape parameters (K x 1)
%
% Outputs :
%   - W and H such that
%
%               V \approx W * H
%
%   - cost : IS divergence though iterations
%
% If you use this code please cite this paper
%
% C. Fevotte, N. Bertin and J.-L. Durrieu. "Nonnegative matrix factorization
% with the Itakura-Saito divergence. With application to music analysis,"
% Neural Computation, vol. 21, no 3, Mar. 2009
%
% Report bugs to Cedric Fevotte
% fevotte -at- telecom-paristech.fr
% Checked 10/02/09

[F,N] = size(V);
K = size(W,2);

cost = zeros(1,n_iter);

% Compute data approximate
V_ap = W*H;

% Compute initial cost value
cost(1) = sum(V(:)./V_ap(:) - log(V(:)./V_ap(:))) - F*N;

for iter=2:n_iter

    for k=1:K

        % Power of component C_k %
        PowC_k = W(:,k) * H(k,:);

        % Power of residual
        PowR_k = V_ap - PowC_k;

        % Wiener gain
        G_k = PowC_k ./ V_ap;

        % Posterior power of component C_k %
        V_k = G_k .* (G_k .* V + PowR_k);

        %% Update H %%

        % ML estimate C_k/F
        C_k = (W(:,k).^-1)'*V_k;

        % MAP estimate (see paper)

        % n = 1 %
        p2 = (alpha_H(k)+1)/H(k,2);
        p1 = F - alpha_H(k) + 1;
        p0 = - C_k(1);
        H(k,1) = (sqrt(p1^2-4*p2*p0) - p1)/(2*p2);

        % n = 2:N-1 %
        for n=2:N-1
            p2 = (alpha_H(k)+1)/H(k,n+1);
            p1 = F + 1;
            p0 = - C_k(n) - (alpha_H(k)+1)* H(k,n-1) ;
            H(k,n) = (sqrt(p1^2-4*p2*p0) - p1)/(2*p2);
        end

        % n = N %
        H(k,N) = (C_k(N) + (alpha_H(k)+1)* H(k,N-1))/(F+alpha_H(k)+1);

        %% Update W %%

        % ML estimate
        D_k = V_k*(H(k,:).^-1)';
        W(:,k) = D_k/N;

        %% Norm-2 normalization %%
        scale = sqrt(sum(W(:,k).^2));
        W(:,k) = W(:,k)/scale;
        H(k,:) = H(k,:)*scale;

        %% Update data approximate %%
        V_ap = PowR_k + W(:,k) * H(k,:);

    end

    % Compute cost value (this is IS divergence, not MAP criterion)
    cost(iter) = sum(sum(V./V_ap - log(V./V_ap))) - F*N;

end

%% Plugin the following piece of code instead of "W(:,k) = D_k/N;"
%% to enforce regularization of the columns of W with an Inverse-Gamma 
%% Markov chain prior, similarly to H.
%
% alpha_W = 10 * ones(K,1);
%
% % f = 1 %
% p2 = (alpha_W(k)+1)/W(2,k);
% p1 = N - alpha_W(k) + 1;
% p0 = - D_k(1);
% W(1,k) = (sqrt(p1^2-4*p2*p0) - p1)/(2*p2);
% 
% % f = 2:F-1 %
% for f=2:F-1
%     p2 = (alpha_W(k)+1)/W(f+1,k);
%     p1 = N + 1;
%     p0 = - D_k(f) - (alpha_W(k)+1)* W(f-1,k);
%     W(f,k) = (sqrt(p1^2-4*p2*p0) - p1)/(2*p2);
% end
% 
% % f = F %
% W(F,k) = (D_k(F) + (alpha_W(k)+1)* W(F-1,k))/(N+alpha_W(k)+1);
