function [W, H, cost] = is_nmf_em(V, W, H, n_iter)

    % Itakura-Saito NMF with SAGE
    %
    % [W, H, cost] = is_nmf_em(V, n_iter, W_ini, H_ini)
    %
    % Inputs:
    %   - V: positive matrix data (F x N)
    %   - n_iter: number of iterations
    %   - W_ini: basis matrix (F x K)
    %   - H_ini: gains matrix (K x N)
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

    h = waitbar(0,'EM/NMF-IS');
    for iter=2:n_iter, waitbar(iter/n_iter,h)

        for k=1:K

            % Power of component C_k %
            PowC_k = W(:,k) * H(k,:);

            % Power of residual
            PowR_k = V_ap - PowC_k;

            % Wiener gain
            G_k = PowC_k ./ V_ap;

            % Posterior power of component C_k %
            V_k = G_k .* (G_k .* V + PowR_k);

            % Update row k of H
            H(k,:) = (W(:,k).^-1)'*V_k/F;

            % Update column k of W        
            W(:,k) = V_k*(H(k,:).^-1)'/N;

            % Norm-2 normalization
            scale = sqrt(sum(W(:,k).^2));
            W(:,k) = W(:,k)/scale;
            H(k,:) = H(k,:)*scale;

            % Update data approximate
            V_ap = PowR_k + W(:,k) * H(k,:);

        end %k

        % Compute cost value
        cost(iter) = sum(V(:)./V_ap(:) - log(V(:)./V_ap(:))) - F*N;

    end
    close(h);   
end