function [W, H, cost] = is_nmf_mu(V, W, H, n_iter)

    % Itakura-Saito NMF with multiplicative updates
    %
    % [W, H, cost] = is_nmf_mu(V, W_ini, H_ini, n_iter)
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

    cost = zeros(1,n_iter);

    % Avoid zeros that may unable fatoration
    eps = 2.2204e-32;   % very small constant  

    % Compute data approximate
    V_ap = W*H+eps;

    % Compute initial cost value
    cost(1) = sum(V(:)./(V_ap(:)+eps) - log(V(:)./(V_ap(:)+eps)+eps)) - F*N;

    h = waitbar(0,'MU/NMF-IS');
    for iter = 2:n_iter, waitbar(iter/n_iter,h)

        % Update W
        W = W .* ((V.*(V_ap+eps).^-2)*H')./((V_ap+eps).^-1*H'+eps);
        V_ap = W*H;

        % Update H
        H = H .* (W'*(V.*(V_ap+eps).^-2))./(W'*(V_ap+eps).^-1+eps);
        V_ap = W*H;

        % Norm-2 normalization
        scale = sqrt(sum(W.^2,1)); 
        W = W .* repmat(scale.^-1,F,1);
        H = H .* repmat(scale',1,N);

        % Compute cost value
        cost(iter) = sum(V(:)./(V_ap(:)+eps) - log(V(:)./(V_ap(:)+eps)+eps)) - F*N;

        % fprintf('MU/NMF-IS: iteration %d of %d, cost = %f\n', iter, n_iter, cost(iter));

    end
    close(h);
end
