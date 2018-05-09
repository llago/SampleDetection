function [Q, W, H, cost, SDRi, ISRi, SIRi, SARi, permi] = ...
    multinmf_inst_mu_eval(V, n_iter, Q, W, H, part, LOG, i_sim, X, nbin, mix_nsamp, fs, wav_file_prefix, switch_Q, switch_W, switch_H)

% Multichannel NMF minimizing Itakura-Saito divergence through multiplicative updates
% with linear instantaneous mixing
%
% [Q, W, H, cost] = multinmf_inst_mu(V, n_iter, Q, W, H, part, LOG, switch_Q, switch_W, switch_H)
%
% Input:
%   - V: positive matrix data       (F x N x n_c)
%   - n_iter: number of iterations
%   - init_Q: mixing matrix         (n_c x n_s)
%   - init_W: basis                 (F x K)
%   - init_H: activation coef.      (K x N)
%   - part : component indices
%   - LOG : log4m object
%   - switch_W, switch_H, switch_Q: (opt) switches (0 or 1) (def = 1)
%
% Output:
%   - Estimated Q, W and H
%   - Cost through iterations betw. data power and fitted variance.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2010 Cedric Fevotte
% (cedric.fevotte -at- telecom-paristech.fr)
%
% This software is distributed under the terms of the GNU Public License
% version 3 (http://www.gnu.org/licenses/gpl.txt)
%
% If you use this code please cite this paper
%
% A. Ozerov and C. Fevotte,
% "Multichannel nonnegative matrix factorization in convolutive mixtures for audio source separation,"
% IEEE Trans. on Audio, Speech and Lang. Proc. special issue on Signal Models and Representations
% of Musical and Environmental Sounds, vol. 18, no. 3, pp. 550-563, March 2010.
% Available: http://www.irisa.fr/metiss/ozerov/Publications/OzerovFevotte_IEEE_TASLP10.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modified by Igor Chame 
% @ apr 2016
% igorchame@poli.ufrj.br
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 14 || isempty(switch_Q)
    switch_Q = 1;
end;

if nargin < 15 || isempty(switch_W)
    switch_W = 1;
end;

if nargin < 16 || isempty(switch_H)
    switch_H = 1;
end;


[F,N,n_c] = size(V);
n_s = size(Q,2);
nsrc = n_s;

%% Definitions %%
V_ap    = zeros(F,N,n_c);
cost    = zeros(1,n_iter+1);

%% Compute app. variance structure V_ap %%
for j=1:n_s
    P_j = W(:,part{j}) * H(part{j},:);
    for i=1:n_c
        V_ap(:,:,i) = V_ap(:,:,i) + Q(i,j) .* P_j;
    end
end

cost(1) = sum(V(:)./V_ap(:) - log(V(:)./V_ap(:))) - F*N*n_c;

% initialize variables to evaluate estimation
results_to_save = [1e0 1e1 2e1 3e1 5e1 1e2 2e2 5e2 7e2 1e3 1e4 1e5 1e6]; % results to save
Q_MU_conv = zeros(nbin, 2, nsrc);
SDRi = zeros(nsrc,n_iter);
ISRi = zeros(nsrc,n_iter);
SIRi = zeros(nsrc,n_iter);
SARi = zeros(nsrc,n_iter);
permi = zeros(nsrc,n_iter);

for iter = 2:(n_iter+1)
    %%% Update Q %%%
    if switch_Q
        for j=1:n_s
            P_j = W(:,part{j}) * H(part{j},:);
            for i=1:n_c
                Q_old  = Q(i,j);
                Q(i,j) = Q(i,j) * sum(sum(V_ap(:,:,i).^-2.*P_j.*V(:,:,i))) / sum(sum(V_ap(:,:,i).^-1 .* P_j));
                V_ap(:,:,i) = V_ap(:,:,i) + (Q(i,j)-Q_old) .* P_j;
            end
        end
    end

    %%% Update W %%%
    if switch_W
        for j=1:n_s
            Kj   = length(part{j});
            Wnum = zeros(F,Kj);
            Wden = zeros(F,Kj);
            for i=1:n_c
                Wnum = Wnum + Q(i,j) * ((V_ap(:,:,i).^-2 .* V(:,:,i)) * H(part{j},:)');
                Wden = Wden + Q(i,j) * (V_ap(:,:,i).^-1 * H(part{j},:)');
            end

            Wj_old = W(:,part{j});
            W(:,part{j}) = W(:,part{j}) .* (Wnum./Wden);

            for i=1:n_c
                V_ap(:,:,i) = V_ap(:,:,i) + Q(i,j) * ((W(:,part{j})-Wj_old)*H(part{j},:));
            end
        end
    end

    %%% Update H %%%
    if switch_H
        for j=1:n_s
            Kj   = length(part{j});
            Hnum = zeros(Kj,N);
            Hden = zeros(Kj,N);
            for i=1:n_c
                Hnum = Hnum + (Q(i,j) * W(:,part{j})') * (V_ap(:,:,i).^-2 .* V(:,:,i));
                Hden = Hden + (Q(i,j) * W(:,part{j})') * V_ap(:,:,i).^-1;
            end

            Hj_old = H(part{j},:);
            H(part{j},:) = H(part{j},:) .* (Hnum./Hden);

            for i=1:n_c
                V_ap(:,:,i) = V_ap(:,:,i) + Q(i,j) * (W(:,part{j})*(H(part{j},:)-Hj_old));
            end

        end

    end

    cost(iter) = sum(V(:)./V_ap(:) - log(V(:)./V_ap(:))) - F*N*n_c;

    %%% Normalize %%%

    %% Scale Q / W %%
    if switch_Q && switch_W
        scale = sum(Q,1);
        Q = Q ./ repmat(scale,n_c,1);
        for j=1:n_s
            W(:,part{j}) = W(:,part{j}) * scale(j);
        end
    end
    %% Scale W / H %%
    if switch_W && switch_H
        scale = sum(W,1);
        W = W ./ repmat(scale ,F,1);
        H = H .* repmat(scale',1,N);
    end

    log_entry = sprintf(' iteration %d of %d, cost = %f', (iter-1), n_iter, cost(iter));;
    LOG.trace('multinmf_inst_mu',log_entry);
    %fprintf('MU update: iteration %d of %d, cost = %f\n', iter, n_iter, cost(iter));

    % Reconstruction of the spatial source images
    for f = 1:nbin
        Q_MU_conv(f,:,:) = Q;
    end;
    Ie_MU = multinmf_recons_im(X, Q_MU_conv, W, H, part);
    ie_MU=istft_multi(Ie_MU,mix_nsamp);
    [SDRi(:,(iter-1)),ISRi(:,(iter-1)),SIRi(:,(iter-1)),SARi(:,(iter-1)),permi(:,(iter-1))]=bss_eval_images(ie_MU,i_sim);
    % write wav file on selected iterations
    if ((iter-1) == n_iter || any(((iter-1))==results_to_save))
        for j=1:nsrc,
            wavwrite(reshape(ie_MU(j,:,:),mix_nsamp,2),fs,[wav_file_prefix 'src' int2str(permi(j,(iter-1))) '_' int2str((iter-1))  '.wav']);
        end
    end

end

