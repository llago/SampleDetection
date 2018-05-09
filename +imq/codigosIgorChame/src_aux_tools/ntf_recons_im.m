function Im = ntf_recons_im(X,Q,W,H,part)

% Reconstructs component images from NTF decomposition through
% "Wiener filtering" of each channel.
%
% Input:
%   - X: truncated STFT of multichannal mixture (F x N x n_c)
% 	- Q: Mixing matrix (I x K)
%   - W: basis                                  (F x K)
%   - H: activation coef.                       (K x N)
%   - part : component indices
%
% Output:
%   - Im : reconstructed source images (F x N x n_s x n_c)

X = permute( X, [3 1 2]);

[I,F,N] = size(X);
H = H';
K = size(W,2);
n_s = size(part,2);


QoW = zeros(I,F,K);
for k=1:K,
    QoW(:,:,k) = Q(:,k)*W(:,k)';
end;

V_ap = zeros(I,F,N);
V_est = zeros(n_s,I,F,N);
for i=1:I,
    V_ap(i,:,:) = squeeze(QoW(i,:,:)) * H';
    for j=1:n_s,
    	V_est(j,i,:,:) = reshape( (squeeze(QoW(i,:,part{j})) * H(:,part{j})'), 1, 1, F, N);
    end;
end;

Im = zeros(F,N,n_s,I);
for j=1:n_s,
    Im(:,:,j,:) = permute( ( squeeze(V_est(j,:,:,:)) ./ V_ap) .* X, [2 3 4 1]);
end;
