function Im = multinmf_recons_im(X,W,H,part)

%
% Reconstructs source images from nmf factorization (conservative)
% 
% Im = nmf_recons_im(X,W,H,part)
% 
%
% Input:
%   - X: truncated STFT of multichannal mixture (F x N x n_c)
%   - W: basis                                  (F x K)
%   - H: activation coef.                       (K x N)
%   - part : component indices
%
% Output:
%   - Im : reconstructed source images (F x N x n_s x n_c)
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

[F,N,n_c] = size(X);
n_s = size(part,2);

P = zeros(F,N,n_s);

for j=1:n_s
    P(:,:,j) = W(:,part{j}) * H(part{j},:);
end

Im = zeros(F, N, n_s, n_c);

if n_c == 1
    % use multinmf_recons_im for n_c = 2
    for j=1:n_s
        Im(:,:,j,1) = (P(:,:,j) ./ sum(P, 3)) .* X(:,:,1);
    end;
else
    error('nmf_recons_im: number of channels must be 1');
end
