function [source_NMF_ind] = assoc_compToSrc_ideal(W,H,S_ref_STFT_V)

    % Ideal association between components and sources by
    % minimizing V distances
    %
    % [source_NMF_ind] = assoc_compToSrc_ideal(W,H,S_ref_STFT_V);
    %
    % Inputs:
    % ------
    % W and H: such that
    %               V \approx W * H
    %
    % S_ref_STFT_V: nsrc x nbin x nfram x nchann matrix containing STFT of source 
    % references
    %
    % Outputs:
    % -------
    % source_NMF_ind: nsrc x 1 vector containing the best association of components
    % with sources
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Modified by Igor Chame 
    % @ apr 2016
    % igorchame@poli.ufrj.br
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Errors %%%
    if nargin<3, error('Not enough input arguments.'); end
    [nbin_est,ncmp]=size(W);
    [ncmp2,nfram_est]=size(H);
    [nsrc,nbin,nfram,nchann]=size(S_ref_STFT_V);
    if nfram_est~=nfram, error('Source references and estimation must have same number of frames.'); end
    if nbin_est~=nbin, error('Source references and estimation must have same number of frequency bins.'); end
    if ncmp~=ncmp2, error('Matrices W and H must have same number of components.'); end
    if nchann~=1, error('Number of Channels must be 1.'); end

    %%% Estimate V of each component %%%
    for k=1:ncmp, V(k,:,:) = W(:,k) * H(k,:); end;

    %%% Calculate similarity between component V spectrum and source references %%%
    similarity = zeros(ncmp,nsrc);
    for k=1:ncmp,
        for j=1:nsrc,    
            similarity(k,j) = sum(sum( (S_ref_STFT_V(j,:,:)).^2 )) ./ ( sum(sum( (S_ref_STFT_V(j,:,:) - V(k,:,:)).^2 )) ); 
        end;
    end;

    %%% Associate components to source by maximizing similarity %%%
    source_NMF_ind = cell(1,nsrc);
    NMF_CompPerSrcNum = ncmp/nsrc;
    for j=1:nsrc,
        % Select NMF_CompPerSrcNum Components that are most alike this source j
        [ sorted_components, sorting_components_ind] = sort(similarity(:,j),'descend');
        maxValues = sorted_components(1:NMF_CompPerSrcNum);
        source_NMF_ind{j} = sorting_components_ind(1:NMF_CompPerSrcNum);
        % Zero similarity of the components that have been classified
        for src_counter=1:nsrc, for k=1:NMF_CompPerSrcNum, similarity(source_NMF_ind{j}(k),src_counter) = 0; end; end;
    end;

return;