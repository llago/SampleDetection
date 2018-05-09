
function [SDRi,ISRi,SIRi,SARi,permi,cost] = ntf_inst_mu_is ( NMF_CompPerSrcNum, nsrc, stft_win_len, max_iterations, mixture_full_path, source_ref_full_path, dataset_name, experiment_name, init_param)

    %
    % ntf_inst_mu_is
    %
    % Example of usage of MU rules for multichannel NMF decomposition in
    %   linear instantaneous mixture
    %
    %
    % input 
    % -----
    %   - NMF_CompPerSrcNum: number of components per source                 
    %   - nsrc: number of sources to estimate in mixture                     
    %   - stft_win_len: length of STFT window                                
    %   - max_iterations: maximum of iterations for optimation algorithm     
    %   - mixture_full_path: full path to mixture file                       (Optional Argument)
    %   - source_ref_full_path: full path to sources references files        (Optional Argument)
    %
    % output
    % ------
    %
    % estimated source images are written in the results_dir
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Modified by Igor Chame 
    % @ apr 2016
    % igorchame@poli.ufrj.br
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%% Errors %%%
    if nargin < 4, error('Not enough input arguments.'); end;
    if nsrc < 2, error('A mixture has to have at least two sources.'); end;
    if not((NMF_CompPerSrcNum>0)&(mod(NMF_CompPerSrcNum,1)==0)), error('Sources have to be associeted positive Integer number of components.'); end;
    if not((stft_win_len>0)&(mod(stft_win_len,1)==0)), error('The length of the SFTF window has to be a positive Integer.'); end;
    if not((max_iterations>0)&(mod(max_iterations,1)==0)) , error('The number of iterations has to be a positive Integer.'); end;

    % Avoid zeros that may unable fatoration
    eps = 2.2204e-6;   % very small constant

    %%% Initial Setup %%%
    % Setting mixture path
    if nargin < 5 || isempty(mixture_full_path)
        fprintf('Setup missing mixture file path\n');
        fprintf('Please select one\n');
        [mixture_file_name,mixture_path] = uigetfile('*.wav',['Select the Mixture WAV File']);
        mixture_full_path = [mixture_path mixture_file_name];
        fprintf(['Mixture file: ' mixture_full_path '\n']);
    end;
    % Setting path to sources references
    if nargin < 6 || isempty(source_ref_full_path) || size(source_ref_full_path,2) ~= nsrc
        fprintf('Setup missing at least one path to source references\n');
        fprintf(['Please select reference files to ' num2str(nsrc) ' sources\n']);
        for j = 1:nsrc
            [source_ref_file_name{j},source_ref_path{j}] = uigetfile('*.wav',['Select the Source ' int2str(j) ' Reference WAV File']);
            source_ref_full_path{j} = [source_ref_path{j} source_ref_file_name{j}];
        fprintf(['Source ' num2str(j) ' ref file: ' source_ref_full_path{j} '\n']);
        end;
    end;
    % Setting dataset name
    if nargin < 7 || isempty(dataset_name), dataset_name = 'temp'; end;
    % Setting experiment name
    if nargin < 8 || isempty(experiment_name), experiment_name = 'notnamed'; end;
    % Setting initialization parameters
    if nargin < 9 || isempty(init_param) || size(init_param,2) ~= 4, 
        init_rnd = 10; init_iterations = 50; after_iterations = 0;
        repeat_total = 10;
    else
        init_rnd = init_param{1}; init_iterations = init_param{2}; after_iterations = init_param{3};
        repeat_total = init_param{4};
    end;

    %%% Execution Time %%%
    execution_time = datestr(now,'yyyymmdd_HHMMSS');

    %%% Set Paths %%%
    % Generating Path to directories
    implementation_name = mfilename; % implementation name is set by the name of the runner
    results_dir = ['data/' experiment_name '/' dataset_name '-' implementation_name '-' execution_time '/']; % experiment result path is relative to this script path
    % Checking if experiment_path exists
    if (exist(results_dir) ~= 7) mkdir(results_dir); end% if it does not, then creat it
    % Set sources directories
    addpath('src_aux_tools'); % auxiliar tools
    addpath('src_aux_tools/log4m'); % logging
    addpath('src_methods'); % methods for source separation

    %%% Setup LOG %%
    LOG = 0;
    LOG = log4m.getLogger([results_dir 'logfile.txt']);
    LOG.setCommandWindowLevel(LOG.ALL);
    LOG.setLogLevel(LOG.INFO);
    % Write log header
    LOG.info('experiment',['% Experiment: ' experiment_name ' %']);
    LOG.info('experiment',['% Implementation: ' implementation_name ' %']);
    LOG.info('experiment',[' Results dir: ' results_dir]);
    LOG.info('experiment',[' Mixture: ' mixture_full_path]);
    for j = 1:nsrc, LOG.info('experiment',[' Source ' num2str(j) '/' num2str(nsrc) ': ' source_ref_full_path{j}]); end;
    LOG.info('experiment',['% Setup %']);
    LOG.info('experiment',['  nsrc: ' num2str(nsrc) ', NMF_CompPerSrcNum: ' num2str(NMF_CompPerSrcNum) ', max_iterations: ' num2str(max_iterations)]);

    %%% Input time-frequency representation %%%
    LOG.info('experiment',['% Mixture %']);
    % fprintf('Input time-frequency representation\n');
    [x, fs]=audioread(mixture_full_path);
    
    nchann = size(x,2);
    % very small normal noise to avoid zeros
    for chan=1:nchann, if min(abs(x(:,chan))) == 0, eps_noise = randn(size(x(:,chan),1),1); x(:,chan) = x(:,chan) + eps*(eps_noise ./ max(abs(eps_noise))); end; end;
    % Guarantees normalization of WAV file
    % for chann = 1:nchann, x(:,chann) = x(:,chann) / max(abs(x(:,chann))); end;
    x = x.';
    mix_nsamp = size(x,2);
    LOG.info('experiment',['  fs: ' num2str(fs,'%8.0u')  ', mix_nsamp: ' num2str(mix_nsamp,'%8.0u') ' (' num2str(mix_nsamp/fs,'%6.1f') 's)' ', mix_nchann: ' num2str(nchann,'%2.0u') '' ]);

    LOG.info('experiment',['% Time-frequency representation %']);
    X=stft_multi(x,stft_win_len);
    nbin = size(X,1);
    nfram = size(X,2);
    LOG.info('experiment',['  nbin: ' num2str(nbin,'%6.0u') ', nfram: ' num2str(nfram,'%6.0u') ', stft_win_len: ' num2str(stft_win_len,'%6.0u') ' (' num2str(round(1000*stft_win_len/fs),'%4.0u') 'ms)' ]);


    LOG.info('experiment','% Initialization %');
    K = NMF_CompPerSrcNum * nsrc;
    source_NMF_ind = cell(1,nsrc);
    V = abs(X).^2; % Power density of mixture STFT
    LOG.info('experiment',['  init_rnd: ' num2str(init_rnd,'%3.0u') ', init_iterations: ' num2str(init_iterations,'%3.0u') ', after_iterations: ' num2str(after_iterations,'%3.0u')]);
    
    % Run Init and Optimization repeat_total times, and select with minimum cost
    cost = zeros(repeat_total,max_iterations); cost_best = 0;
    for repeat=1:repeat_total,
        LOG.info('experiment',['% Optimization Problem (' num2str(repeat, '%2.0u') '/' num2str(repeat_total, '%2.0u') ') %']);
        tic;
        %%% Random initialization of multichannel NMF parameters %%%%
        % Initialize W and H with random variables after Monte Carlo
        [Q_init, W_init, H_init] = init_ntf_mu_is_pow(X, K, init_rnd, init_iterations, after_iterations);
        elapsedTime = toc; LOG.info('experiment',['  Time elapsed Initialization: ' num2str(elapsedTime,'%9.3f') 's']);

        %%% Optimization Problem %%%
        % run max_iterations iterations of multichannel NMF MU rules 
        tic;
        [Q, W, H, cost(repeat,:)] = ntf_mu_is(V, max_iterations, Q_init, W_init, H_init);
        elapsedTime = toc;  LOG.info('experiment',['  Time elapsed Optimization: ' num2str(elapsedTime,'%9.3f') 's']);
        LOG.info('experiment',['  Min cost: ' num2str(cost(repeat,end),'%9.2f')]);
    
        if ( cost_best == 0 ) || ( cost(repeat,end) < cost_best(1,end) )
            LOG.info('experiment','  !New best!');
            Q_best = Q;
            W_best = W;
            H_best = H;
            cost_best = cost(repeat,end);
        end;
    end;
    Q = Q_best;
    W = W_best;
    H = H_best;
    clear Q_init W_init H_init Q_init W_best H_best;
    h = figure; semilogy(1:max_iterations,cost); title('Divergencia de Itakura-Saito'); ylabel('$C(\cdot)$','Interpreter','latex'); xlabel('Iterations'); grid on;
    saveas( h, [results_dir 'cost.fig'], 'fig');
    saveas( h, [results_dir 'cost.png'], 'png');
    h = figure; loglog(1:max_iterations,cost); ylabel('$C(\cdot)$','Interpreter','latex'); xlabel('$s_{fat}$','Interpreter','latex'); grid on;
    saveas( h, [results_dir 'cost.eps'], 'epsc');



    %%% Reading sources references%%%
    % to evaluate estimations and to associate components
    s_ref=zeros(nsrc,mix_nsamp,nchann);
    S_ref_STFT_power = zeros(nsrc,nbin,nfram,nchann);
    for j=1:nsrc,
        s_ref_temp = audioread(source_ref_full_path{j});
        s_ref(j,:,:)=reshape(s_ref_temp(1:mix_nsamp,1:nchann),1,mix_nsamp,nchann);
        
        S_ref_STFT = stft_multi(s_ref_temp(1:mix_nsamp,1:nchann).',stft_win_len);
        S_ref_STFT_power(j,:,:,:) = reshape( (abs(S_ref_STFT).^2), 1, nbin, nfram, nchann);
    end;
    clear s_ref_temp S_ref_STFT; % clear temp variables

    %%% Associacao de componentes com fontes %%%
    LOG.info('experiment','% Components Association with Sources %');
    tic;
    [source_NMF_ind] = assoc_compToSrc_multichan_ideal(Q,W,H,S_ref_STFT_power);
    elapsedTime = toc;  LOG.info('experiment',['  Time elapsed Association: ' num2str(elapsedTime,'%9.3f') 's']);

    %%% Reconstruction of the spatial source %%%
    LOG.info('experiment','% Synthesis %');
    % fprintf('Reconstruction of the spatial source components\n');
    % for k =1:K, component_NMF_ind{k} = k; end;
    tic;
    S_est_STFT = ntf_recons_im(X,Q,W,H,source_NMF_ind);
    s_est = istft_multi(S_est_STFT,mix_nsamp);
    elapsedTime = toc;  LOG.info('experiment',['  Time elapsed Synthesis: ' num2str(elapsedTime,'%9.3f') 's']);

    %%% Evaluate Separation %%%
    LOG.info('experiment','% Evaluate Separation %');
    tic;
    [SDRi,ISRi,SIRi,SARi,permi]=bss_eval_images(s_est,s_ref);
    elapsedTime = toc;  LOG.info('experiment',['  Time elapsed Evaluation: ' num2str(elapsedTime,'%9.3f') 's']);
    LOG.info('experiment','[  SRC, SDRi, ISRi, SIRi, SARi]');
    for j=1:nsrc, LOG.info('experiment',[ '[ ' num2str(j,'%02.0u') ' , ' num2str(SDRi(j),'%+5.1f') ' , ' num2str(ISRi(j),'%+5.1f') ' , ' num2str(SIRi(j),'%+5.1f') ' , ' num2str(SARi(j),'%+5.1f') ' ]']); end;

    %%% Save WAV file of source estimations %%%
    % Save estimation of highest SDR source
    for j=1:nsrc,
        % Guarantees normalization of channels
        % for chann = 1:size(s_est,3), s_est(j,:,chann) = s_est(j,:,chann) ./ max(abs(s_est(j,:,chann))); end;
        audiowrite([results_dir 'src_est_' num2str(j,'%02.0u') '.wav'],reshape(s_est(permi(j),:,:),mix_nsamp,nchann),fs); 
    end;

    %%% Plot estimated W and H %%%
    LOG.info('experiment','% Plot estimated W and H %');
    fprintf('Plot estimated W and H\n');
    h = figure;
    plot_ind = 1;
    for k = 1:NMF_CompPerSrcNum
        for j = 1:nsrc
            subplot(NMF_CompPerSrcNum, nsrc, plot_ind);
            plot(max(log10(max(W(:,source_NMF_ind{j}(k)), 1e-40)), -10));
            title(sprintf('Source_%d, log10(W_%d)', j, k));
            plot_ind = plot_ind + 1;
        end;
    end;
    saveas( h, [results_dir 'w.fig'], 'fig');
    saveas( h, [results_dir 'w.png'], 'png');
    h = figure;
    plot_ind = 1;
    for k = 1:NMF_CompPerSrcNum
        for j = 1:nsrc
            subplot(NMF_CompPerSrcNum, nsrc, plot_ind);
            plot(H(source_NMF_ind{j}(k),:));
            title(sprintf('Source_%d, H_%d', j, k));
            plot_ind = plot_ind + 1;
        end;
    end;
    saveas( h, [results_dir 'h.fig'], 'fig');
    saveas( h, [results_dir 'h.png'], 'png');

    %%% Saving Workspace %%%
    LOG.info('experiment','% Saving workspace %');
    % save selected variables into a .MAT file
    %save example.mat A B -v7.3;
    clear X x s_ref S_ref_STFT_power LOG V chann plot_ind j k f h; % clear some information to reduce file size
    save ([results_dir 'workspace.mat']);

    % Rename result_dir to include mean of SIR for all sources
    % movefile(results_dir,[results_dir(1:end-1) '_' num2str(mean(SIRi), '%+5.1f')]);
end
