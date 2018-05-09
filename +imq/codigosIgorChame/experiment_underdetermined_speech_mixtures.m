function experiment_underdetermined_speech_mixtures ()

    %
    % input 
    % -----
    %
    % ...
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

    %%% Experiment Time %%%
    experiment_time = datestr(now,'yyyymmdd_HHMMSS');

    %% Experiment name
    experiment_name = [mfilename '-' experiment_time ]; % experiment name is set by the name of the runner

    %% Naming Datasets
    dataset_name{1} = 'SiSEC2008_female3';
    dataset_name{2} = 'SiSEC2008_male3';
    dataset_name{3} = 'SiSEC2010_female4';
    dataset_name{4} = 'SiSEC2010_male4';
    dataset_name{5} = 'SiSEC2010_femalemale_wChannDiff';
    dataset_name{6} = 'SiSEC2010_femalemale_noChannDiff';

    %% Path to Bases
    % depends on the computer running experiment
    % basic_path = '/Volumes/Mac HD/igorchame/Google Drive/Documents/1.Faculdade/UFRJ/2015.1_projetoFinal/';
    basic_path = '/home/igor.chame/Documents/';
    
    % Setting mixture path
    mixture_full_path{1} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2008_dev/dev1_female3_inst_mix.wav' ];
    mixture_full_path{2} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2008_dev/dev1_male3_inst_mix.wav' ];
    mixture_full_path{3} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2010_dev/dev2_female4_inst_mix.wav' ];
    mixture_full_path{4} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2010_dev/dev2_male4_inst_mix.wav' ];
    mixture_full_path{5} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2010_dev/SiSEC2010_femalemale_1_wChannDiff_mix.wav' ];
    mixture_full_path{6} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2010_dev/SiSEC2010_femalemale_1_noChannDiff_mix.wav' ];


    % Setting path to sources references
    source_ref_full_path{1}{1} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2008_dev/dev1_female3_inst_sim_1.wav' ]; 
    source_ref_full_path{1}{2} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2008_dev/dev1_female3_inst_sim_2.wav' ]; 
    source_ref_full_path{1}{3} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2008_dev/dev1_female3_inst_sim_3.wav' ];
    
    source_ref_full_path{2}{1} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2008_dev/dev1_male3_inst_sim_1.wav' ]; 
    source_ref_full_path{2}{2} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2008_dev/dev1_male3_inst_sim_2.wav' ]; 
    source_ref_full_path{2}{3} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2008_dev/dev1_male3_inst_sim_3.wav' ];

    source_ref_full_path{3}{1} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2010_dev/dev2_female4_inst_sim_1.wav' ]; 
    source_ref_full_path{3}{2} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2010_dev/dev2_female4_inst_sim_2.wav' ]; 
    source_ref_full_path{3}{3} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2010_dev/dev2_female4_inst_sim_3.wav' ];
    source_ref_full_path{3}{4} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2010_dev/dev2_female4_inst_sim_4.wav' ];
    
    source_ref_full_path{4}{1} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2010_dev/dev2_male4_inst_sim_1.wav' ]; 
    source_ref_full_path{4}{2} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2010_dev/dev2_male4_inst_sim_2.wav' ]; 
    source_ref_full_path{4}{3} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2010_dev/dev2_male4_inst_sim_3.wav' ]; 
    source_ref_full_path{4}{4} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2010_dev/dev2_male4_inst_sim_4.wav' ];

    source_ref_full_path{5}{1} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2010_dev/SiSEC2010_femalemale_1_wChannDiff_sim_1.wav' ]; 
    source_ref_full_path{5}{2} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2010_dev/SiSEC2010_femalemale_1_wChannDiff_sim_2.wav' ]; 
    
    source_ref_full_path{6}{1} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2010_dev/SiSEC2010_femalemale_1_noChannDiff_sim_1.wav' ]; 
    source_ref_full_path{6}{2} = [ basic_path 'bases/SiSEC_2013_Tasks/1.Underdetermined_speech_and_music_mixtures/SiSEC2010_dev/SiSEC2010_femalemale_1_noChannDiff_sim_2.wav' ];
    
    
    %% Getting how many sources are present by the number of references
    % nsrc = size(source_ref_full_path,2); 
    nsrc(1) = 3;
    nsrc(2) = 3;
    nsrc(3) = 4;
    nsrc(4) = 4;
    nsrc(5) = 2;
    nsrc(6) = 2;

    %% Setting number of components per source
    NMF_CompPerSrcNum{1} = [ 10];
    NMF_CompPerSrcNum{2} = [ 10];
    NMF_CompPerSrcNum{3} = [ 10];
    NMF_CompPerSrcNum{4} = [ 10];
    NMF_CompPerSrcNum{5} = [ 10];
    NMF_CompPerSrcNum{6} = [ 10];

    %% Setting number of max iterations for optimization
    max_iterations = 5000;

    %% Setting STFT window length to be a power of two
    stft_win_len = 1024;

    %% Initialization parameters
    init_rnd = 10; init_iterations = 20; after_iterations = 0; repeat_total = 10;
    init_param = cell(1,4); init_param{1} = init_rnd; init_param{2} = init_iterations; init_param{3} = after_iterations; init_param{4} = repeat_total;

    %% Different parameters
    run_param_max = size(NMF_CompPerSrcNum{1},2);

    %% Different implementations that will be executed
    diff_implement = 8;
    implementation_name = cell(1,diff_implement);
    implementation_name{1} = 'nmf_inst_mu_is';
    implementation_name{2} = 'nmf_inst_mu_kl';
    implementation_name{3} = 'nmf_inst_em_is';
    implementation_name{4} = 'ntf_inst_mu_is';
    implementation_name{5} = 'ntf_inst_mu_kl';
    implementation_name{6} = 'ntf_inst_em_is';
    implementation_name{7} = 'ntfclus_inst_mu_is';
    implementation_name{8} = 'ntfclus_inst_mu_kl';
    addpath('src_implementations'); % implementations

    %% Result tables for one parameters group
    SDR     = zeros( run_param_max, size(dataset_name,2), diff_implement, max(nsrc));
    ISR     = zeros( run_param_max, size(dataset_name,2), diff_implement, max(nsrc));
    SIR     = zeros( run_param_max, size(dataset_name,2), diff_implement, max(nsrc));
    SAR     = zeros( run_param_max, size(dataset_name,2), diff_implement, max(nsrc));
    perm    = zeros( run_param_max, size(dataset_name,2), diff_implement, max(nsrc));
    cost    = zeros( run_param_max, size(dataset_name,2), diff_implement, repeat_total, max(max_iterations));

    %%% Executing experiments in datasets and with parameters %%%
    for run_param=1:run_param_max,
        for dataset=1:size(dataset_name,2),

            %% nmf_inst_mu_is
            implement = 1; 
            [SDR(run_param,dataset,implement,1:nsrc(dataset)),ISR(run_param,dataset,implement,1:nsrc(dataset)),SIR(run_param,dataset,implement,1:nsrc(dataset)),SAR(run_param,dataset,implement,1:nsrc(dataset)),perm(run_param,dataset,implement,1:nsrc(dataset)),cost(run_param,dataset,implement,:,:)] = ...
                nmf_inst_mu_is( NMF_CompPerSrcNum{dataset}(run_param), nsrc(dataset), stft_win_len, max_iterations, mixture_full_path{dataset}, source_ref_full_path{dataset}, dataset_name{dataset}, experiment_name, init_param);
            close all; 

            %% nmf_inst_mu_kl
            implement = 2; 
            [SDR(run_param,dataset,implement,1:nsrc(dataset)),ISR(run_param,dataset,implement,1:nsrc(dataset)),SIR(run_param,dataset,implement,1:nsrc(dataset)),SAR(run_param,dataset,implement,1:nsrc(dataset)),perm(run_param,dataset,implement,1:nsrc(dataset)),cost(run_param,dataset,implement,:,:)] = ...
                nmf_inst_mu_kl( NMF_CompPerSrcNum{dataset}(run_param), nsrc(dataset), stft_win_len, max_iterations, mixture_full_path{dataset}, source_ref_full_path{dataset}, dataset_name{dataset}, experiment_name, init_param);
            close all; 

            %% nmf_inst_em_is
            implement = 3; 
            [SDR(run_param,dataset,implement,1:nsrc(dataset)),ISR(run_param,dataset,implement,1:nsrc(dataset)),SIR(run_param,dataset,implement,1:nsrc(dataset)),SAR(run_param,dataset,implement,1:nsrc(dataset)),perm(run_param,dataset,implement,1:nsrc(dataset)),cost(run_param,dataset,implement,:,:)] = ...
                nmf_inst_em_is( NMF_CompPerSrcNum{dataset}(run_param), nsrc(dataset), stft_win_len, max_iterations, mixture_full_path{dataset}, source_ref_full_path{dataset}, dataset_name{dataset}, experiment_name, init_param);
            close all; 

            %% ntf_inst_mu_is
            implement = 4; 
            [SDR(run_param,dataset,implement,1:nsrc(dataset)),ISR(run_param,dataset,implement,1:nsrc(dataset)),SIR(run_param,dataset,implement,1:nsrc(dataset)),SAR(run_param,dataset,implement,1:nsrc(dataset)),perm(run_param,dataset,implement,1:nsrc(dataset)),cost(run_param,dataset,implement,:,:)] = ...
               ntf_inst_mu_is( NMF_CompPerSrcNum{dataset}(run_param), nsrc(dataset), stft_win_len, max_iterations, mixture_full_path{dataset}, source_ref_full_path{dataset}, dataset_name{dataset}, experiment_name, init_param);
            close all; 

            %% ntf_inst_mu_kl
            implement = 5; 
            [SDR(run_param,dataset,implement,1:nsrc(dataset)),ISR(run_param,dataset,implement,1:nsrc(dataset)),SIR(run_param,dataset,implement,1:nsrc(dataset)),SAR(run_param,dataset,implement,1:nsrc(dataset)),perm(run_param,dataset,implement,1:nsrc(dataset)),cost(run_param,dataset,implement,:,:)] = ...
                ntf_inst_mu_kl( NMF_CompPerSrcNum{dataset}(run_param), nsrc(dataset), stft_win_len, max_iterations, mixture_full_path{dataset}, source_ref_full_path{dataset}, dataset_name{dataset}, experiment_name, init_param);
            close all; 

            %% ntf_inst_em_is
            implement = 6; 
            [SDR(run_param,dataset,implement,1:nsrc(dataset)),ISR(run_param,dataset,implement,1:nsrc(dataset)),SIR(run_param,dataset,implement,1:nsrc(dataset)),SAR(run_param,dataset,implement,1:nsrc(dataset)),perm(run_param,dataset,implement,1:nsrc(dataset)),cost(run_param,dataset,implement,:,:)] = ...
                ntf_inst_em_is( NMF_CompPerSrcNum{dataset}(run_param), nsrc(dataset), stft_win_len, max_iterations, mixture_full_path{dataset}, source_ref_full_path{dataset}, dataset_name{dataset}, experiment_name, init_param);
            close all; 

            %% ntfclus_inst_mu_is
            implement = 7; 
            [SDR(run_param,dataset,implement,1:nsrc(dataset)),ISR(run_param,dataset,implement,1:nsrc(dataset)),SIR(run_param,dataset,implement,1:nsrc(dataset)),SAR(run_param,dataset,implement,1:nsrc(dataset)),perm(run_param,dataset,implement,1:nsrc(dataset)),cost(run_param,dataset,implement,:,:)] = ...
                ntfclus_inst_mu_is( NMF_CompPerSrcNum{dataset}(run_param), nsrc(dataset), stft_win_len, max_iterations, mixture_full_path{dataset}, source_ref_full_path{dataset}, dataset_name{dataset}, experiment_name, init_param);
            close all; 

            %% ntfclus_inst_mu_kl
            implement = 8; 
            [SDR(run_param,dataset,implement,1:nsrc(dataset)),ISR(run_param,dataset,implement,1:nsrc(dataset)),SIR(run_param,dataset,implement,1:nsrc(dataset)),SAR(run_param,dataset,implement,1:nsrc(dataset)),perm(run_param,dataset,implement,1:nsrc(dataset)),cost(run_param,dataset,implement,:,:)] = ...
                ntfclus_inst_mu_kl( NMF_CompPerSrcNum{dataset}(run_param), nsrc(dataset), stft_win_len, max_iterations, mixture_full_path{dataset}, source_ref_full_path{dataset}, dataset_name{dataset}, experiment_name, init_param);
            close all; 
        end;
        size_sdr = size(SDR(run_param,:,:,:)); create_experiment_latex_table(experiment_name, NMF_CompPerSrcNum{dataset}(run_param), dataset_name, implementation_name, reshape(SDR(run_param,:,:,:),[size_sdr(2:end) 1]), reshape(ISR(run_param,:,:,:),[size_sdr(2:end) 1]), reshape(SIR(run_param,:,:,:),[size_sdr(2:end) 1]), reshape(SAR(run_param,:,:,:),[size_sdr(2:end) 1]), reshape(perm(run_param,:,:,:),[size_sdr(2:end) 1]), nsrc);
    end;

    %% Ploting max and min costs per number of component
    for dataset=1:size(dataset_name,2),
        % for implement=1:diff_implement,
        %     plot_experiment_cost_maxmin( experiment_name, dataset_name{dataset}, implementation_name{implement}, NMF_CompPerSrcNum{dataset}, squeeze(cost(:,dataset,implement,:,end)));
        %     plot_experiment_sdr_perrun( experiment_name, dataset_name{dataset}, implementation_name{implement}, NMF_CompPerSrcNum{dataset}, squeeze(SDR(:,dataset,implement,:)));
        % end;
        for run_param=1:run_param_max,
            plot_experiment_avg_sdrperimplement(experiment_name, dataset_name{dataset}, implementation_name, NMF_CompPerSrcNum{dataset}(run_param), squeeze(SDR(run_param,dataset,:,:)));
        end;
    end;

end

