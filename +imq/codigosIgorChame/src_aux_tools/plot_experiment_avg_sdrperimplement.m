function plot_experiment_avg_sdrperimplement( experiment_name, dataset_name, implementation_name, NMF_CompPerSrcNum, sdr)
  
    %
    %
    % input 
    % -----
    %
    %
    % output
    % ------
    %
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Modified by Igor Chame 
    % @ apr 2016
    % igorchame@poli.ufrj.br
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    results_dir = ['data/' experiment_name '/']; % experiment result path is relative to this script path
    filename = ['rslt_' dataset_name '_param' num2str(NMF_CompPerSrcNum,'%2.0u')];

    h = figure; 
    barh(mean(sdr,2)); 
    hold on  
    ylabel('Implementation'); 
    xlabel('SDR (dBs)'); 
    set(gca, 'YTickLabel',strrep(implementation_name, '_', ' '))
    grid on;

    for implement=1:size(implementation_name,2), max_sdr(implement) = max(sdr(implement,:)); min_sdr(implement) = min(sdr(implement,:)); end;
    plot(max_sdr',(1:size(implementation_name,2)),'--','Color','black'); 
    plot(min_sdr',(1:size(implementation_name,2)),'Color','black'); 
    hold off
    saveas( h, [results_dir filename '.eps'], 'epsc');
    saveas( h, [results_dir filename '.fig'], 'fig');
    %saveas( h, [results_dir filename '.png'], 'png');

    close all
end 
