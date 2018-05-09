function plot_experiment_sdr_perrun( experiment_name, dataset_name, implementation_name, NMF_CompPerSrcNum, sdr)
  
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
    filename = [dataset_name '-' implementation_name '-sdr'];

    % Remove zero columns
    sdr( :, all(~sdr,1) ) = [];

    h = figure; 
    hold on  
    bar(NMF_CompPerSrcNum,mean(sdr,2)); 
    plot(NMF_CompPerSrcNum,sdr,'*'); 
    ylabel('SDR (dBs)'); 
    xlabel('$M_{source}$','Interpreter','latex'); 
    grid on;
    
    ax = gca;
    ax.YAxisLocation = 'origin';

    for param=1:size(NMF_CompPerSrcNum,2), max_sdr(param) = max(sdr(param,:)); min_sdr(param) = min(sdr(param,:)); end;
    plot(NMF_CompPerSrcNum,max_sdr,'--','Color','black'); 
    plot(NMF_CompPerSrcNum,min_sdr,'Color','black'); 
    hold off
    saveas( h, [results_dir filename '.eps'], 'epsc');
    saveas( h, [results_dir filename '.fig'], 'fig');
    %saveas( h, [results_dir filename '.png'], 'png');

    close all
end 
