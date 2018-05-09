function plot_experiment_cost_maxmin( experiment_name, dataset_name, implementation_name, NMF_CompPerSrcNum, cost)
  
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
    filename = [dataset_name '-' implementation_name '-cost' ];

    h = figure; 
    hold on  
    plot(NMF_CompPerSrcNum,cost,'*'); 
    ylabel('$C(\cdot)$','Interpreter','latex'); 
    xlabel('$M_{source}$','Interpreter','latex'); 
    grid on;

    for param=1:size(NMF_CompPerSrcNum,2), max_cost(param) = max(cost(param,:)); min_cost(param) = min(cost(param,:)); end;
    plot(NMF_CompPerSrcNum,max_cost,'--','Color','black'); 
    plot(NMF_CompPerSrcNum,min_cost,'Color','black'); 
    hold off
    saveas( h, [results_dir filename '.eps'], 'epsc');
    saveas( h, [results_dir filename '.fig'], 'fig');
    %saveas( h, [results_dir filename '.png'], 'png');

    close all
end 
