function create_experiment_latex_table (experiment_name, run_param, dataset_name, implementation_name, SDR, ISR, SIR, SAR, perm, nsrc_data)

    %
    % Create a Latex table with the results of 
    % the experiment in multiple datasets
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

    for dataset=1:size(dataset_name,2)

        filename = [ 'rslt_' dataset_name{dataset} '_param' num2str(run_param,'%2.0u') '.txt'];
        nsrc = nsrc_data(dataset);
        file = fopen([results_dir  filename],'w');

        %%% Top %%%
        fprintf(file,'\\begin{table}\n');
        fprintf(file,'    \\centering\n');
        fprintf(file,'    \\begin{tabular}{'); for src=1:nsrc, fprintf(file,' | c '); end; fprintf(file,'| }\n');
        fprintf(file,'        \\hline\n');

        %%% Header %%%
        fprintf(file,'         \\multicolumn{%u}{|c|}{\\textbf{%s}} \\\\ \\hline \n',nsrc+1,strrep(dataset_name{dataset}, '_', '\_'));
        fprintf(file,'         ');  for j =1:nsrc, fprintf(file,'& $s_{%1.0u}$ ',j); end; fprintf(file,'\\\\  \n');
        fprintf(file,'         '); for j =1:nsrc, fprintf(file,'& Inst%1.0u ',j); end; fprintf(file,'\\\\ \\hline \n');
        
        for implementation=1:size(implementation_name,2), 
            fprintf(file,'          & \\multicolumn{%u}{c|}{%s} \\\\ \\hline \n',nsrc,strrep(implementation_name{implementation}, '_', '\_'));
            
            fprintf(file,'         SDR '); for j =1:nsrc, fprintf(file,'& %8.1f ',SDR(dataset,implementation,j))); end; fprintf(file,'\\\\ \\hline \n');
            fprintf(file,'         SIR '); for j =1:nsrc, fprintf(file,'& %8.1f ',SIR(dataset,implementation,j))); end; fprintf(file,'\\\\ \\hline \n');
            fprintf(file,'         SAR '); for j =1:nsrc, fprintf(file,'& %8.1f ',SAR(dataset,implementation,j))); end; fprintf(file,'\\\\ \\hline \n');   
        end;

        %%% Footer %%%
        fprintf(file,'    \\end{tabular} \n');
        fprintf(file,'    \\caption{Title} \n');
        fprintf(file,'    \\label{tbl:label} \n');
        fprintf(file,'\\end{table} \n');

        fclose(file);
    end
end