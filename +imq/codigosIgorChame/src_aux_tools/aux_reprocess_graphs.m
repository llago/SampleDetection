% Select first estimate
yd = get(gco,'YData')';
xd = get(gco,'XData');

% Select second estimate
yd(:,2) = get(gco,'YData')';

experiment_name = '1';
dataset_name = 'SiSEC2008_tamy-que_pena_tanto_faz';
implementation_name = 'nmf_inst_mu_kl';

plot_experiment_sdr_perrun (experiment_name, dataset_name, implementation_name, xd, yd);
