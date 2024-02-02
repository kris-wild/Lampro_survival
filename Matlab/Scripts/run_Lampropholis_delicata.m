close all; 
global pets

pets = {'Lampropholis_delicata'};
check_my_pet(pets); 


% Read about how to set estimation and output options (estim_options) on the online
% manual: http://www.debtheory.org/wiki/index.php?title=Run_file

estim_options('default'); 
estim_options('max_step_number',5e3); 
estim_options('max_fun_evals',5e3);  

estim_options('pars_init_method', 2);
estim_options('results_output', -3);
estim_options('method', 'nm');

estim_pars; 
load('results_Lampropholis_delicata.mat');
[stat, txtStat] = statistics_st(metaPar.model, par);
printstat_st(stat, txtStat)