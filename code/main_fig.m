%% Produce figures 6 & 7 from paper
clear
%% Figure 6

%Load Shepp
load('shepp_logan.mat');

main_param_fig_6;

%main_lippmann;
main_seagle;

%% Figure 7
clear uem;
load('FoamDielExtTM_3_GHz.mat');
main_param_fig_7;
main_lippmann;
%main_seagle;