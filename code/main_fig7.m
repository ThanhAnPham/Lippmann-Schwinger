%% Produce figures 7 from paper

% This script reproduces the figure 7 from [1]
%
% Dataset from the Fresnel Institute [2]
% Reference
% [1] E. Soubies, T.-A. Pham, and M. Unser,
%     “Efficient inversion of multiple-scattering model for optical diffraction tomography,”
%     Opt. Express 25, 21786–21800 (2017).
% [2] Jean-Michel Geffrin, Pierre Sabouroux and Christelle Eyraud, “Free space experimental scattering database
%     continuation: experimental set-up and measurement precision,” Inverse Probl. 21, (2005).

clear
addpath('../');
set(0,'DefaultTextInterpreter','LaTex');
%% Figure 7


load('FoamDielExtTM_3_GHz_conj_1.mat');
main_param_fig_7;
main_lippmann;