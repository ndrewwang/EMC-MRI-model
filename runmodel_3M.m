%% MATLAB COMSOL Link
clc; clear; close all
%% COMSOL PARAMETERS
L_param = 0.0080; %[m]
cavg = 3000; %[mol/m^3]
i = 60E-6; %[A]
A = pi*(2E-3)^2; %[m^2]
ipulse = i/A; %[A/m^2]
pos_eval_pts = 100;
current_profile_array = {'0' '7200' '1';  ...
                        '7200' '72000' '0'};
t_eval_string = sprintf('range(0,20,64000)');

%% TRANSPORT & THERMODYNAMIC PROPERTIES
kappa_fn = '(48.93*y^1.5 - 284.8*y^2.5 + 817.7*y^4)^2';
chi_fn = '1 - 18.38*y^0.5 + 155.3*y - 450.6*y^1.5 + 1506*y^2.5';
partfrac_fn = '1.326e-12*c^3 - 1.453e-8*c^2 + 0.0001023*c';
V0_fn = '(1E-3)*(104.105)/(1007.1 - 0.5*-8.1212*c^1.5 - 9*-4.0131E-5*c^10)';
Ve_fn = '(1E-3)*(151.905 - (114.2 + 1.5*-8.1212*c^0.5 + 10*-4.0131E-5*c^9))/(1007.1 - 0.5*-8.1212*c^1.5 - 9*-4.0131E-5*c^10)';
tp0_fn = '0.4107 - 1.487*y + 2.547*y^2';
D_fn = '5.378E-9*y^2 -2.996E-9*y + 4.998E-10';
%% RUN MODEL
clc
output = GalPol_Model(L_param,cavg,ipulse,current_profile_array,t_eval_string,pos_eval_pts,Ve_fn,V0_fn,D_fn,tp0_fn,chi_fn,lambda_fn,kappa_fn,partfrac_fn);
filename = 'comsol_output.csv';
output_table = readtable(filename);