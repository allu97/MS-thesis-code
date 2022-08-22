clearvars;clc;
load('QuarriesREE.mat')
load('SamplesREE.mat')
%%
% Percent of NaN values in quarries REE data
XC_nan_p = sum(isnan(XC_ree),'all')/numel(XC_ree)
% Percent of NaN values in samples REE data
X_nan_p = sum(isnan(X_ree),'all')/numel(X_ree)
%% samples vs quarries REE
idx_XC = isnan(XC_ree);
idx_X = isnan(X_ree);
XC_ree(idx_XC) = 0;
X_ree(idx_X) = 0;

%%
% Check if ANY value in column contains NaN
% Samples
cols_any_NaN_samples = any(isnan(X_ree),1);

var_any_samples_after = X_ree_var(~cols_any_NaN_samples);

% Quarries
cols_any_NaN_quarries = any(isnan(XC_ree),1);

var_any_quarries_after = X_ree_var(~cols_any_NaN_quarries);

% Check if ALL values in column contains NaN
% Samples
cols_all_NaN_samples = all(isnan(X_ree),1);

var_all_samples_after = X_ree_var(~cols_all_NaN_samples);

% Quarries
% Heponiemi data is cleaned from NaN value (replaced with 0),
% so we are only interested in Hailniemi and Hämeenkylä quarries
cols_all_NaN_quarries = all(isnan(XC_ree(52:end,:)),1);

var_all_quarries_after = X_ree_var(~cols_all_NaN_quarries);

% Indices from quarries
new1_X_ree = X_ree(:,~cols_all_NaN_quarries);
new1_XC_ree = XC_ree(:,~cols_all_NaN_quarries);

% Indices from samples
new2_X_ree = X_ree(:,~cols_all_NaN_samples);
new2_XC_ree = XC_ree(:,~cols_all_NaN_samples);

% Combine both indices
cols_all_NaN_both = cols_all_NaN_quarries | cols_all_NaN_samples;
new3_X_ree = X_ree(:,~cols_all_NaN_both);
new3_XC_ree = XC_ree(:,~cols_all_NaN_both);
var_all_both_after = X_ree_var(~cols_all_NaN_both);

vars_removed_1 = X_ree_var(cols_all_NaN_both);
% Sample indices remove more variables than quarries
% These are removed from both data
% 'Mg' 'Ti' 'V' 'Cr' 'Mn' 'Se' 'Ag' 'Cd' 'Sb' 'Pr' 'Nd' 'W' 'Hg' 'Bi'

% Next let's remove variables with only few samples.
% Lets say if the variable has more than 50% of NaN values, remove it

% Number of observations in samples
n = size(new3_X_ree,1);
% Number of observations in quarries (no Heponiemi)
m = size(new3_XC_ree(52:end,:),1);

nans_new3_X = isnan(new3_X_ree);
nans_new3_XC = isnan(new3_XC_ree(52:end,:));

cols_X_over_50 = sum(nans_new3_X)/n > 0.50;
cols_XC_over_50 = sum(nans_new3_XC)/m > 0.50;

combine_50 = cols_X_over_50 | cols_XC_over_50;

new4_X_ree = new3_X_ree(:,~combine_50);
new4_XC_ree = new3_XC_ree(:,~combine_50);

vars_removed_2 = var_all_both_after(combine_50);
var_over_50 = var_all_both_after(~combine_50);
% These have over 50% of NaN values
% 'P' 'S' 'Co' 'Ni' 'Cu' 'Zn' 'As' 'Sn' 'La' 'Ta' 'Th' 'U'
%%
threshold = 0.4;
load('Hailniemi.mat','-regexp','[^j]$')
% Heponiemi REE
[hn_ree_new_data, hn_ree_new_vars, hn_ree_vars_removed,idx1_ree] = ...
    lessNaNs(hn_ree_data,hn_ree_var,threshold); 
% Heponiemi no REE
[hn_new_data, hn_new_vars, hn_vars_removed,idx1] = ...
    lessNaNs(hn_data,hn_var,threshold); 

load('Hameenkyla.mat','-regexp','[^j]$')
[hk_ree_new_data, hk_ree_new_vars, hk_ree_vars_removed,idx2_ree] = ...
    lessNaNs(hk_ree_data,hk_ree_var,threshold); 
[hk_new_data, hk_new_vars, hk_vars_removed,idx2] = ...
    lessNaNs(hk_data,hk_var,threshold); 

load('Sample1.mat','-regexp','[^j]$')
[s1_ree_new_data, s1_ree_new_vars, s1_ree_vars_removed,idx3_ree] = ...
    lessNaNs(s1_ree_data,s1_ree_var,threshold); 
[s1_new_data, s1_new_vars, s1_vars_removed,idx3] = ...
    lessNaNs(s1_data,s1_var,threshold); 

load('Sample2.mat','-regexp','[^j]$')
[s2_ree_new_data, s2_ree_new_vars, s2_ree_vars_removed,idx4_ree] = ...
    lessNaNs(s2_ree_data,s2_ree_var,threshold); 
[s2_new_data, s2_new_vars, s2_vars_removed,idx4] = ...
    lessNaNs(s2_data,s2_var,threshold); 

idx_ree = idx1_ree | idx2_ree | idx3_ree | idx4_ree;
idx = idx1 | idx2 | idx3 | idx4;

function [newData, newVars, vars_removed, idx_removed] = lessNaNs(data,vars,threshold)
% Check if all observations in column have NaN
% idx_removed = all(isnan(data),1);
% 
% newVars = vars(~idx_removed);
% newData = data(:,~idx_removed);

% Next let's remove variables with only few samples.
% Lets say if the variable has more than <threshold> of NaN values,
% remove it.

% Number of observations data
n = size(data,1);

nans_data = isnan(data);

idx_removed = sum(nans_data)/n > threshold;

newData = data(:,~idx_removed);

vars_removed = vars(idx_removed);
newVars = vars(~idx_removed);

end