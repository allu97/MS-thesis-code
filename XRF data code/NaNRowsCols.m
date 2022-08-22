clearvars;clc
%% Hailniemi
load('Hailniemi.mat')
[hn_data1, hn_var1, hn_removed1, hn_idx_removed1] = ...
    allNaNs(hn_data,hn_var);

[hn_ree_data1, hn_ree_var1, hn_ree_removed1, hn_ree_idx_removed1] = ...
    allNaNs(hn_ree_data,hn_ree_var);

Hn_table = createNaNTable(hn_data,hn_var,hn_obj);

[hn_data2, hn_var2, hn_removed2, hn_idx_removed2] = ...
    lessNaNs(hn_data,hn_var,0.89);

[hn_ree_data2, hn_ree_var2, hn_ree_removed2, hn_ree_idx_removed2] = ...
    lessNaNs(hn_ree_data,hn_ree_var,0.89);

clear hn_data hn_var
Hn_ree_table = createNaNTable(hn_ree_data,hn_ree_var,hn_ree_obj);
clear hn_ree_data hn_ree_var

Hn_table1 = createNaNTable(hn_data1,hn_var1,hn_obj);
Hn_ree_table1 = createNaNTable(hn_ree_data1,hn_ree_var1,hn_ree_obj);

Hn_table2 = createNaNTable(hn_data2,hn_var2,hn_obj);
Hn_ree_table2 = createNaNTable(hn_ree_data2,hn_ree_var2,hn_ree_obj);


% writetable(Hn_table,'Hailniemi_NaN.xlsx','Sheet','GeoChem','WriteRowNames',true);
% writetable(Hn_ree_table,'Hailniemi_NaN.xlsx','Sheet','GeoChemREE','WriteRowNames',true);
% writetable(Hn_table1,'Hailniemi_NaN.xlsx','Sheet','GeoChem_step1','WriteRowNames',true);
% writetable(Hn_ree_table1,'Hailniemi_NaN.xlsx','Sheet','GeoChemREE_step1','WriteRowNames',true);
%% Hameenkyla
load('Hameenkyla2.mat')
[hk_data1, hk_var1, hk_removed1, hk_idx_removed1] = ...
    allNaNs(hk_data,hk_var);

[hk_ree_data1, hk_ree_var1, hk_ree_removed1, hk_ree_idx_removed1] = ...
    allNaNs(hk_ree_data,hk_ree_var);

Hk_table = createNaNTable(hk_data,hk_var,hk_obj);

[hk_data2, hk_var2, hk_removed2, hk_idx_removed2] = ...
    lessNaNs(hk_data,hk_var,0.89);

[hk_ree_data2, hk_ree_var2, hk_ree_removed2, hk_ree_idx_removed2] = ...
    lessNaNs(hk_ree_data,hk_ree_var,0.89);

clear hk_data hk_var
Hk_ree_table = createNaNTable(hk_ree_data,hk_ree_var,hk_ree_obj);
clear hk_ree_data hk_ree_var

Hk_table1 = createNaNTable(hk_data1,hk_var1,hk_obj);
Hk_ree_table1 = createNaNTable(hk_ree_data1,hk_ree_var1,hk_ree_obj);

Hk_table2 = createNaNTable(hk_data2,hk_var2,hk_obj);
Hk_ree_table2 = createNaNTable(hk_ree_data2,hk_ree_var2,hk_ree_obj);


% These don't work, NaN values missing and decimal control needs
% improvement
% writetable(Hk_table,'Hameenkyla_NaN.xlsx','Sheet','GeoChem','WriteRowNames',true);
% writetable(Hk_ree_table,'Hameenkyla_NaN.xlsx','Sheet','GeoChemREE','WriteRowNames',true);
% writetable(Hk_table1,'Hameenkyla_NaN.xlsx','Sheet','GeoChem_step1','WriteRowNames',true);
% writetable(Hk_ree_table1,'Hameenkyla_NaN.xlsx','Sheet','GeoChemREE_step1','WriteRowNames',true);

%% Heponiemi
load('HeponiemiREE2.mat')

[h_ree_data1, h_ree_var1, h_ree_removed1, h_ree_idx_removed1] = ...
    allNaNs(h_ree_data,h_ree_var);

H_ree_table = createNaNTable(h_ree_data,h_ree_var,h_ree_obj);

[h_ree_data2, h_ree_var2, h_ree_removed2, h_ree_idx_removed2] = ...
    lessNaNs(h_ree_data,h_ree_var,0.89);

clear h_ree_data h_ree_var

H_ree_table1 = createNaNTable(h_ree_data1,h_ree_var1,h_ree_obj);

H_ree_table2 = createNaNTable(h_ree_data2,h_ree_var2,h_ree_obj);
%% Heponiemi 2
load('Heponiemi2.mat')

[hepo_data1, hepo_var1, hepo_removed1, hepo_idx_removed1] = ...
    allNaNs(hepo_data,hepo_var);

[hepo_ree_data1, hepo_ree_var1, hepo_ree_removed1, hepo_ree_idx_removed1] = ...
    allNaNs(h_ree_data,h_ree_var);

hepo_table = createNaNTable(hepo_data,hepo_var,hepo_obj);

[hepo_data2, hepo_var2, hepo_removed2, hepo_idx_removed2] = ...
    lessNaNs(hepo_data,hepo_var,0.89);

[hepo_ree_data2, hepo_ree_var2, h_ree_removed2, hepo_ree_idx_removed2] = ...
    lessNaNs(h_ree_data,h_ree_var,0.89);

clear hepo_data hepo_var
hepo_ree_table = createNaNTable(h_ree_data,h_ree_var,hepo_ree_obj);
clear h_ree_data h_ree_var

hepo_table1 = createNaNTable(hepo_data1,hepo_var1,hepo_obj);
hepo_ree_table1 = createNaNTable(hepo_ree_data1,hepo_ree_var1,hepo_ree_obj);

hepo_table2 = createNaNTable(hepo_data2,hepo_var2,hepo_obj);
hepo_ree_table2 = createNaNTable(hepo_ree_data2,hepo_ree_var2,hepo_ree_obj);


%% Sample 1
load('Sample1.mat')
[s1_data1, s1_var1, s1_removed1, s1_idx_removed1] = ...
    allNaNs(s1_data,s1_var);

[s1_ree_data1, s1_ree_var1, s1_ree_removed1, s1_ree_idx_removed1] = ...
    allNaNs(s1_ree_data,s1_ree_var);

s1_table = createNaNTable(s1_data,s1_var,s1_obj);

[s1_data2, s1_var2, s1_removed2, s1_idx_removed2] = ...
    lessNaNs(s1_data,s1_var,0.89);

[s1_ree_data2, s1_ree_var2, s1_ree_removed2, s1_ree_idx_removed2] = ...
    lessNaNs(s1_ree_data,s1_ree_var,0.89);

clear s1_data s1_var
s1_ree_table = createNaNTable(s1_ree_data,s1_ree_var,s1_ree_obj);
clear s1_ree_data s1_ree_var

s1_table1 = createNaNTable(s1_data1,s1_var1,s1_obj);
s1_ree_table1 = createNaNTable(s1_ree_data1,s1_ree_var1,s1_ree_obj);

s1_table2 = createNaNTable(s1_data2,s1_var2,s1_obj);
s1_ree_table2 = createNaNTable(s1_ree_data2,s1_ree_var2,s1_ree_obj);


%% Sample 2
load('Sample2.mat')
[s2_data1, s2_var1, s2_removed1, s2_idx_removed1] = ...
    allNaNs(s2_data,s2_var);

[s2_ree_data1, s2_ree_var1, s2_ree_removed1, s2_ree_idx_removed1] = ...
    allNaNs(s2_ree_data,s2_ree_var);

s2_table = createNaNTable(s2_data,s2_var,s2_obj);

[s2_data2, s2_var2, s2_removed2, s2_idx_removed2] = ...
    lessNaNs(s2_data,s2_var,0.89);

[s2_ree_data2, s2_ree_var2, s2_ree_removed2, s2_ree_idx_removed2] = ...
    lessNaNs(s2_ree_data,s2_ree_var,0.89);

clear s2_data s2_var
s2_ree_table = createNaNTable(s2_ree_data,s2_ree_var,s2_ree_obj);
clear s2_ree_data s2_ree_var

s2_table1 = createNaNTable(s2_data1,s2_var1,s2_obj);
s2_ree_table1 = createNaNTable(s2_ree_data1,s2_ree_var1,s2_ree_obj);

s2_table2 = createNaNTable(s2_data2,s2_var2,s2_obj);
s2_ree_table2 = createNaNTable(s2_ree_data2,s2_ree_var2,s2_ree_obj);
%%
tot_ree = H_ree_table(end,:).Variables+...
    Hk_ree_table(end,:).Variables+...
    Hn_ree_table(end,:).Variables+...
    s1_ree_table(end,:).Variables+...
    s2_ree_table(end,:).Variables;
tot_ree(2,:) = tot_ree(1,:)/211;

tot = hepo_table(end,:).Variables+...
    Hk_table(end,:).Variables+...
    Hn_table(end,:).Variables+...
    s1_table(end,:).Variables+...
    s2_table(end,:).Variables;
tot(2,:) = tot(1,:)/210;
function [newData, newVars, vars_removed, idx_removed] = allNaNs(data,vars)
% Check if all observations in column have NaN
idx_removed = all(isnan(data),1);

newVars = vars(~idx_removed);
newData = data(:,~idx_removed);

vars_removed = vars(idx_removed);

end

function [newData, newVars, vars_removed, idx_removed] = lessNaNs(data,vars,threshold)

% Remove variables with only few samples.
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

function Table = createNaNTable(data,var,obj)
% Percent
% data(end+1,:) = sum(isnan(data))/size(data,1);
% Absolute value
data(end+1,:) = sum(isnan(data));
obj(end+1) = 'col_NaNs';

% data(:,end+1) = sum(isnan(data),2)/size(data,2);
data(:,end+1) = sum(isnan(data),2);
var(end+1) = cellstr('row_NaNs');
Table = array2table(data,'VariableNames',var,'RowNames',obj);

end