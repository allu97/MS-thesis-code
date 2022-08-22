load('HeponiemiREE.mat');
XC_ree = h_ree_data;
Class_ree = ones(51,1);
XC_ree_obj = h_ree_obj;

load('Hailniemi.mat')
% REE data
[n,m] = size(hn_ree_data);
XC_ree = [XC_ree;hn_ree_data];

Class_ree = [Class_ree;ones(n,1)*2];

XC_ree_obj = [XC_ree_obj;hn_ree_obj];

var_ree = hn_var;
% No REE
[n,m] = size(hn_data);
XC = hn_data;

Class = ones(n,1);

XC_obj = hn_obj;

load('Hameenkyla.mat')
% REE data
[n,m] = size(hk_ree_data);

XC_ree = [XC_ree; hk_ree_data];

Class_ree = [Class_ree; ones(n,1)*3];

XC_ree_obj = [XC_ree_obj; hk_ree_obj];

var = hn_var;
% No REE data
[n,m] = size(hk_data);

XC = [XC; hk_data];

Class = [Class; ones(n,1)*2];

XC_obj = [XC_obj; hk_obj];

%%
% REE data
load('DataREENaNReplaced2.mat')
XC_ree = hepo;
Class_ree = ones(41,1);
XC_ree_obj = h_ree_obj;

[n,m] = size(hn);
XC_ree = [XC_ree;hn];

Class_ree = [Class_ree;ones(n,1)*2];

XC_ree_obj = [XC_ree_obj;hn_ree_obj];

var_ree = hn_var;

[n,m] = size(hk);

XC_ree = [XC_ree; hk];

Class_ree = [Class_ree; ones(n,1)*3];

XC_ree_obj = [XC_ree_obj; hk_ree_obj];

%% Step 1: First NaN removal
% Replace NaNs with zero
idx1 = isnan(XC_ree);
XC_ree(idx1) = 0;
%% Step 2: Remove variables with all NaN/0 values
% V, Se, Ag, Cd, Pr, W, Au, Hg
idx2 = [9 18 25 26 32 35 36 37];
XC_ree(:,idx2) = '';
var_ree(idx2) = '';
%% Step 3: Remove the variables as described by Paavo and Satu-Pia
% Mg, P, S, Ti, Cr, Mn, Co, Ni, Cu, As, Sb, Nd, Bi, U
idx3 = [1 4 5 8 9 10 12 13 14 16 24 28 31 33];
XC_ree(:,idx3) = '';
var_ree(idx3) = '';
%% Step 4: Remove variables that GeoREE does not favor
% 50 kV measurement does not favor Cr,Mn,Fe,Ni,Co,Cu,Zn
% Warning: Conflicts with step 3!
idx4 = strcmp('Cr',var_ree);
idx4 = idx4 | strcmp('Mn',var_ree);
idx4 = idx4 | strcmp('Fe',var_ree);
idx4 = idx4 | strcmp('Ni',var_ree);
idx4 = idx4 | strcmp('Co',var_ree);
idx4 = idx4 | strcmp('Cu',var_ree);
idx4 = idx4 | strcmp('Zn',var_ree);

XC_ree(:,idx4) = '';
var_ree(idx4) = '';
%% Real Step4: Remove outliers Batch 1
hepo_idx = [15 22 24 34 37];
% for documenting print object names
XC_ree_obj(hepo_idx)

hn_idx = [49 58 63 64 69 70 78 84];
XC_ree_obj(hn_idx)

hk_idx = [93 98 109 114 117 120];
XC_ree_obj(hk_idx)

idx4 = [hepo_idx hn_idx hk_idx];
XC_ree(idx4,:) = '';
Class_ree(idx4) = '';
XC_ree_obj(idx4) = [];
%% Step 5: Remove outliers Batch 2
hepo_idx = [1 6 14 26 36];
XC_ree_obj(hepo_idx)

hn_idx = [54 63 70 74];
XC_ree_obj(hn_idx)

hk_idx = [78 82 84 99 100 102];
XC_ree_obj(hk_idx)

idx5 = [hepo_idx hn_idx hk_idx];
XC_ree(idx5,:) = '';
Class_ree(idx5) = '';
XC_ree_obj(idx5) = [];
%% Step 6: Remove outliers Batch 3
hepo_idx = [1 18 26];
XC_ree_obj(hepo_idx)

hn_idx = [33 39 54];
XC_ree_obj(hn_idx)

hk_idx = [66 77 90 92];
XC_ree_obj(hk_idx)

idx6 = [hepo_idx hn_idx hk_idx];
XC_ree(idx6,:) = '';
Class_ree(idx6) = '';
XC_ree_obj(idx6) = [];
%% Step 7: Remove outliers Batch 4
hepo_idx = [2 19 23 25];
XC_ree_obj(hepo_idx)

hn_idx = [38 53 58];
XC_ree_obj(hn_idx)

hk_idx = [63 67 69];
XC_ree_obj(hk_idx)

idx6 = [hepo_idx hn_idx hk_idx];
XC_ree(idx6,:) = '';
Class_ree(idx6) = '';
XC_ree_obj(idx6) = [];
%% REAL REAL Step 4: Remove sample 24 and Mo, Sn
idx4 = [12 13];
XC_ree(:,idx4) = '';
var_ree(idx4) = '';

XC_ree(24,:) = '';
XC_ree_obj(24) = [];
Class_ree(24) = '';
% Save as QuarriesREENaNReplaced3Step4.mat
%% Attempt 3 Step 5: Remove outliers Batch 1
hepo_idx = [15 22];
XC_ree_obj(hepo_idx)

hn_idx = [62 68 69];
XC_ree_obj(hn_idx)

hk_idx = [97 113 116];
XC_ree_obj(hk_idx)

idx5 = [hepo_idx hn_idx hk_idx];
XC_ree(idx5,:) = '';
Class_ree(idx5) = '';
XC_ree_obj(idx5) = [];
%% Attempt 3 Step 6: Remove outliers Batch 2

hn_idx = [46 58 78 81];
XC_ree_obj(hn_idx)

hk_idx = [90 102 108 111];
XC_ree_obj(hk_idx)

idx6 = [hn_idx hk_idx];
XC_ree(idx6,:) = '';
Class_ree(idx6) = '';
XC_ree_obj(idx6) = [];
%%
% GeoChem data
load('DataGNaNReplaced.mat')
XC = hepo;
[n,m] = size(hepo);
Class = ones(n,1);
XC_obj = hepo_obj;

[n,m] = size(hn);
XC = [XC;hn];

Class = [Class;ones(n,1)*2];

XC_obj = [XC_obj;hn_obj];

var = hn_var;

[n,m] = size(hk);

XC = [XC; hk];

Class = [Class; ones(n,1)*3];

XC_obj = [XC_obj; hk_obj];

%% Step 1: First NaN removal
% Replace NaNs with zero
idx1 = isnan(XC);
XC(idx1) = 0;
%% Step 2: Remove REE (Ba, La, Ce, Pr, Nd, Ta)
idx2 = [29:34];
XC(:,idx2) = '';
var(idx2) = '';
%% Step 3: Remove variables many NaN/0 values
% Mg, P, S, Co, Se, Mo, Ag, Cd, Sn, Sb, W, Au, Bi, Th, U

idx3 = [1 4 5 13 18 24 25 26 27 28 29 30 33 34 35];
XC(:,idx3) = '';
var(idx3) = '';
%% Step 4: Remove the variables as described by Paavo and Satu-Pia
% Mg, Cr, Sb, Nd, Bi
idx4 = [9 29];
XC(:,idx4) = '';
var(idx4) = '';
%% Real Step 4: Remove outliers from data
idx4 = [22 23 66 67 76 111 113 126];
XC_obj(idx4)

XC(idx4,:) = '';
Class(idx4) = '';
XC_obj(idx4) = [];
