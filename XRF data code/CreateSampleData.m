load('Sample1.mat')

[n,m] = size(s1_ree_data);

class_s_ree = ones(n,1);

X_ree = s1_ree_data;

X_ree_obj = s1_ree_obj;

X_ree_var = s1_var;

X_var = s1_var;

% No REE
[n,m] = size(s1_data);

class_s = ones(n,1);

X = s1_data;

X_obj = s1_obj;

load('Sample2.mat')

[n,m] = size(s2_ree_data);

class_s_ree = [class_s_ree; ones(n,1)*2];

X_ree = [X_ree;s2_ree_data];

X_ree_obj = [X_ree_obj;s2_ree_obj];

% No REE
[n,m] = size(s2_data);

class_s = [class_s;ones(n,1)*2];

X = [X;s2_data];

X_obj = [X_obj;s2_obj];

%%
% REE data
load('DataREENaNReplaced2.mat')

[n,m] = size(s1);

class_s_ree = ones(n,1);

X_ree = s1;

X_ree_obj = s1_ree_obj;

X_ree_var = hn_var;

[n,m] = size(s2);

class_s_ree = [class_s_ree; ones(n,1)*2];

X_ree_obj = [X_ree_obj;s2_ree_obj];

X_ree = [X_ree;s2];

subclass_s_ree = [ones(10,1); ones(10,1)*2; ones(10,1)*3; ones(11,1)*4];
subclass_s_ree = [subclass_s_ree; ones(n,1)*5];


%% Step 1: First NaN removal
% Replace NaNs with zero
idx1 = isnan(X_ree);
X_ree(idx1) = 0;
%% Step 2: Remove variables with all NaN/0 values
% V, Se, Ag, Cd, Pr, W, Au, Hg
idx2 = [9 18 25 26 32 35 36 37];
X_ree(:,idx2) = '';
X_ree_var(idx2) = '';
%% Step 3: Remove the variables as described by Paavo and Satu-Pia
% Mg, P, S, Ti, Cr, Mn, Co, Ni, Cu, As, Sb, Nd, Bi, U
idx3 = [1 4 5 8 9 10 12 13 14 16 24 28 31 33];
X_ree(:,idx3) = '';
X_ree_var(idx3) = '';
%% Step 4: Remove variables that GeoREE does not favor
% 50 kV measurement does not favor Cr,Mn,Fe,Ni,Co,Cu,Zn
% Warning: Conflicts with step 3!
idx4 = strcmp('Cr',X_ree_var);
idx4 = idx4 | strcmp('Mn',X_ree_var);
idx4 = idx4 | strcmp('Fe',X_ree_var);
idx4 = idx4 | strcmp('Ni',X_ree_var);
idx4 = idx4 | strcmp('Co',X_ree_var);
idx4 = idx4 | strcmp('Cu',X_ree_var);
idx4 = idx4 | strcmp('Zn',X_ree_var);

X_ree(:,idx4) = '';
X_ree_var(idx4) = '';
%% Real Step 4: Remove outliers Batch 1
s1_idx = [12 17 28 40];
X_ree_obj(s1_idx)

s2_idx = [56 61 62];
X_ree_obj(s2_idx)

idx4 = [s1_idx s2_idx];
X_ree(idx4,:) = '';
class_s_ree(idx4) = '';
X_ree_obj(idx4) = [];
%% Step 5: Remove outliers Batch 2
s1_idx = [4 10 21 29];
X_ree_obj(s1_idx)

s2_idx = [59 69];
X_ree_obj(s2_idx)

idx5 = [s1_idx s2_idx];
X_ree(idx5,:) = '';
class_s_ree(idx5) = '';
X_ree_obj(idx5) = [];
%% Step 6: Remove outliers Batch 3
s1_idx = [12 19];
X_ree_obj(s1_idx)

s2_idx = [40 46 52 53];
X_ree_obj(s2_idx)

idx6 = [s1_idx s2_idx];
X_ree(idx6,:) = '';
class_s_ree(idx6) = '';
X_ree_obj(idx6) = [];
%% Step 7: Remove outliers Batch 4
s1_idx = [18];
X_ree_obj(s1_idx)

s2_idx = [59 62];
X_ree_obj(s2_idx)

idx7 = [s1_idx s2_idx];
X_ree(idx7,:) = '';
class_s_ree(idx7) = '';
X_ree_obj(idx7) = [];
%% Step 7.5 Remove outliers from Th
s1_idx = [3 16 18 19 20];
X_ree_obj(s1_idx)

X_ree(s1_idx,:) = '';
class_s_ree(s1_idx) = '';
X_ree_obj(s1_idx) = [];
%% REAL REAL Step 4: Remove Mo, Sn
idx4 = [12 13];
X_ree(:,idx4) = '';
X_ree_var(idx4) = '';
%% Attempt 3 Step 5: Remove outliers Batch 1
s1_idx = [12 28];
X_ree_obj(s1_idx)

s2_idx = 61;
X_ree_obj(s2_idx)

idx5 = [s1_idx s2_idx];

X_ree(idx5,:) = '';
class_s_ree(idx5) = '';
X_ree_obj(idx5) = [];
%% Attempt 3 Step 6: Remove outliers Batch 2
s1_idx = [22 23];
X_ree_obj(s1_idx)

s2_idx = 73;
X_ree_obj(s2_idx)

idx6 = [s1_idx s2_idx];

X_ree(idx6,:) = '';
class_s_ree(idx6) = '';
X_ree_obj(idx6) = [];
%% GeoChem data
load('DataGNaNReplaced.mat')

[n,m] = size(s1);

class_s = ones(n,1);

X = s1;

X_obj = s1_obj;

X_var = hn_var;

[n,m] = size(s2);

class_s = [class_s; ones(n,1)*2];

X_obj = [X_obj;s2_obj];

X = [X;s2];

subclass_s = [ones(10,1); ones(10,1)*2; ones(10,1)*3; ones(11,1)*4];
subclass_s = [subclass_s; ones(n,1)*5];

%% Step 1: First NaN removal
% Replace NaNs with zero
idx1 = isnan(X);
X(idx1) = 0;
%% Step 2: Remove REE
idx2 = [29:34];
X(:,idx2) = '';
X_var(idx2) = '';
%% Step 3: Remove variables with all NaN/0 values
% Mg, P, S, Co, Se, Mo, Ag, Cd, Sn, Sb, W, Au, Bi, Th, U
idx3 = [1 4 5 13 18 24 25 26 27 28 29 30 33 34 35];
X(:,idx3) = '';
X_var(idx3) = '';
%% Step 4: Remove the variables as described by Paavo and Satu-Pia
% Mg, Cr, Sb, Nd, Bi
idx4 = [9 29];
X(:,idx4) = '';
X_var(idx4) = '';
%% Real step 4: Remove outliers from data
idx4 = [8 9 31 42 54 85];
X_obj(idx4)

X(idx4,:) = '';
class_s(idx4) = '';
subclass_s(idx4) = '';
X_obj(idx4) = [];
%% For REE fingerprint 
% Step 1: load samples from SamplesREENaNsReplacedStep1
load('SamplesREENaNsReplacedStep1.mat')
%% Step 2: Remove samples that have no measurements or just one or few
% Mg, Ti, V, Cr, Mn, As, Se, Ag, Cd, Sb, Pr, Nd, W, Au, Hg, Bi
idx2 = [1 8 9 10 11 17 18 25 26 28 32 33 35 36 37 39];
X_ree(:,idx2) = '';
X_ree_var(idx2) = '';
