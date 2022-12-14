%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: C:\Users\Allu\Desktop\Töissä\Masters Thesis\Data\XRF\NaStA_PaavoHärmä_XRF_Mittaukset_Heponiemi_1_chemistry-ExportData-2020-07-03.xls
%    Worksheet: NaStA_PaavoHärmä_XRF_Mittaukset
%
% Auto-generated by MATLAB on 01-Mar-2022 13:25:11

%% Set up the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 91);

% Specify sheet and range
opts.Sheet = "NaStA_PaavoHärmä_XRF_Mittaukset";
opts.DataRange = "A2:CM82";

% Specify column names and types
opts.VariableNames = ["Var1", "Var2", "Var3", "Var4", "SampleID", "Var6", "Var7", "MgConcentration", "MgError1s", "AlConcentration", "AlError1s", "SiConcentration", "SiError1s", "PConcentration", "PError1s", "SConcentration", "SError1s", "KConcentration", "KError1s", "CaConcentration", "CaError1s", "TiConcentration", "TiError1s", "VConcentration", "VError1s", "CrConcentration", "CrError1s", "MnConcentration", "MnError1s", "FeConcentration", "FeError1s", "CoConcentration", "CoError1s", "NiConcentration", "NiError1s", "CuConcentration", "CuError1s", "ZnConcentration", "ZnError1s", "AsConcentration", "AsError1s", "SeConcentration", "SeError1s", "RbConcentration", "RbError1s", "SrConcentration", "SrError1s", "YConcentration", "YError1s", "ZrConcentration", "ZrError1s", "NbConcentration", "NbError1s", "MoConcentration", "MoError1s", "AgConcentration", "AgError1s", "CdConcentration", "CdError1s", "SnConcentration", "SnError1s", "SbConcentration", "SbError1s", "BaConcentration", "BaError1s", "LaConcentration", "LaError1s", "CeConcentration", "CeError1s", "PrConcentration", "PrError1s", "NdConcentration", "NdError1s", "TaConcentration", "TaError1s", "WConcentration", "WError1s", "AuConcentration", "AuError1s", "HgConcentration", "HgError1s", "PbConcentration", "PbError1s", "BiConcentration", "BiError1s", "ThConcentration", "ThError1s", "UConcentration", "UError1s", "LEConcentration", "LEError1s"];
opts.SelectedVariableNames = ["SampleID", "MgConcentration", "MgError1s", "AlConcentration", "AlError1s", "SiConcentration", "SiError1s", "PConcentration", "PError1s", "SConcentration", "SError1s", "KConcentration", "KError1s", "CaConcentration", "CaError1s", "TiConcentration", "TiError1s", "VConcentration", "VError1s", "CrConcentration", "CrError1s", "MnConcentration", "MnError1s", "FeConcentration", "FeError1s", "CoConcentration", "CoError1s", "NiConcentration", "NiError1s", "CuConcentration", "CuError1s", "ZnConcentration", "ZnError1s", "AsConcentration", "AsError1s", "SeConcentration", "SeError1s", "RbConcentration", "RbError1s", "SrConcentration", "SrError1s", "YConcentration", "YError1s", "ZrConcentration", "ZrError1s", "NbConcentration", "NbError1s", "MoConcentration", "MoError1s", "AgConcentration", "AgError1s", "CdConcentration", "CdError1s", "SnConcentration", "SnError1s", "SbConcentration", "SbError1s", "BaConcentration", "BaError1s", "LaConcentration", "LaError1s", "CeConcentration", "CeError1s", "PrConcentration", "PrError1s", "NdConcentration", "NdError1s", "TaConcentration", "TaError1s", "WConcentration", "WError1s", "AuConcentration", "AuError1s", "HgConcentration", "HgError1s", "PbConcentration", "PbError1s", "BiConcentration", "BiError1s", "ThConcentration", "ThError1s", "UConcentration", "UError1s", "LEConcentration", "LEError1s"];
opts.VariableTypes = ["char", "char", "char", "char", "string", "char", "char", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "categorical", "double", "categorical", "double", "double", "double", "double", "double", "categorical", "double", "categorical", "double", "categorical", "double", "double", "double", "categorical", "double", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "categorical", "double", "categorical", "double", "categorical", "double", "categorical", "double", "double", "double", "categorical", "double", "double", "double", "categorical", "double", "categorical", "double", "categorical", "double", "categorical", "double", "categorical", "double", "categorical", "double", "double", "double", "categorical", "double", "categorical", "double", "categorical", "double", "double", "double"];

% Replace cats with double
cat_idx = contains(opts.VariableTypes,'categorical');
opts.VariableTypes(cat_idx) = {'double'};

% Specify variable properties
opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var4", "SampleID", "Var6", "Var7"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var4", "SampleID", "Var6", "Var7", "MgConcentration", "VConcentration", "CrConcentration", "CoConcentration", "NiConcentration", "CuConcentration", "AsConcentration", "SeConcentration", "AgConcentration", "CdConcentration", "SnConcentration", "SbConcentration", "LaConcentration", "PrConcentration", "NdConcentration", "TaConcentration", "WConcentration", "AuConcentration", "HgConcentration", "BiConcentration", "ThConcentration", "UConcentration"], "EmptyFieldRule", "auto");

% Import the data
Hepo = readtable("C:\Users\Allu\Desktop\Töissä\Masters Thesis\Data\XRF\NaStA_PaavoHärmä_XRF_Mittaukset_Heponiemi_1_chemistry-ExportData-2020-07-03.xls", opts, "UseExcel", false);

% Error idx from variables
Errors_idx = contains(opts.SelectedVariableNames,'Error');
Errors_idx(1) = '';
% REE index
REE = contains(Hepo.SampleID,'R'); 
%%
Hepo_ree = Hepo(REE,:);
% Get object/sample names
Hepo_REE_obj = Hepo_ree.SampleID;
% Fetch only the numbers of the samples
[~,remain] = strtok(Hepo_REE_obj,{'1','2','3','4'});
[token1,remain1] = strtok(remain,'R');
[~,remain2] = strtok(remain,{'1','2','3','4','5','6','7','8','9'});
% Write objects as hn_ree_<square>_<number>
hepo_ree_obj = 'h_r_'+token1+'_'+remain2;

% Table to matrix
h_ree_data = Hepo_ree(:,2:end).Variables;
% Remove errors
h_ree_data(:,Errors_idx) = '';

h_ree_var = Hepo_ree.Properties.VariableNames(2:end);
% Remove errors
h_ree_var(Errors_idx) = '';
% Erase "Concentration" from variable names
h_ree_var = erase(h_ree_var,"Concentration");
%% Data without REE
% Data without REE
Hepo_ = Hepo(~REE,:);

% Get object/sample names
Hepo_obj = Hepo_.SampleID;
% Fetch only the numbers of the samples
[~,remain] = strtok(Hepo_obj,{'1','2','3','4'});
[token1,remain1] = strtok(remain,'G');
[~,remain2] = strtok(remain,{'1','2','3','4','5','6','7','8','9'});
% Write objects as h_<square>_<number>
hepo_obj = 'h_'+token1+'_'+remain2;

% Table to matrix
hepo_data = Hepo_(:,2:end).Variables;
% Remove errors
hepo_data(:,Errors_idx) = '';

hepo_var = Hepo_.Properties.VariableNames(2:end);
% Remove errors
hepo_var(Errors_idx) = '';
% Erase "Concentration" from variable names
hepo_var = erase(hepo_var,"Concentration");
%% Clear temporary variables
clear opts