% Script for cleaning the FTIR spectra data file
% This fetches only the A, B, BG and C samples
load('FTIRDataUpdated.mat')

% Remove from end of file unnecessary samples
Reflectances(1126:end,:) = '';
SampleNames(1126:end) = [];
% Remove 'Blokki'
Reflectances(824:870,:) = '';
SampleNames(824:870,:) = [];

% Rename SampleNames
[token,remain] = strtok(SampleNames,'csv');
% Remove date, fetch the number in parantheses
[~,remain2] = strtok(remain,'(');
% remove .csv
[token3,~] = strtok(remain2,'.');

SampleNames = token+token3;

% Create class labels
% 1 = A (Wyborgite, Baltic Brown)
% 2 = B (Black granite)
% 3 = BG (Baltic Green) 
% 4 = C (Pyterlite, Carmen Red)
Class = [ones(1,299) ones(1,303)*2 ones(1,221)*3 ones(1,255)*4]';

%% Create a plot of average spectras

idx_A = Class == 1;
idx_B = Class == 2;
idx_BG = Class == 3;
idx_C = Class == 4;

avg_A = mean(Reflectances(idx_A,:));
avg_B = mean(Reflectances(idx_B,:));
avg_BG = mean(Reflectances(idx_BG,:));
avg_C = mean(Reflectances(idx_C,:));
figure; hold on; grid on
plot(Wavelenghts,avg_A,'Color',[0.9290 0.6940 0.1250],'linewidth',2);
plot(Wavelenghts,avg_B,'k','linewidth',2);
plot(Wavelenghts,avg_BG,'Color',[0.5 0.5 0.5],'linewidth',2);
plot(Wavelenghts,avg_C,'r','linewidth',2);
h = title('Average of full range');
h.Interpreter = "latex"; h.FontSize = 12;
h = xlabel('Wavenumber $cm^{-1}$');
h.Interpreter = "latex"; h.FontSize = 12;
h = ylabel('Reflectance \%');
h.Interpreter = "latex"; h.FontSize = 12;
h = legend('Wyborgite','Black granite','Baltic green','Pyterlite','Location','northwest');
h.Interpreter = "latex"; h.FontSize = 12;

% If the x axis values need to be reversed (descending order)
set(gca, 'xdir', 'reverse')
%%
figure; hold on; grid on
plot(Wavelenghts(1449:end),avg_A(1449:end),'Color',[0.9290 0.6940 0.1250],'linewidth',2);
plot(Wavelenghts(1449:end),avg_B(1449:end),'k','linewidth',2);
plot(Wavelenghts(1449:end),avg_BG(1449:end),'Color',[0.5 0.5 0.5],'linewidth',2);
plot(Wavelenghts(1449:end),avg_C(1449:end),'r','linewidth',2);
h = title('Average of range 650-1300 $cm^{-1}$');
h.Interpreter = "latex"; h.FontSize = 12;
h = xlabel('Wavenumber $cm^{-1}$');
h.Interpreter = "latex"; h.FontSize = 12;
h = ylabel('Reflectance \%');
h.Interpreter = "latex"; h.FontSize = 12;
h = legend('Wyborgite','Black granite','Baltic green','Pyterlite','Location','northeast');
h.Interpreter = "latex"; h.FontSize = 12;

% If the x axis values need to be reversed (descending order)
set(gca, 'xdir', 'reverse')

%% Create data for Simca

Wavelenghts = Wavelenghts(1449:end);
Reflectances = Reflectances(:,1449:end);
Wavelenghts = round(Wavelenghts);
Wavelenghts_cell = num2cell(Wavelenghts);

test_indx = false(length(SampleNames),1);
% A index A-1, ..., A-5 and A-15
test_indx([1:20 121:224],:) = true;
% B index B-10, ..., B-15
test_indx(320:429,:) = true;
% BG index BG-10, ..., BG-12
test_indx(623:682,:) = true;
% C index C-1, ..., C-5 and C-17
test_indx([824:838 954:1028],:) = true;

X = Reflectances(test_indx,:);
XC = Reflectances(~test_indx,:);

Class_XC = Class(~test_indx);
Class_X = Class(test_indx);

obj_XC = SampleNames(~test_indx);
obj_X = SampleNames(test_indx);

%%
figure; hold on; grid on
% A samples
plot(Wavelenghts(1449:end),mean(Reflectances(1:20,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(21:40,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(41:60,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(61:80,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(81:100,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(101:120,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(121:144,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(145:164,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(165:184,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(185:204,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(205:224,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(225:244,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(245:264,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(265:279,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(280:299,1449:end)),'linewidth',2);

h = title('Average of range 650-1300 $cm^{-1}$');
h.Interpreter = "latex"; h.FontSize = 12;
h = xlabel('Wavenumber $cm^{-1}$');
h.Interpreter = "latex"; h.FontSize = 12;
h = ylabel('Reflectance \%');
h.Interpreter = "latex"; h.FontSize = 12;
h = legend('A1', 'A10', 'A11', 'A12', 'A13', 'A14', 'A15', 'A2', 'A3', 'A4', ...
    'A5', 'A6', 'A7', 'A8', 'A9','Location','NW');
h.Interpreter = "latex"; h.FontSize = 12;

set(gca, 'xdir', 'reverse')

figure; hold on; grid on
% B samples
plot(Wavelenghts(1449:end),mean(Reflectances(300:319,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(320:339,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(340:359,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(360:379,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(380:399,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(400:419,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(420:439,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(440:459,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(460:480,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(481:501,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(502:521,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(522:541,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(542:562,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(563:582,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(583:602,1449:end)),'linewidth',2);

h = title('Average of range 650-1300 $cm^{-1}$');
h.Interpreter = "latex"; h.FontSize = 12;
h = xlabel('Wavenumber $cm^{-1}$');
h.Interpreter = "latex"; h.FontSize = 12;
h = ylabel('Reflectance \%');
h.Interpreter = "latex"; h.FontSize = 12;
h = legend('B1', 'B10', 'B11', 'B12', 'B13', 'B14', 'B15', 'B2', 'B3', 'B4', ...
    'B5', 'B6', 'B7', 'B8', 'B9','Location','NW');
h.Interpreter = "latex"; h.FontSize = 12;

set(gca, 'xdir', 'reverse')

figure; hold on; grid on
% BG samples
plot(Wavelenghts(1449:end),mean(Reflectances(603:622,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(623:642,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(643:662,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(663:682,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(683:702,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(703:722,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(723:743,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(744:763,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(764:783,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(784:803,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(804:823,1449:end)),'linewidth',2);

h = title('Average of range 650-1300 $cm^{-1}$');
h.Interpreter = "latex"; h.FontSize = 12;
h = xlabel('Wavenumber $cm^{-1}$');
h.Interpreter = "latex"; h.FontSize = 12;
h = ylabel('Reflectance \%');
h.Interpreter = "latex"; h.FontSize = 12;
h = legend('BG1', 'BG10', 'BG11', 'BG12', 'BG2', 'BG4', ...
    'BG5', 'BG6', 'BG7', 'BG8', 'BG9','Location','NW');
h.Interpreter = "latex"; h.FontSize = 12;

set(gca, 'xdir', 'reverse')

figure; hold on; grid on
% C samples
plot(Wavelenghts(1449:end),mean(Reflectances(824:838,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(839:858,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(859:878,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(879:898,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(899:918,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(919:938,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(939:953,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(954:968,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(969:983,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(984:998,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(999:1013,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(999:1013,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(1014:1028,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(1029:1043,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(1044:1058,1449:end)),'linewidth',2);
plot(Wavelenghts(1449:end),mean(Reflectances(1059:end,1449:end)),'linewidth',2);

h = title('Average of range 650-1300 $cm^{-1}$');
h.Interpreter = "latex"; h.FontSize = 12;
h = xlabel('Wavenumber $cm^{-1}$');
h.Interpreter = "latex"; h.FontSize = 12;
h = ylabel('Reflectance \%');
h.Interpreter = "latex"; h.FontSize = 12;
h = legend('C1', 'C10', 'C11', 'C12','C13','C14','C16','C17', 'C2', 'C4', ...
    'C5', 'C6', 'C7','C9','Location','NW');
h.Interpreter = "latex"; h.FontSize = 12;

set(gca, 'xdir', 'reverse')

%% Plot individual samples
% add space after samples: A1-1, B1, BG1, C1-1
sample_str = contains(SampleNames,'A1-1 ');
low_idx = find(sample_str==1, 1 );
max_idx = find(sample_str==1, 1, 'last' );

plotFTIRsample(Reflectances,Wavelenghts,SampleNames,low_idx,max_idx)

%% First derivative
% A fist-order derivative is the rate of change of absorbance with respect
% to wavelength. It passes through zero at the same wavelength as max of
% the absorbance band.
[n,m] = size(Reflectances);
temp = [Reflectances zeros(n,1)];
Ref_d1 = diff(Reflectances,1,2);
Ref_d1_g = gradient(Reflectances) ./ gradient(Wavelenghts)'*-1;

figure;plot(Ref_d1')
figure;plot(Ref_d1_g')
%%
figure; hold on; grid on
plot(Wavelenghts,mean(Ref_d1(idx_A,:)),'Color',[0.9290 0.6940 0.1250],'linewidth',2)
plot(Wavelenghts,mean(Ref_d1(idx_B,:)),'k','linewidth',2)
plot(Wavelenghts,mean(Ref_d1(idx_BG,:)),'Color',[0.5 0.5 0.5],'linewidth',2)
plot(Wavelenghts,mean(Ref_d1(idx_C,:)),'r','linewidth',1)

h = title('Average of $1^{st}$ derivative');
h.Interpreter = "latex"; h.FontSize = 12;
h = xlabel('Wavenumber $cm^{-1}$');
h.Interpreter = "latex"; h.FontSize = 12;
h = ylabel('Diff');
h.Interpreter = "latex"; h.FontSize = 12;
h = legend('Wyborgite','Black granite','Baltic green','Pyterlite','Location','northwest');
h.Interpreter = "latex"; h.FontSize = 12;
set(gca, 'xdir', 'reverse')
%% Second derivative
% Again two ways to compute. diff will exclude two variables, gradient is
% more modern way and will include every variable.
Ref_d2 = diff(Reflectances,2,2);
Ref_d2_g = gradient(Ref_d1_g) ./ gradient(Wavelenghts)'*-1;
figure;plot(Ref_d2')
figure;plot(Ref_d2_g')
%% mRMR
[idx, scores] = fscmrmr(Reflectances, Class);

fig = figure;
bar(scores(idx))
xlabel('Predictor Rank')
ylabel('Predictor Importance Score')
title("Predictor Importace Ranking")
xticks(1:length(idx))
xticklabels(Wavelenghts_cell)
%% loads plots
figure
t = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
nexttile
bar(Wavelenghts,model4.cloads(:,1));
title("Principal component 1");
xlabel("Wavelength")
ylabel("Loads")
set(gca, 'xdir', 'reverse')
nexttile

bar(Wavelenghts,model4.cloads(:,2));
title("Principal component 2");
xlabel("Wavelength")
ylabel("Loads")
set(gca, 'xdir', 'reverse')

%% Cut off 650-900 wavelengths
% These wavelengths had much noise and could be seen on feature weights
% Use the original data with all samples
Wavelenghts(217:end) = '';
Reflectances(:,217:end) = '';
% save data as FTIRDataSmall.mat
function plotFTIRsample(Reflectances,Wavelengths,SampleNames,lower_limit,upper_limit)
% lower_limit = 603; upper_limit = 622;
figure; hold on; grid on
plot(Wavelengths(1449:end),mean(Reflectances(lower_limit:upper_limit,1449:end)),'linewidth',2);
s = "mean";
legend(s);
idx = 1;
for i = lower_limit:upper_limit
    idx = idx + 1;
    plot(Wavelengths(1449:end),Reflectances(i,1449:end))
    s(idx) = SampleNames(i);
    legend(s)
end
title(SampleNames(i))
end