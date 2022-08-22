close all;clearvars;clc
load('FTIRforSIMCA.mat')

Yc = zeros(length(Class_XC),4);
Yc(Class_XC==1,1) = ones(sum(Class_XC==1),1);
Yc(Class_XC==2,2) = ones(sum(Class_XC==2),1);
Yc(Class_XC==3,3) = ones(sum(Class_XC==3),1);
Yc(Class_XC==4,4) = ones(sum(Class_XC==4),1);

Yt = zeros(length(Class_X),4);
Yt(Class_X==1,1) = ones(sum(Class_X==1),1);
Yt(Class_X==2,2) = ones(sum(Class_X==2),1);
Yt(Class_X==3,3) = ones(sum(Class_X==3),1);
Yt(Class_X==4,4) = ones(sum(Class_X==4),1);
%% Preprocess the data
% Standardize calibration data and save mean and std
XC_mean = mean(XC);
XC_c = XC - XC_mean;
XC_std = std(XC);
XC_s = XC_c ./ XC_std;

% Same for class labels
YC_mean = mean(Yc);
YC_c = Yc - YC_mean;
YC_std = std(Yc);
YC_s = YC_c ./ YC_std;

% Standardize test set with calibration set parameters
Xt_c = X - XC_mean;
Xt_s = Xt_c ./ XC_std;

Yt_c = Yt - YC_mean;
Yt_s = Yt_c ./ YC_std;

%% Calibration
n = min(size(XC));
ncomp = n; % number of components used
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(XC_s, YC_s, ncomp, 'cv',5);

% Plot MSE
figure,hold on, plot(0:ncomp,MSE(2,:),'-bo'),plot(0:ncomp,MSE(1,:),'-ro');
xlabel('Number of components'),ylabel('Estimated Mean Squared Prediction Error')
legend({'PLS: MSE in Y','PLS: MSE in X'})
% Plot biplot
% figure,biplot(XL(:,1:2),XS(:,1:2))
% Plot variance captured by X and Y
figure,hold on
plot(1:ncomp,100*cumsum(PCTVAR(1,:))/sum(PCTVAR(1,:)),'r')
plot(1:ncomp,100*cumsum(PCTVAR(2,:))/sum(PCTVAR(2,:)),'k')
xlabel('Number of components'),ylabel('Explained Variance')
legend({'PLS: Explained Variance in X', 'PLS: Explained Variance in Y'})
% Variable importance
[varsPLS, idxPLS] = sort(abs(BETA),'descend');

%% Testing
ncomp = 90; % number of components used 
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(XC_s, YC_s, ncomp);

yfitPLS = [ones(size(Xt_s,1),1) Xt_s]*BETA;

for m = 1:length(yfitPLS(:,1))
    % Loop through every sample and take the max value which indicates the
    % possible class label
    a = yfitPLS(m,:);
    clear idx; clear maxV;
    [maxV, idx] = max(a);
    for j = 1:length(yfitPLS(m,:))
        if idx == j
            yfitPLS(m,j) = 1;
            % Set the class to 1 if it was the index of max value
        else
            yfitPLS(m,j) = 0;
            % Else the class label is zero
        end
    end
end
% Compute Rsquared
TSS = sum((Yt - mean(Yt)).^2);
RSS = sum((Yt - yfitPLS).^2);
rsquared = 1-RSS/TSS;
% Create confusion matrix plot
Yfit = zeros(size(yfitPLS,1),1);
for i = 1:size(yfitPLS,1)
    if yfitPLS(i,1) == 1
        Yfit(i,1) = 1;
    elseif yfitPLS(i,2) == 1
        Yfit(i,1) = 2;
    elseif yfitPLS(i,3) == 1
        Yfit(i,1) = 3;
    else
        Yfit(i,1) = 4;
    end
end

YT = zeros(size(Yt,1),1);
for i = 1:size(yfitPLS,1)
    if Yt(i,1) == 1
        YT(i,1) = 1;
    elseif Yt(i,2) == 1
        YT(i,1) = 2;
    elseif Yt(i,3) == 1
        YT(i,1) = 3;
    else
        YT(i,1) = 4;
    end
end

C = confusionmat(YT,Yfit);
% figure,confusionchart(C,'RowSummary','row-normalized','ColumnSummary','column-normalized');
figure,confusionchart(C,'RowSummary','total-normalized','ColumnSummary','total-normalized');

% A column-normalized column summary displays the number of correctly and 
% incorrectly classified observations for each predicted class as 
% percentages of the number of observations of the corresponding predicted
% class.

% A row-normalized row summary displays the number of correctly and 
% incorrectly classified observations for each true class as
% percentages of the number of observations of the corresponding true class.

% Normalized PLS weights
