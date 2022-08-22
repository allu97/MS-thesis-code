close all;clearvars;clc
% load('QuarriesREENaNReplaced2Step7.mat')
load('QuarriesREENaNReplaced3Step6.mat')
%% Divided data to train and test
c = cvpartition(Class_ree,'Holdout',0.25);

testData = XC_ree(c.test,:);
testClass = Class_ree(c.test);

trainData = XC_ree(~c.test,:);
trainClass = Class_ree(~c.test);

c1_data = trainData(trainClass == 1,:);
c2_data = trainData(trainClass == 2,:);
c3_data = trainData(trainClass == 3,:);

[Z1,m1,s1] = zscore(c1_data);
[Z2,m2,s2] = zscore(c2_data);
[Z3,m3,s3] = zscore(c3_data);
%% Determine parameters for aSIMCA
% How many PCs can be used
PC_v = 1:min(size(XC_ree))-1;

% Decision thresold for residual distance
% alpha_v = [.9 .925 .95 .975 .990 .995 .999 .9995 .9999];
% alpha_v = .95:.005:.995;
alpha_v = 1:0.1:5;


% [coeff,score,latent,tsquared,explained,mu] = pca(Z1,'NumComponents',PCs,'centered',false);
% % Note that even when you specify a reduced component space, 
% % pca computes the T-squared values in the full space, using all components.
% % The T-squared value in the reduced space corresponds to the
% % Mahalanobis distance in the reduced space.
% tsqreduced = mahal(score,score);
% % Calculate the T-squared values in the discarded space by taking 
% % the difference of the T-squared values in the full space and 
% % Mahalanobis distance in the reduced space.
% tsqdiscarded = tsquared - tsqreduced;

%% Fit train data (same threshold method as in the tool)
% [~,~,spec1] = aSimca_iter2(Z1,m1,s1,trainData,trainClass == 1,PC_v,alpha_v);
% [~,~,spec2] = aSimca_iter2(Z2,m2,s2,trainData,trainClass == 2,PC_v,alpha_v);
% [~,~,spec3] = aSimca_iter2(Z3,m3,s3,trainData,trainClass == 3,PC_v,alpha_v);

%% Fit test data with Q threshold
% [sens1,spec1] = aSimca_iter(Z1,m1,s1,testData,Class1,PC_v,alpha_v,1);
% [sens2,spec2] = aSimca_iter(Z2,m2,s2,testData,Class2,PC_v,alpha_v,2);
% [sens3,spec3] = aSimca_iter(Z3,m3,s3,testData,Class3,PC_v,alpha_v,3);
%% Fit test data (same threshold method as in the tool)
% [acc1,sens1,~] = aSimca_iter2(Z1,m1,s1,testData,Class1,PC_v,alpha_v);
% [acc2,sens2,~] = aSimca_iter2(Z2,m2,s2,testData,Class2,PC_v,alpha_v);
% [acc3,sens3,~] = aSimca_iter2(Z3,m3,s3,testData,Class3,PC_v,alpha_v);
%%
c = cvpartition(Class_ree,'KFold',10);
numFolds = c.NumTestSets;

acc1 = zeros(length(alpha_v),length(PC_v),numFolds);
acc2 = zeros(size(acc1));
acc3 = zeros(size(acc1));

sens1 = zeros(size(acc1));
sens2 = zeros(size(acc1));
sens3 = zeros(size(acc1));

spec1 = zeros(size(acc1));
spec2 = zeros(size(acc1));
spec3 = zeros(size(acc1));

[n,m] = size(XC_ree);
for i = 1:numFolds
    testData = XC_ree(c.test(i),:);
    testClass = Class_ree(c.test(i));

    trainData = XC_ree(~c.test(i),:);
    trainClass = Class_ree(~c.test(i));

    c1_data = trainData(trainClass == 1,:);
    c2_data = trainData(trainClass == 2,:);
    c3_data = trainData(trainClass == 3,:);

    % Autoscaling
    [Z1,m1,s1] = zscore(c1_data);
    [Z2,m2,s2] = zscore(c2_data);
    [Z3,m3,s3] = zscore(c3_data);
    [~,~,spec1(:,:,i)] = aSimca_iter2(Z1,m1,s1,trainData,trainClass == 1,PC_v,alpha_v);
    [~,~,spec2(:,:,i)] = aSimca_iter2(Z2,m2,s2,trainData,trainClass == 2,PC_v,alpha_v);
    [~,~,spec3(:,:,i)] = aSimca_iter2(Z3,m3,s3,trainData,trainClass == 3,PC_v,alpha_v);

    Class1 = testClass==1;
    Class2 = testClass==2;
    Class3 = testClass==3;

    [acc1(:,:,i),sens1(:,:,i),~] = aSimca_iter2(Z1,m1,s1,testData,Class1,PC_v,alpha_v);
    [acc2(:,:,i),sens2(:,:,i),~] = aSimca_iter2(Z2,m2,s2,testData,Class2,PC_v,alpha_v);
    [acc3(:,:,i),sens3(:,:,i),~] = aSimca_iter2(Z3,m3,s3,testData,Class3,PC_v,alpha_v);
    
    % Mean centering
    m1 = mean(c1_data);
    m2 = mean(c2_data);
    m3 = mean(c3_data);
    
    s1 = ones(1,m);
    
    [~,~,spec1(:,:,i)] = aSimca_iter2(c1_data,m1,s1,trainData,trainClass == 1,PC_v,alpha_v);
    [~,~,spec2(:,:,i)] = aSimca_iter2(c2_data,m2,s1,trainData,trainClass == 2,PC_v,alpha_v);
    [~,~,spec3(:,:,i)] = aSimca_iter2(c3_data,m3,s1,trainData,trainClass == 3,PC_v,alpha_v);

    Class1 = testClass==1;
    Class2 = testClass==2;
    Class3 = testClass==3;

    [acc1(:,:,i),sens1(:,:,i),~] = aSimca_iter2(c1_data,m1,s1,testData,Class1,PC_v,alpha_v);
    [acc2(:,:,i),sens2(:,:,i),~] = aSimca_iter2(c2_data,m2,s1,testData,Class2,PC_v,alpha_v);
    [acc3(:,:,i),sens3(:,:,i),~] = aSimca_iter2(c3_data,m3,s1,testData,Class3,PC_v,alpha_v);
    
    % No pre-treatment
    m1 = zeros(1,m);
    
    s1 = ones(1,m);
    
    [~,~,spec1(:,:,i)] = aSimca_iter2(c1_data,m1,s1,trainData,trainClass == 1,PC_v,alpha_v);
    [~,~,spec2(:,:,i)] = aSimca_iter2(c2_data,m1,s1,trainData,trainClass == 2,PC_v,alpha_v);
    [~,~,spec3(:,:,i)] = aSimca_iter2(c3_data,m1,s1,trainData,trainClass == 3,PC_v,alpha_v);

    Class1 = testClass==1;
    Class2 = testClass==2;
    Class3 = testClass==3;

    [acc1(:,:,i),sens1(:,:,i),~] = aSimca_iter2(c1_data,m1,s1,testData,Class1,PC_v,alpha_v);
    [acc2(:,:,i),sens2(:,:,i),~] = aSimca_iter2(c2_data,m1,s1,testData,Class2,PC_v,alpha_v);
    [acc3(:,:,i),sens3(:,:,i),~] = aSimca_iter2(c3_data,m1,s1,testData,Class3,PC_v,alpha_v);
end
acc1 = mean(acc1,3);
acc2 = mean(acc2,3);
acc3 = mean(acc3,3);

sens1 = mean(sens1,3);
sens2 = mean(sens2,3);
sens3 = mean(sens3,3);

spec1 = mean(spec1,3);
spec2 = mean(spec2,3);
spec3 = mean(spec3,3);
%% Return max accuracy and its index
[val1, idx] = max(acc1(:));
% [val1, idx] = max(sens1(:));
[row1, col1] = ind2sub(size(acc1),idx);
sprintf("Class 1: Acc: %0.5f, sens: %0.5f, spec: %0.5f, alpha: %0.5f, PCs: %d", ...
    acc1(row1,col1),sens1(row1,col1), spec1(row1,col1), alpha_v(row1), PC_v(col1))

[val2, idx] = max(acc2(:));
% [val2, idx] = max(sens2(:));
[row2, col2] = ind2sub(size(acc2),idx);
sprintf("Class 2: Acc: %0.5f, sens: %0.5f, spec: %0.5f, alpha: %0.5f, PCs: %d", ...
    acc2(row2,col2), sens2(row2,col2), spec2(row2,col2), alpha_v(row2), PC_v(col2))

[val3, idx] = max(acc3(:));
% [val3, idx] = max(sens3(:));
[row3, col3] = ind2sub(size(acc3),idx);
sprintf("Class 3: Acc: %0.5f, sens: %0.5f, spec: %0.5f, alpha: %0.5f, PCs: %d", ...
    acc3(row3,col3), sens3(row3,col3), spec3(row3,col3), alpha_v(row3), PC_v(col3))

%% Plotting sensitivity and specificity
f = figure; subplot(1,2,1); f.Position(3) = f.Position(3)*2;
imagesc(interp2(sens1),'XData',PC_v,'YData',alpha_v); 
c = colorbar;
c.Label.String = 'Sensitivity';
title('Sensitivity surface modeling C1');xlabel('Number of PCs'); ylabel('alpha'); 

subplot(1,2,2)
imagesc(interp2(spec1),'XData',PC_v,'YData',alpha_v); 
c = colorbar;
c.Label.String = 'Specificity';
title('Specificity surface modeling C1');xlabel('Number of PCs'); ylabel('alpha'); 

f = figure; subplot(1,2,1); f.Position(3) = f.Position(3)*2;
imagesc(interp2(sens2),'XData',PC_v,'YData',alpha_v); 
c = colorbar;
c.Label.String = 'Sensitivity';
title('Sensitivity surface modeling C2');xlabel('Number of PCs'); ylabel('alpha'); 

subplot(1,2,2)
imagesc(interp2(spec2),'XData',PC_v,'YData',alpha_v); 
c = colorbar;
c.Label.String = 'Specificity';
title('Specificity surface modeling C2');xlabel('Number of PCs'); ylabel('alpha'); 

f = figure; subplot(1,2,1); f.Position(3) = f.Position(3)*2;
imagesc(interp2(sens3),'XData',PC_v,'YData',alpha_v); 
c = colorbar;
c.Label.String = 'Sensitivity';
title('Sensitivity surface modeling C3');xlabel('Number of PCs'); ylabel('alpha'); 

subplot(1,2,2)
imagesc(interp2(spec3),'XData',PC_v,'YData',alpha_v); 
c = colorbar;
c.Label.String = 'Specificity';
title('Specificity surface modeling C3');xlabel('Number of PCs'); ylabel('alpha'); 
%%
load('SamplesREENaNsReplaced3Step6.mat')

c1_data = XC_ree(Class_ree == 1,:);
c2_data = XC_ree(Class_ree == 2,:);
c3_data = XC_ree(Class_ree == 3,:);

[Z1,m1,s1] = zscore(c1_data);
[Z2,m2,s2] = zscore(c2_data);
[Z3,m3,s3] = zscore(c3_data);

% Autoscaling
si2norm_c = zeros(size(X_ree,1),3);
[si2norm_c(:,1),s0lim(1)] = Simca_classify(Z1,m1,s1,X_ree,PC_v(col1),alpha_v(row1));
[si2norm_c(:,2),s0lim(2)] = Simca_classify(Z2,m2,s2,X_ree,PC_v(col2),alpha_v(row2));
[si2norm_c(:,3),s0lim(3)] = Simca_classify(Z3,m3,s3,X_ree,PC_v(col3),alpha_v(row3));

Class_idx = zeros(size(X_ree,1),3);
Class_idx(si2norm_c(:,1) < s0lim(1),1) = 1;
Class_idx(si2norm_c(:,2) < s0lim(2),2) = 2;
Class_idx(si2norm_c(:,3) < s0lim(3),3) = 3;

function [si2norm,s0lim] = Simca_classify(Z,mean,std,testData,PCs,alpha)
[coeff,score,~] = pca(Z,'NumComponents',PCs,'centered',false);
Ns = length(score);
Nv = length(coeff);
Nf = (Ns-PCs-1)*(Nv-PCs);

% Residual variances in classes, s0^2
s02 = sum(sum((Z-score*coeff').^2))/Nf;
% Decision threshold tinv(alpha,Nf)
s0lim = alpha*s02;

Ztest = (testData-mean)./std;
% Test fit
si2t=sum(((Ztest-(Ztest*coeff)*coeff').^2/(Nv-PCs)),2);
si2norm=si2t./s0lim;

end

function [acc, sens, spec] = aSimca_iter2(Z,mean,std,testData,testClass,PC_v,alpha_v)
% This uses normalized residual variance in class
acc = zeros(length(alpha_v),length(PC_v));
sens = zeros(length(alpha_v),length(PC_v));
spec = zeros(length(alpha_v),length(PC_v));

i = 0;
for alpha = alpha_v
    i = i+1;
    j = 0;
    for PCs = PC_v
        j = j+1;
        [coeff,score,~] = pca(Z,'NumComponents',PCs,'centered',false);
        Ns = length(score);
        Nv = length(coeff);
        Nf = (Ns-PCs-1)*(Nv-PCs);

        % Residual variances in classes, s0^2
        s02 = sum(sum((Z-score*coeff').^2))/Nf;
        % Decision threshold 
        s0lim = alpha*s02;

        Ztest = (testData-mean)./std;
        % Test fit
        si2t=sum(((Ztest-(Ztest*coeff)*coeff').^2/(Nv-PCs)),2);
        si2norm=si2t./s0lim;
        % Get class
        [acc(i,j), sens(i,j), spec(i,j)] = tp_tn_rate(testClass,si2norm < s0lim);
    end
end
end

function [sens, spec] = aSimca_iter(Z,mean,std,testData,testClass,PC_v,alpha_v,class)
acc = zeros(length(alpha_v),length(PC_v));
sens = zeros(length(alpha_v),length(PC_v));
spec = zeros(length(alpha_v),length(PC_v));
Q = zeros(length(alpha_v),length(PC_v));
i = 0;
for alpha = alpha_v
    i = i+1;
    j = 0;
    for PCs = PC_v
        j = j+1;
        Q(i,j) = qThreshold(Z,PCs,alpha);
%         Q_test = zeros(size(testData,1),1);
        Q_test = qFit(Z,mean,std,testData,PCs);
        Class_idx = zeros(size(Q_test));
        Class_idx(Q_test < Q(i,j)) = class;
        % Create Positive/Negative response for each class
        [acc(i,j), sens(i,j), spec(i,j)] = tp_tn_rate(testClass,Class_idx==class);
    end
end
end

function [accuracy, sensitivity, specificity] = tp_tn_rate(Class,jj)
C = confusionmat(Class,jj);

TP = C(1,1);
FN = C(1,2);
FP = C(2,1);
TN = C(2,2);

% true positive rate / recall
sensitivity = TP/(TP+FN); 
% true negative rate
specificity = TN/(TN+FP);
accuracy = (TP+TN)/(TP+FP+FN+TN);
% precision = TP/(TP+FP);
% F1 = (2*precision*sensitivity)/(precision+sensitivity)
end

function Q_test = qFit(Z,mean,std,testData,PCs)
% Contributions, Q-values
coeff = pca(Z,'NumComponents',PCs,'centered',false);
[m,~]=size(coeff);
I=eye(m);
Ztest = (testData-mean)./std;
qcontr = ((I - coeff*coeff')*Ztest').^2;
Q_test = sum(qcontr);
end

function Q_alpha = qThreshold(Z,PCs,alpha)
% Calculates the decision threshold Q_alpha as proposed by Jackson and
% Mudholkar using number of PCs and decision threshold alpha
Z_cov = cov(Z);
[~,S,~] = svd(Z_cov);
[m,n] = size(Z);

temp = diag(S);
if n < m
    emod = temp(PCs+1:n,:);
else
    emod = temp(PCs+1:m,:);
end
% Sums of eigenvalues raised to the kth power for the secondary components
% that are excluded from the PCA model.
th1 = sum(emod);
th2 = sum(emod.^2);
th3 = sum(emod.^3);
h0 = 1 - ((2*th1*th3)/(3*th2^2));

% Standard normal deviate corresponding to the upper percentile
z_alpha = norminv(alpha);
Q_alpha = th1*(((z_alpha*sqrt(2*th2*h0^2)/th1) + 1 + th2*h0*(h0-1)/th1^2)^(1/h0));
end