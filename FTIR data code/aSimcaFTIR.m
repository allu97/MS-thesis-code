% close all;clearvars;clc

XC = Reflectances;
% XC = x1;
% XC = Ref_d1_g;
% XC = Ref_d2_g;
% XC = Reflectances_w;
Class_XC = Class;
%% Determine parameters for aSIMCA
% How many PCs can be used
PC_v = 5:min(size(XC))-200;
% Decision thresold for residual distance
% alpha_v = 1:0.1:5;
% alpha_v = [.905:.01:0.995];
alpha_v = [.705:.05:0.995];

%%
c = cvpartition(Class_XC,'KFold',10);
numFolds = c.NumTestSets;

acc1 = zeros(length(alpha_v),length(PC_v),numFolds);
acc2 = zeros(size(acc1));
acc3 = zeros(size(acc1));
acc4 = zeros(size(acc1));

sens1 = zeros(size(acc1));
sens2 = zeros(size(acc1));
sens3 = zeros(size(acc1));
sens4 = zeros(size(acc1));

spec1 = zeros(size(acc1));
spec2 = zeros(size(acc1));
spec3 = zeros(size(acc1));
spec4 = zeros(size(acc1));

[n,m] = size(XC);

for i = 1:numFolds
    testData = XC(c.test(i),:);
    testClass = Class_XC(c.test(i));

    trainData = XC(~c.test(i),:);
    trainClass = Class_XC(~c.test(i));

    c1_data = trainData(trainClass == 1,:);
    c2_data = trainData(trainClass == 2,:);
    c3_data = trainData(trainClass == 3,:);
    c4_data = trainData(trainClass == 4,:);

    % Autoscaling 
%     [Z1,m1,s1] = zscore(c1_data);
%     [Z2,m2,s2] = zscore(c2_data);
%     [Z3,m3,s3] = zscore(c3_data);
%     [Z4,m4,s4] = zscore(c4_data);
%     [~,~,spec1(:,:,i)] = aSimca_iter2(Z1,m1,s1,trainData,trainClass == 1,PC_v,alpha_v);
%     [~,~,spec2(:,:,i)] = aSimca_iter2(Z2,m2,s2,trainData,trainClass == 2,PC_v,alpha_v);
%     [~,~,spec3(:,:,i)] = aSimca_iter2(Z3,m3,s3,trainData,trainClass == 3,PC_v,alpha_v);
%     [~,~,spec4(:,:,i)] = aSimca_iter2(Z4,m4,s4,trainData,trainClass == 4,PC_v,alpha_v);
% 
%     Class1 = testClass==1;
%     Class2 = testClass==2;
%     Class3 = testClass==3;
%     Class4 = testClass==4;
% 
%     [acc1(:,:,i),sens1(:,:,i),~] = aSimca_iter2(Z1,m1,s1,testData,Class1,PC_v,alpha_v);
%     [acc2(:,:,i),sens2(:,:,i),~] = aSimca_iter2(Z2,m2,s2,testData,Class2,PC_v,alpha_v);
%     [acc3(:,:,i),sens3(:,:,i),~] = aSimca_iter2(Z3,m3,s3,testData,Class3,PC_v,alpha_v);
%     [acc4(:,:,i),sens4(:,:,i),~] = aSimca_iter2(Z4,m4,s4,testData,Class4,PC_v,alpha_v);

    % Mean centering 
    m1 = mean(c1_data);
    m2 = mean(c2_data);
    m3 = mean(c3_data);
    m4 = mean(c4_data);
    
    s1 = ones(1,m);
    
    [~,~,spec1(:,:,i)] = aSimca_iter2(c1_data,m1,s1,trainData,trainClass == 1,PC_v,alpha_v);
    [~,~,spec2(:,:,i)] = aSimca_iter2(c2_data,m2,s1,trainData,trainClass == 2,PC_v,alpha_v);
    [~,~,spec3(:,:,i)] = aSimca_iter2(c3_data,m3,s1,trainData,trainClass == 3,PC_v,alpha_v);
    [~,~,spec4(:,:,i)] = aSimca_iter2(c4_data,m4,s1,trainData,trainClass == 4,PC_v,alpha_v);

    Class1 = testClass==1;
    Class2 = testClass==2;
    Class3 = testClass==3;
    Class4 = testClass==4;

    [acc1(:,:,i),sens1(:,:,i),~] = aSimca_iter2(c1_data,m1,s1,testData,Class1,PC_v,alpha_v);
    [acc2(:,:,i),sens2(:,:,i),~] = aSimca_iter2(c2_data,m2,s1,testData,Class2,PC_v,alpha_v);
    [acc3(:,:,i),sens3(:,:,i),~] = aSimca_iter2(c3_data,m3,s1,testData,Class3,PC_v,alpha_v);
    [acc4(:,:,i),sens4(:,:,i),~] = aSimca_iter2(c4_data,m4,s1,testData,Class4,PC_v,alpha_v);
    
    % No pre-treatment 
%     m1 = zeros(1,m);
%     
%     s1 = ones(1,m);
%     
%     [~,~,spec1(:,:,i)] = aSimca_iter2(c1_data,m1,s1,trainData,trainClass == 1,PC_v,alpha_v);
%     [~,~,spec2(:,:,i)] = aSimca_iter2(c2_data,m1,s1,trainData,trainClass == 2,PC_v,alpha_v);
%     [~,~,spec3(:,:,i)] = aSimca_iter2(c3_data,m1,s1,trainData,trainClass == 3,PC_v,alpha_v);
%     [~,~,spec4(:,:,i)] = aSimca_iter2(c4_data,m1,s1,trainData,trainClass == 4,PC_v,alpha_v);
% 
%     Class1 = testClass==1;
%     Class2 = testClass==2;
%     Class3 = testClass==3;
%     Class4 = testClass==4;
% 
%     [acc1(:,:,i),sens1(:,:,i),~] = aSimca_iter2(c1_data,m1,s1,testData,Class1,PC_v,alpha_v);
%     [acc2(:,:,i),sens2(:,:,i),~] = aSimca_iter2(c2_data,m1,s1,testData,Class2,PC_v,alpha_v);
%     [acc3(:,:,i),sens3(:,:,i),~] = aSimca_iter2(c3_data,m1,s1,testData,Class3,PC_v,alpha_v);
%     [acc4(:,:,i),sens4(:,:,i),~] = aSimca_iter2(c4_data,m1,s1,testData,Class4,PC_v,alpha_v);
end
acc1 = mean(acc1,3);
acc2 = mean(acc2,3);
acc3 = mean(acc3,3);
acc4 = mean(acc4,3);

sens1 = mean(sens1,3);
sens2 = mean(sens2,3);
sens3 = mean(sens3,3);
sens4 = mean(sens4,3);

spec1 = mean(spec1,3);
spec2 = mean(spec2,3);
spec3 = mean(spec3,3);
spec4 = mean(spec4,3);
%% Return max accuracy and its index
[val1, idx] = max(acc1(:));
% [val1, idx] = max(sens1(:));
% [val1, idx] = max(spec1(:));
[row1, col1] = ind2sub(size(acc1),idx);
% [row1, col1] = ind2sub(size(spec1),idx);
sprintf("Class 1: Acc: %0.5f, sens: %0.5f, spec: %0.5f, alpha: %0.5f, PCs: %d", ...
    acc1(row1,col1),sens1(row1,col1), spec1(row1,col1), alpha_v(row1), PC_v(col1))

[val2, idx] = max(acc2(:));
% [val2, idx] = max(sens2(:));
% [val2, idx] = max(spec2(:));
[row2, col2] = ind2sub(size(acc2),idx);
% [row2, col2] = ind2sub(size(spec2),idx);
sprintf("Class 2: Acc: %0.5f, sens: %0.5f, spec: %0.5f, alpha: %0.5f, PCs: %d", ...
    acc2(row2,col2), sens2(row2,col2), spec2(row2,col2), alpha_v(row2), PC_v(col2))

[val3, idx] = max(acc3(:));
% [val3, idx] = max(sens3(:));
% [val3, idx] = max(spec3(:));
[row3, col3] = ind2sub(size(acc3),idx);
% [row3, col3] = ind2sub(size(spec3),idx);
sprintf("Class 3: Acc: %0.5f, sens: %0.5f, spec: %0.5f, alpha: %0.5f, PCs: %d", ...
    acc3(row3,col3), sens3(row3,col3), spec3(row3,col3), alpha_v(row3), PC_v(col3))

[val4, idx] = max(acc4(:));
% [val4, idx] = max(sens4(:));
% [val4, idx] = max(spec4(:));
[row4, col4] = ind2sub(size(acc4),idx);
% [row4, col4] = ind2sub(size(spec4),idx);
sprintf("Class 4: Acc: %0.5f, sens: %0.5f, spec: %0.5f, alpha: %0.5f, PCs: %d", ...
    acc4(row4,col4), sens4(row4,col4), spec4(row4,col4), alpha_v(row4), PC_v(col4))

%% Plotting sensitivity and specificity
% f = figure; subplot(1,3,1); f.Position(3) = f.Position(3)*2;
% imagesc(interp2(sens1),'XData',PC_v,'YData',alpha_v); 
% hold on; plot(col1,alpha_v(row1),'r+');
% c = colorbar;
% c.Label.String = 'Sensitivity';
% title('Sensitivity surface modeling C1');xlabel('Number of PCs'); ylabel('alpha'); 
% 
% subplot(1,3,2)
% imagesc(interp2(spec1),'XData',PC_v,'YData',alpha_v);
% hold on; plot(col1,alpha_v(row1),'r+');
% c = colorbar;
% c.Label.String = 'Specificity';
% title('Specificity surface modeling C1');xlabel('Number of PCs'); ylabel('alpha'); 
% 
% subplot(1,3,3)
% imagesc(interp2(acc1),'XData',PC_v,'YData',alpha_v); 
% hold on; plot(col1,alpha_v(row1),'r+');
% c = colorbar;
% c.Label.String = 'Accuracy';
% title('Accuracy surface modeling C1');xlabel('Number of PCs'); ylabel('alpha'); 
% 
% 
% f = figure; subplot(1,3,1); f.Position(3) = f.Position(3)*2;
% imagesc(interp2(sens2),'XData',PC_v,'YData',alpha_v); 
% hold on; plot(col2,alpha_v(row2),'r+');
% c = colorbar;
% c.Label.String = 'Sensitivity';
% title('Sensitivity surface modeling C2');xlabel('Number of PCs'); ylabel('alpha'); 
% 
% subplot(1,3,2)
% imagesc(interp2(spec2),'XData',PC_v,'YData',alpha_v); 
% hold on; plot(col2,alpha_v(row2),'r+');
% c = colorbar;
% c.Label.String = 'Specificity';
% title('Specificity surface modeling C2');xlabel('Number of PCs'); ylabel('alpha'); 
% 
% subplot(1,3,3)
% imagesc(interp2(acc2),'XData',PC_v,'YData',alpha_v); 
% hold on; plot(col2,alpha_v(row2),'r+');
% c = colorbar;
% c.Label.String = 'Accuracy';
% title('Accuracy surface modeling C2');xlabel('Number of PCs'); ylabel('alpha'); 
% 
% f = figure; subplot(1,3,1); f.Position(3) = f.Position(3)*2;
% imagesc(interp2(sens3),'XData',PC_v,'YData',alpha_v); 
% hold on; plot(col3,alpha_v(row3),'r+');
% c = colorbar;
% c.Label.String = 'Sensitivity';
% title('Sensitivity surface modeling C3');xlabel('Number of PCs'); ylabel('alpha'); 
% 
% subplot(1,3,2)
% imagesc(interp2(spec3),'XData',PC_v,'YData',alpha_v); 
% hold on; plot(col3,alpha_v(row3),'r+');
% c = colorbar;
% c.Label.String = 'Specificity';
% title('Specificity surface modeling C3');xlabel('Number of PCs'); ylabel('alpha'); 
% 
% subplot(1,3,3)
% imagesc(interp2(acc3),'XData',PC_v,'YData',alpha_v); 
% hold on; plot(col3,alpha_v(row3),'r+');
% c = colorbar;
% c.Label.String = 'Accuracy';
% title('Accuracy surface modeling C3');xlabel('Number of PCs'); ylabel('alpha'); 
% 
% f = figure; subplot(1,3,1); f.Position(3) = f.Position(3)*2;
% imagesc(interp2(sens4),'XData',PC_v,'YData',alpha_v); 
% hold on; plot(col4,alpha_v(row4),'r+');
% c = colorbar;
% c.Label.String = 'Sensitivity';
% title('Sensitivity surface modeling C4');xlabel('Number of PCs'); ylabel('alpha'); 
% 
% subplot(1,3,2)
% imagesc(interp2(spec4),'XData',PC_v,'YData',alpha_v);
% hold on; plot(col4,alpha_v(row4),'r+');
% c = colorbar;
% c.Label.String = 'Specificity';
% title('Specificity surface modeling C4');xlabel('Number of PCs'); ylabel('alpha'); 
% 
% subplot(1,3,3)
% imagesc(interp2(acc4),'XData',PC_v,'YData',alpha_v); 
% hold on; plot(col4,alpha_v(row4),'r+');
% c = colorbar;
% c.Label.String = 'Accuracy';
% title('Accuracy surface modeling C4');xlabel('Number of PCs'); ylabel('alpha'); 
%% better version

% t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
% nexttile
% imagesc(interp2(sens1),'XData',PC_v,'YData',alpha_v); 
% hold on; plot(col1,alpha_v(row1),'r+');
% c = colorbar;
% title('Sensitivity surface modeling C1');xlabel('Number of PCs'); ylabel('alpha'); 
% 
% nexttile
% imagesc(interp2(spec1),'XData',PC_v,'YData',alpha_v);
% hold on; plot(col1,alpha_v(row1),'r+');
% c = colorbar;
% title('Specificity surface modeling C1');xlabel('Number of PCs'); ylabel('alpha'); 
% 
% nexttile
% imagesc(interp2(acc1),'XData',PC_v,'YData',alpha_v); 
% hold on; plot(col1,alpha_v(row1),'r+');
% c = colorbar;
% title('Accuracy surface modeling C1');xlabel('Number of PCs'); ylabel('alpha'); 

%%
c = cvpartition(Class,'Holdout',0.20);
trainData = XC(~c.test,:);
testData = XC(c.test,:);
trainClass = Class(~c.test);
testClass = Class(c.test);

c1_data = trainData(trainClass == 1,:);
c2_data = trainData(trainClass == 2,:);
c3_data = trainData(trainClass == 3,:);
c4_data = trainData(trainClass == 4,:);

si2norm_c = zeros(size(testData,1),4);

% Autoscaling
% [Z1,m1,s1] = zscore(c1_data);
% [Z2,m2,s2] = zscore(c2_data);
% [Z3,m3,s3] = zscore(c3_data);
% [Z4,m4,s4] = zscore(c4_data);
% 
% [si2norm_c(:,1),s0lim(1)] = Simca_classify(Z1,m1,s1,testData,PC_v(col1),alpha_v(row1));
% [si2norm_c(:,2),s0lim(2)] = Simca_classify(Z2,m2,s2,testData,PC_v(col2),alpha_v(row2));
% [si2norm_c(:,3),s0lim(3)] = Simca_classify(Z3,m3,s3,testData,PC_v(col3),alpha_v(row3));
% [si2norm_c(:,4),s0lim(4)] = Simca_classify(Z4,m4,s4,testData,PC_v(col4),alpha_v(row4));
% 
% [si2norm_c(:,1),s0lim(1)] = Simca_classify(Z1,m1,s1,testData,5,0.96);
% [si2norm_c(:,2),s0lim(2)] = Simca_classify(Z2,m2,s2,testData,3,0.97);
% [si2norm_c(:,3),s0lim(3)] = Simca_classify(Z3,m3,s3,testData,5,0.95);
% [si2norm_c(:,4),s0lim(4)] = Simca_classify(Z4,m4,s4,testData,5,.95);

% Class_idx = zeros(size(X,1),4);
% Class_idx(si2norm_c(:,1) <= 1,1) = 1;
% Class_idx(si2norm_c(:,2) <= 1,2) = 2;
% Class_idx(si2norm_c(:,3) <= 1,3) = 3;
% Class_idx(si2norm_c(:,4) <= 1,4) = 4;

% Mean centering
m1 = mean(c1_data);
m2 = mean(c2_data);
m3 = mean(c3_data);
m4 = mean(c4_data);

s1 = ones(1,m);

[si2norm_c(:,1),s0lim(1)] = Simca_classify(c1_data,m1,s1,testData,PC_v(col1),alpha_v(row1));
[si2norm_c(:,2),s0lim(2)] = Simca_classify(c2_data,m2,s1,testData,PC_v(col2),alpha_v(row2));
[si2norm_c(:,3),s0lim(3)] = Simca_classify(c3_data,m3,s1,testData,PC_v(col3),alpha_v(row3));
[si2norm_c(:,4),s0lim(4)] = Simca_classify(c4_data,m4,s1,testData,PC_v(col4),alpha_v(row4));

% No pre-treatment
% m1 = zeros(1,m);
% 
% s1 = ones(1,m);
% 
% [si2norm_c(:,1),s0lim(1)] = Simca_classify(c1_data,m1,s1,testData,PC_v(col1),alpha_v(row1));
% [si2norm_c(:,2),s0lim(2)] = Simca_classify(c2_data,m1,s1,testData,PC_v(col2),alpha_v(row2));
% [si2norm_c(:,3),s0lim(3)] = Simca_classify(c3_data,m1,s1,testData,PC_v(col3),alpha_v(row3));
% [si2norm_c(:,4),s0lim(4)] = Simca_classify(c4_data,m1,s1,testData,PC_v(col4),alpha_v(row4));

[~,Class_idx] = min(si2norm_c,[],2);

figure; hold on
for i = 1:4
    classes = num2str(i);
    plot(si2norm_c(:,i)','DisplayName',classes);
    legend
end
xlabel('Object number')
ylabel('Normalized distance')
set(gca,'yscale','log')
% hline(1);hold off

acc_test = sum(testClass == Class_idx)/length(testClass)

function [si2norm,s0lim] = Simca_classify(Z,mean,std,testData,PCs,alpha)
[coeff,score,~] = pca(Z,'NumComponents',PCs,'centered',false);
Ns = length(score);
Nv = length(coeff);
Nf = (Ns-PCs-1)*(Nv-PCs);

% Residual variances in classes, s0^2
s02 = sum(sum((Z-score*coeff').^2))/Nf;
% Decision threshold tinv(alpha,Nf)
s0lim = tinv(alpha,Nf)*s02;

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
        s0lim = tinv(alpha,Nf)*s02;

        Ztest = (testData-mean)./std;
        % Test fit
        si2t=sum(((Ztest-(Ztest*coeff)*coeff').^2/(Nv-PCs)),2);
        si2norm=si2t./s0lim;
%         [~,Class_idx] = min(si2norm,[],2);
        % Get class
%         [acc(i,j), sens(i,j), spec(i,j)] = tp_tn_rate(testClass,si2norm < s0lim);
        [acc(i,j), sens(i,j), spec(i,j)] = tp_tn_rate(testClass,si2norm < 1);
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