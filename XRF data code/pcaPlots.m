
%% Standardize
Z = zscore(X); % Standardized data
%% PCA matlab
[coeff,score,latent,tsquared,explained,mu] = pca(Z,'NumComponents',5);

figure,plot(1:length(tsquared),tsquared),xlabel('sample'),ylabel('T2 square scores')

figure,stairs(cumsum(explained))
%% Biplot
class1 = find(class==1);
class2 = find(class==2);
class3 = find(class==3);
c1 = 1;
c2 = 2;
% Biplot with scores and 
figure
biplot(coeff(:,[c1 c2]),'Scores',score(class1,[c1 c2]),'VarLabels',var)
hold on
biplot(coeff(:,[c1 c2]),'Scores',score(class2,[c1 c2]),'Color','k')
biplot(coeff(:,[c1 c2]),'Scores',score(class3,[c1 c2]),'Color','b')

%% outlier detection
threshold = 25;
tOver = find(tsquared > threshold);
figure,biplot(coeff(:,[c1 c2]),'scores',score(tOver,[c1 c2]));
