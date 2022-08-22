% c = cvpartition(Class,'Holdout',0.25);
% testData = Reflectances(c.test,:);
% testData = Ref_d1_g(c.test,:);
% testData = Ref_d2_g(c.test,:);
testData = x1(c.test,:);
testClass = Class(c.test);

% trainData = Reflectances(~c.test,:);
% trainData = Ref_d1_g(~c.test,:);
% trainData = Ref_d2_g(~c.test,:);
trainData = x1(~c.test,:);
trainClass = Class(~c.test);
%%
load('testresults.mat', 'classes')
s1_indx = testClass == 1;
s2_indx = testClass == 2;
s3_indx = testClass == 3;
s4_indx = testClass == 4;
for i = 1:4
    [acc(i), sens(i), spec(i)] = tp_tn_rate(testClass==i,classes(:,2)==i|classes(:,3)==i);
end
fprintf('\nCorrectly labeled in Class 1: %f',sum(classes(s1_indx,2:3) == 1,'all')/size(classes(s1_indx,2),1)*100)
fprintf('\nCorrectly labeled in Class 2: %f',sum(classes(s2_indx,2:3) == 2,'all')/size(classes(s2_indx,2),1)*100)
fprintf('\nCorrectly labeled in Class 3: %f',sum(classes(s3_indx,2:3) == 3,'all')/size(classes(s3_indx,2),1)*100)
fprintf('\nCorrectly labeled in Class 4: %f',sum(classes(s4_indx,2:3) == 4,'all')/size(classes(s4_indx,2),1)*100)
fprintf('\n')
% [acc, sens, spec] = tp_tn_rate(testClass,classes);
%% testi
Cc = classes(:,2);
Cc(Cc==0) = classes(classes(:,2)==0,3);
C = confusionmat(testClass,Cc);
confusionchart(C)
fprintf('\nCorrectly labeled in Class 1: %f',C(1,1)/size(classes(s1_indx,2),1)*100);
fprintf('\nCorrectly labeled in Class 2: %f',C(2,2)/size(classes(s2_indx,2),1)*100);
fprintf('\nCorrectly labeled in Class 3: %f',C(3,3)/size(classes(s3_indx,2),1)*100);
fprintf('\nCorrectly labeled in Class 4: %f',C(4,4)/size(classes(s4_indx,2),1)*100);
fprintf('\n')


function [accuracy, sensitivity, specificity] = tp_tn_rate(Class,jj)
C = confusionmat(Class,jj);

% TP = C(1,1);
% FN = C(1,2);
% FP = C(2,1);
% TN = C(2,2);
TN = C(1,1);
FP = C(1,2);
FN = C(2,1);
TP = C(2,2);

% true positive rate / recall
sensitivity = TP/(TP+FN); 
% true negative rate
specificity = TN/(TN+FP);
accuracy = (TP+TN)/(TP+FP+FN+TN);
% precision = TP/(TP+FP);
% F1 = (2*precision*sensitivity)/(precision+sensitivity)
end