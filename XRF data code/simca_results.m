% load('TestResults.mat')
s1_indx = class_s == 1;
s2_indx = class_s == 2;
fprintf('Percent of Sample 1 inside Heponiemi Class: %f',sum(classes(s1_indx,2) == 1)/size(classes(s1_indx,2),1)*100)
fprintf('\nPercent of Sample 1 inside Hailniemi Class: %f',sum(classes(s1_indx,2) == 2)/size(classes(s1_indx,2),1)*100)
fprintf('\nPercent of Sample 1 inside Hämeenkylä Class: %f',sum(classes(s1_indx,2) == 3)/size(classes(s1_indx,2),1)*100)
fprintf('\n')
fprintf('\nPercent of Sample 2 inside Heponiemi Class: %f',sum(classes(s2_indx,2) == 1)/size(classes(s2_indx,2),1)*100)
fprintf('\nPercent of Sample 2 inside Hailniemi Class: %f',sum(classes(s2_indx,2) == 2)/size(classes(s2_indx,2),1)*100)
fprintf('\nPercent of Sample 2 inside Hämeenkylä Class: %f',sum(classes(s2_indx,2) == 3)/size(classes(s2_indx,2),1)*100)
fprintf('\n')
fprintf('\nPercent of Sample 1 closest to Heponiemi Class: %f',sum(classes(s1_indx,2:3) == 1,'all')/size(classes(s1_indx,2),1)*100)
fprintf('\nPercent of Sample 1 closest to Hailniemi Class: %f',sum(classes(s1_indx,2:3) == 2,'all')/size(classes(s1_indx,2),1)*100)
fprintf('\nPercent of Sample 1 closest to Hämeenkylä Class: %f',sum(classes(s1_indx,2:3) == 3,'all')/size(classes(s1_indx,2),1)*100)
fprintf('\n')
fprintf('\nPercent of Sample 2 closest to Heponiemi Class: %f',sum(classes(s2_indx,2:3) == 1,'all')/size(classes(s2_indx,2),1)*100)
fprintf('\nPercent of Sample 2 closest to Hailniemi Class: %f',sum(classes(s2_indx,2:3) == 2,'all')/size(classes(s2_indx,2),1)*100)
fprintf('\nPercent of Sample 2 closest to Hämeenkylä Class: %f',sum(classes(s2_indx,2:3) == 3,'all')/size(classes(s2_indx,2),1)*100)
