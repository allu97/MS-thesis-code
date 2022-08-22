load('FTIRDataSamples.mat')
idx_A = Class == 1;
idx_B = Class == 2;
idx_BG = Class == 3;
idx_C = Class == 4;
%% Fisher-Snedecor variable F
% Both intra and extra-spectral variance have been calulated using the
% Bessel's correction n-1
% Interspectral variance (InterV) of all spectra 
InterVarianz = sum((Reflectances-mean(Reflectances)).^2,2)./(size(Reflectances,2)-1);
% Intraspectral variance (IntraV) of all spectra within one class 
IntraVarianz_A = sum((Reflectances(:,idx_A)-mean(Reflectances(:,idx_A))).^2,2)./(sum(idx_A)-1);
IntraVarianz_B = sum((Reflectances(:,idx_B)-mean(Reflectances(:,idx_B))).^2,2)./(sum(idx_B)-1);
IntraVarianz_BG = sum((Reflectances(:,idx_BG)-mean(Reflectances(:,idx_BG))).^2,2)./(sum(idx_BG)-1);
IntraVarianz_C = sum((Reflectances(:,idx_C)-mean(Reflectances(:,idx_C))).^2,2)./(sum(idx_C)-1);

F_A = IntraVarianz_A./InterVarianz;
F_B = IntraVarianz_B./InterVarianz;
F_BG = IntraVarianz_BG./InterVarianz;
F_C = IntraVarianz_C./InterVarianz;

% Compute Fisher criterion Fc
% Fc is the critical value provided by the Fisher-Snedecor tables by
% freedom degree number of interSP and by freedom degree number of IntraSP.
% Fc depend on the percentage of risk, generally set to 5%.
% If F < Fc the null hypothesis H0 is not verified and the classification is
% possible. In this case, IntraSP is significantly lower than the InterSP.
%% Plot individual samples
% add space after samples: A1-1, B1, BG1, C1-1
sample_str = contains(SampleNames,'A1');
low_idx = find(sample_str==1, 1 );
max_idx = find(sample_str==1, 1, 'last' );

plotFTIRsample(Reflectances,Wavelenghts,SampleNames,low_idx,max_idx)

function plotFTIRsample(Reflectances,Wavelengths,SampleNames,lower_limit,upper_limit)
% lower_limit = 603; upper_limit = 622;
figure; hold on; grid on
plot(Wavelengths,mean(Reflectances(lower_limit:upper_limit,:)),'linewidth',2);
s = "mean";
legend(s);
idx = 1;
for i = lower_limit:upper_limit
    idx = idx + 1;
    plot(Wavelengths,Reflectances(i,:))
    s(idx) = SampleNames(i);
    legend(s)
end
title(SampleNames(i))
end