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
avg_A = mean(Reflectances(idx_A,1476:1664));
avg_B = mean(Reflectances(idx_B,1476:1664));
avg_BG = mean(Reflectances(idx_BG,1476:1664));
avg_C = mean(Reflectances(idx_C,1476:1664));
figure; hold on; grid on
plot(Wavelenghts(1476:1664),avg_A,'Color',[0.9290 0.6940 0.1250],'linewidth',2);
plot(Wavelenghts(1476:1664),avg_B,'k','linewidth',2);
plot(Wavelenghts(1476:1664),avg_BG,'Color',[0.5 0.5 0.5],'linewidth',2);
plot(Wavelenghts(1476:1664),avg_C,'r','linewidth',2);
h = title('Average of cut range');
h.Interpreter = "latex"; h.FontSize = 12;
h = xlabel('Wavenumber $cm^{-1}$');
h.Interpreter = "latex"; h.FontSize = 12;
h = ylabel('Reflectance \%');
h.Interpreter = "latex"; h.FontSize = 12;
h = legend('Wyborgite','Black granite','Baltic green','Pyterlite','Location','northwest');
h.Interpreter = "latex"; h.FontSize = 12;

% If the x axis values need to be reversed (descending order)
set(gca, 'xdir', 'reverse')