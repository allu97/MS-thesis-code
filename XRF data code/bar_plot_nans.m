clearvars;clc
load('missingDataSummary.mat')
ree_summary = [H_ree_table{end,1:42};
    Hn_ree_table{end,1:42};Hk_ree_table{end,1:42};
    s1_ree_table{end,1:42};s2_ree_table{end,1:42}];
g_summary = [hepo_table{end,1:42};
    Hn_table{end,1:42};Hk_table{end,1:42};
    s1_table{end,1:42};s2_table{end,1:42}];
hk_var_bold = {'\textbf{Ba}', ...
    '\textbf{La}', ...
    '\textbf{Ce}', ...
    '\textbf{Pr}', ...
    '\textbf{Nd}', ...
    '\textbf{Ta}'};
hk_var_ree_bold = hk_var;
hk_var_ree_bold(29:34) = hk_var_bold;
% X = categorical(hk_var);
% X = reordercats(X,hk_var);
X = categorical(hk_var_ree_bold);
X = reordercats(X,hk_var_ree_bold);

% set(gca, 'TickLabelInterpreter', 'latex')
% set(gca,'fontweigth','bold')
%% REE data
f = figure;
t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile
b = barh(X,ree_summary,'stacked');
legend('Heponiemi','Hailniemi','Hämeenkylä',...
'Sample 1','Sample 2','location','east','FontSize',11,'FontName','times')
title('GEOChemREE missing value summary','FontName','Times')
xticks(0:25:225)
% align numbers to left (negative) or right (positive)
k = -8;
% limit for how many NaNs for sample is printed
lim = 9;

xtips1 = b(1).YEndPoints + k;
ytips1 = b(1).XEndPoints;
labels1 = string(b(1).YData);
labels1(str2double(labels1)<lim) = '';
text(xtips1,ytips1,labels1,'VerticalAlignment','middle','color','w')

xtips2 = b(2).YEndPoints + k;
ytips2 = b(2).XEndPoints;
labels2 = string(b(2).YData);
labels2(str2double(labels2)<lim) = '';
text(xtips2,ytips2,labels2,'VerticalAlignment','middle')

xtips3 = b(3).YEndPoints + k;
ytips3 = b(3).XEndPoints;
labels3 = string(b(3).YData);
labels3(str2double(labels3)<lim) = '';
text(xtips3,ytips3,labels3,'VerticalAlignment','middle')

xtips4 = b(4).YEndPoints + k;
ytips4 = b(4).XEndPoints;
labels4 = string(b(4).YData);
labels4(str2double(labels4)<lim) = '';
text(xtips4,ytips4,labels4,'VerticalAlignment','middle','color','w')

xtips5 = b(5).YEndPoints + k;
ytips5 = b(5).XEndPoints;
labels5 = string(b(5).YData);
labels5(str2double(labels5)<lim) = '';
text(xtips5,ytips5,labels5,'VerticalAlignment','middle')

% text(xtips5+15,ytips5,num2cell(tot_ree(1,1:42)),'VerticalAlignment','middle','FontWeight','Bold')
f.Position = [50 50 f.Position(3) f.Position(4)*2];
set(gca, 'TickLabelInterpreter', 'latex')
%%
f = figure;
t = tiledlayout(1,1,'TileSpacing','Compact','Padding','Compact');
nexttile
b = barh(X,g_summary,'stacked');
legend('Heponiemi','Hailniemi','Hämeenkylä',...
'Sample 1','Sample 2','location','east','FontSize',11,'FontName','times')
title('GEOChem missing value summary','FontName','Times')
xticks(0:25:225)
% align numbers to left (negative) or right (positive)
k = -8;
% limit for how many NaNs for sample is printed
lim = 9;

xtips1 = b(1).YEndPoints + k;
ytips1 = b(1).XEndPoints;
labels1 = string(b(1).YData);
labels1(str2double(labels1)<lim) = '';
text(xtips1,ytips1,labels1,'VerticalAlignment','middle','color','w')

xtips2 = b(2).YEndPoints + k;
ytips2 = b(2).XEndPoints;
labels2 = string(b(2).YData);
labels2(str2double(labels2)<lim) = '';
text(xtips2,ytips2,labels2,'VerticalAlignment','middle')

xtips3 = b(3).YEndPoints + k;
ytips3 = b(3).XEndPoints;
labels3 = string(b(3).YData);
labels3(str2double(labels3)<lim) = '';
text(xtips3,ytips3,labels3,'VerticalAlignment','middle')

xtips4 = b(4).YEndPoints + k;
ytips4 = b(4).XEndPoints;
labels4 = string(b(4).YData);
labels4(str2double(labels4)<lim) = '';
text(xtips4,ytips4,labels4,'VerticalAlignment','middle','color','w')

xtips5 = b(5).YEndPoints + k;
ytips5 = b(5).XEndPoints;
labels5 = string(b(5).YData);
labels5(str2double(labels5)<lim) = '';
text(xtips5,ytips5,labels5,'VerticalAlignment','middle')

% text(xtips5+15,ytips5,num2cell(tot_ree(1,1:42)),'VerticalAlignment','middle','FontWeight','Bold')
f.Position = [50 50 f.Position(3) f.Position(4)*2];
set(gca, 'TickLabelInterpreter', 'latex')