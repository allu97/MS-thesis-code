load('FTIRDataSamples.mat')
Reflectances = Reflectances(:,1449:end);

%% A samples
A_idx = [4, 28, 30, 53, 54, 62, 63, 65, 75, 79, 104, 107, 114,...
    117, 172, 173, 179, 192, 268, 271, 276, 277, 293, 297];
SampleNames(A_idx)'
%%
B_idx = [312, 320, 340, 341, 343, 349, 351, 370, 371, 372, 375,...
    390, 391, 392, 398, 399, 410, 416, 421, 423, 435, 466, 471,...
    507, 512, 537, 548, 556, 560, 571, 579, 581, 582, 598];
SampleNames(B_idx)'
%%
BG_idx = [606, 603, 621, 626, 634, 651, 653, 654, 659, 673, 674,...
    678, 684, 693, 694, 696, 702, 724, 728, 730, 740, 759, 760,...
    766, 772, 776, 783, 814, 820];
SampleNames(BG_idx)'
%%
C_idx = [838, 846, 848, 877, 881, 885, 906, 908, 909, 913, 914,...
    916, 917, 922, 950, 951, 954, 955, 967, 1000, 1001, 1006,...
    1023, 1035, 1064];
SampleNames(C_idx)'

idx1 = [A_idx B_idx BG_idx C_idx];
Reflectances(idx1,:) = '';
Class(idx1) = '';
SampleNames(idx1) = [];
%%
A_idx = [29, 33, 38, 42, 45, 50, 51, 52, 67, 95, 110, 137, 156,...
    165, 206, 227, 243, 245, 248, 249, 250, 251, 252, 258, 259];
SampleNames(A_idx)'
%%
B_idx = [295, 303, 326, 334, 339, 340, 369, 377, 379, 382, 383,...
    387, 392, 402, 411, 415, 419, 424, 448, 455, 465, 469, 470,...
    473, 507, 512, 533, 538];
SampleNames(B_idx)'
%%
BG_idx = [553, 564, 569, 571, 588, 592, 593, 613, 617, 623, 624,...
    625, 626, 627, 648, 653, 654, 679, 690, 714, 731, 732, 733];
SampleNames(BG_idx)'
%%
C_idx = [740, 757, 782, 783, 786, 802, 813, 815, 823, 824, 825,...
    826, 846, 847, 849, 862, 872, 874, 893, 897, 898, 917, 918,...
    924, 962, 963];
SampleNames(C_idx)'

idx2 = [A_idx B_idx BG_idx C_idx];
Reflectances(idx2,:) = '';
Class(idx2) = '';
SampleNames(idx2) = [];
%%
idx3 = [658, 701, 720 829, 853];
SampleNames(idx3)'
Reflectances(idx3,:) = '';
Class(idx3) = '';
SampleNames(idx3) = [];
%% For 1st derivative (mean centered)
% From log unconditional probability
idx1 = [62 63 192 654 688 913 944];
SampleNames(idx1)'
Reflectances(idx1,:) = '';
Ref_d1(idx1,:) = '';
Class(idx1) = '';
SampleNames(idx1) = [];
%% Step 2 FIXED for 1st derivative
% From log up
idx2 = [655 979 1027];
SampleNames(idx2)'
Reflectances(idx2,:) = '';
Ref_d1(idx2,:) = '';

Class(idx2) = '';
SampleNames(idx2) = [];
%% Use only Log UP for original reflectances
% Auto and mean center from log up (15 k-folds)
idx1 = [4 53 62 63 107 114 192 245 297 341 621 626 ...
    634 653 654 660 731 759 760 876 877 894 908 909 913];
SampleNames(idx1)'
Reflectances(idx1,:) = '';
Ref_d1(idx1,:) = '';
Class(idx1) = '';
SampleNames(idx1) = [];
%% Step 2 of using only Log UP for reflectance
idx2 = [27 166 167 173 686 714 741 819 850 873 942];
SampleNames(idx2)'
Reflectances(idx2,:) = '';
Ref_d1(idx2,:) = '';
Class(idx2) = '';
SampleNames(idx2) = [];
%% Step 3 of using only Log UP for reflectance
idx2 = [634 639 821];
SampleNames(idx2)'
Reflectances(idx2,:) = '';
Ref_d1(idx2,:) = '';
Class(idx2) = '';
SampleNames(idx2) = [];
%% Step 4 of using only Log UP for reflectance
idx2 = [748 984 1001];
SampleNames(idx2)'
Reflectances(idx2,:) = '';
Ref_d1(idx2,:) = '';
Class(idx2) = '';
SampleNames(idx2) = [];
%% Step 5 of using only Log UP for reflectance
idx2 = [592];
SampleNames(idx2)'
Reflectances(idx2,:) = '';
Ref_d1(idx2,:) = '';
Class(idx2) = '';
SampleNames(idx2) = [];
%% Step 6 Log UP for 1st derivative
idx2 = [637 665 904 945 992];
SampleNames(idx2)'
Reflectances(idx2,:) = '';
Ref_d1(idx2,:) = '';
Class(idx2) = '';
SampleNames(idx2) = [];
%% Step 1 For FTIRSmall with Log UP
idx = [4 28 53 62 63 107 192 621 651 653 654 660 731 759 760 877 894 908 909 913 967];
SampleNames(idx)'
Reflectances(idx,:) = '';
Class(idx) = '';
SampleNames(idx) = [];
%% Step 2 For FTIRSmall with Log UP
idx = [108 166 173 428 626 854 924 965 1013];
SampleNames(idx)'
Reflectances(idx,:) = '';
Class(idx) = '';
SampleNames(idx) = [];
%% Step 3 For FTIRSmall with Log UP
idx = [1017 1021];
SampleNames(idx)'
Reflectances(idx,:) = '';
Class(idx) = '';
SampleNames(idx) = [];