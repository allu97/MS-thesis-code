clearvars;clc
hn = replaceNaNs('Hailniemi.mat','Hailniemi',3);
hk = replaceNaNs('Hameenkyla.mat','Hameenkyla',3);
hepo = replaceNaNs('HeponiemiREE.mat','Heponiemi_raw',3);
s1 = replaceNaNs('Sample1.mat','Sample1',3);
s2 = replaceNaNs('Sample2.mat','Sample2',3);

%% GeoChemREE data
hn = replaceNaNs2('Hailniemi.mat','Hailniemi_ree');
hk = replaceNaNs2('Hameenkyla.mat','Hameenkyla_ree');
hepo = replaceNaNs2('HeponiemiREE2.mat','Heponiemi_raw');
s1 = replaceNaNs2('Sample1.mat','Sample1_ree');
s2 = replaceNaNs2('Sample2.mat','Sample2_ree');
%% GeoChem data
hn = replaceNaNs2('Hailniemi.mat','Hailniemi_g');
hk = replaceNaNs2('Hameenkyla.mat','Hameenkyla_g');
hepo = replaceNaNs2('Heponiemi2.mat','Hepo_');
s1 = replaceNaNs2('Sample1.mat','Sample1_g');
s2 = replaceNaNs2('Sample2.mat','Sample2_g');
function X = replaceNaNs(filename,tablename,denominator)
% This function replaces NaN with uncertainty value found or "Error"
% NaN is replaced with uncertainty divided by denominator
File = load(filename,tablename);

% Data = File.tablename;
Data = eval(strcat('File.',tablename));
[n,m] = size(Data);

X = zeros(n-1,ceil(m/2));

for i = 1:n
    for j = 2:2:m-1
        x = Data(i,j).Variables;
        if isnan(x)
            x = Data(i,j+1).Variables/denominator;
        end
        X(i,j/2) = x;
    end
end
end

function X = replaceNaNs2(filename,tablename)
% This function replaces NaN values with error values obtained from
% variables that are still measured.
% The replaced value is 3*mean of error of obtained values.
% This is the same value as 99% LOD for the sample variable??
File = load(filename,tablename);
Data = eval(strcat('File.',tablename));
values = Data(:,2:end).Variables;
[n,m] = size(values);

X = zeros(n,ceil(m/2));
% These need to be replace, so get the indeces
nans = isnan(values);
% Get the indeces of values still present in NaN columns
not_nans = ~nans;
not_nans(:,~any(nans)) = false;
% Shift each row one position to right to get error idx
nans_err = circshift(not_nans,1,2);
% Get mean of error measurements for each column 
err_values = zeros(n,m);
err_values(nans_err) = values(nans_err);
% mean_err = 3*(sum(err_values)./sum(nans_err));
mean_err = ceil(3*(sum(err_values)./sum(nans_err)));


for i = 1:n
    for j = 2:2:m
        x = values(i,j-1);
        if isnan(x)
            % Replace NaN with mean*3 of error values
            x = mean_err(j);
        end
        X(i,j/2) = x;
    end
end
end