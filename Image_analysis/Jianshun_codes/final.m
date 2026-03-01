%% Main Script to Run Analysis on Multiple Datasets
% This script processes multiple experimental runs by calling runAnalyse.
% Uncomment the desired sections or modify the loop to process specific datasets.

clear; clc;

% Define a list of datasets to process
% Each entry: {sequence, date, run, folderPath, ROIcenter, measurementName}
% datasets = {'TwoDGas', '2026/02/04', '0000', 'C:\Users\Joschka\Desktop\Remy\Codes\SR_droplets\Image_analysis\Jianshun_codes\Data_analysed\', ...
%         [1490 2080], 'Stripe90'};

    % Example: 'StructuralPhaseTransition', '2025/10/31', '0002', ...
    % 'C:\Users\Jianshun Gao\Documents\coding\Calculations\Data\2025\10\31\0002\',
    % [1490 2080], 'BECToDroplets'
    % };

% Uncomment and add datasets as needed:
% datasets = [datasets];
% ; {'TwoDGas', '2026/02/04', '0000', ...
%     'C:\Users\Joschka\Desktop\Remy\Codes\SR_droplets\Image_analysis\Jianshun_codes\Data_analysed\', ...
%     [], 'Default'}];

% Alternatively, process a single dataset by setting parameters manually
year = "2026";
month = "02";
day = "11";
shot = "0000";

dataSources = {struct('sequence', 'StructuralPhaseTransition', ...
                      'date', char(year + "/" + month + "/" + day), ...
                      'runs', shot)};
folderPath = "C:\Users\Joschka\Desktop\Remy\Codes\SR_droplets\Image_analysis\Jianshun_codes\Data_analysed\" + year + "\" + month + "\" + day + "\" + shot + "\";
ROIcenter = [1500 1950];
measurementName = "Stripes_90°";

% Call the analysis function
runAnalyse(dataSources, folderPath, ROIcenter, measurementName);

% After processing, optionally clean up temporary folders
% deleteFullODImagesFoldersQuick("C:\Users\Jianshun Gao\Documents\DyData");

%% Batch processing (example)
% for i = 1:size(datasets,1)
%     ds = datasets(i,:);
%     dataSources = {struct('sequence', ds{1}, 'date', ds{2}, 'runs', ds{3})};
%     folderPath = ds{4};
%     ROIcenter = ds{5};
%     measurementName = ds{6};
%     runAnalyse(dataSources, folderPath, ROIcenter, measurementName);
%     deleteFullODImagesFoldersQuick("C:\Users\Jianshun Gao\Documents\DyData");
% end