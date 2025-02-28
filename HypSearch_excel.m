clear;
clc;
close all;

% List of .mat files containing StoreDice structs
FileDir = 'D:\UCSD_Acads\ProfGal_Research\test_hypsearch';
fileList = dir(fullfile(FileDir, '*.mat'));

% Initialize the Excel file name
excelFileName = fullfile(FileDir, 'HypParamSearch_StoreDiceData.xlsx');

% Loop through each .mat file
NUM_COMBOS_RUN = 60;
NUM_SUBS = 9;
NUM_PAIRS = 3;
%cumu_AvgDice = zeros(1, NUM_COMBOS_RUN);
AllSubjectsDice = zeros(NUM_COMBOS_RUN, NUM_SUBS * NUM_PAIRS);
count_subs = 0;
for i = 1:length(fileList) + 1
    if (i == length(fileList) + 1)
        outputTable = table();
        outputTable.Normalize = data.normalize(:); % (:) is transpose
        outputTable.PCA = data.pca(:);
        outputTable.kNN = data.kNN(:);
        outputTable.SelfTune = data.self_tune(:);
        %overallAvgDice = cumu_AvgDice / length(fileList);  % EACH FILE SHOULD CORRESPOND TO ONE SUBJECT ONLY!!!
        
        overallAvgDice = mean(AllSubjectsDice, 2);
        overallVarDice = var(AllSubjectsDice, 0, 2);
        overallMinDice = min(AllSubjectsDice, [], 2);
        overallMaxDice = max(AllSubjectsDice, [], 2);

        outputTable.OverallAvgDice = overallAvgDice;
        outputTable.overallVarDice = overallVarDice;
        outputTable.overallMinDice = overallMinDice;
        outputTable.overallMaxDice = overallMaxDice;
        sheetName = 'Final';
    
        % Write the table to the Excel file, specifying the sheet name
        writetable(outputTable, excelFileName, 'Sheet', sheetName);
    else
        % Load the .mat file
        matData = load(fullfile(FileDir, fileList(i).name));
        
        % Assume the StoreDice variable is in the loaded data
        % Extract the name of the field (subject name) from the struct
        subjectNames = fieldnames(matData.StoreDice);
        
        % Loop through each subject in the current mat file
        
        for j = 1:length(subjectNames)
            count_subs = count_subs + 1;
            subject = subjectNames{j};
            data = matData.StoreDice.(subject);
    
            % Processing data now
            data.VarDice = getVarianceDice(data.Dice(:));
            [data.MinDice, data.MaxDice] = getMinMaxDice(data.Dice(:));
    
            % Initialize an empty table to store the data
            numEntries = length(data.normalize);  % Assuming all fields have the same length
            outputTable = table();
    
            % Populate the table with data from the struct fields
            outputTable.Normalize = data.normalize(:);
            outputTable.PCA = data.pca(:);
            outputTable.kNN = data.kNN(:);
            outputTable.SelfTune = data.self_tune(:);
            outputTable.AvgDice = data.avgDice(:);
            outputTable.VarDice = data.VarDice(:);
            outputTable.MinDice = data.MinDice(:);
            outputTable.MaxDice = data.MaxDice(:);

            %cumu_AvgDice = cumu_AvgDice + data.avgDice;
            for j2 = 1:numEntries
                AllSubjectsDice(j2, (count_subs-1)*NUM_PAIRS+1 : (count_subs)*NUM_PAIRS) = data.Dice{j2};
            end
            
    
            % Define the sheet name as the subject name to avoid overwriting
            sheetName = subject;
    
            % Write the table to the Excel file, specifying the sheet name
            writetable(outputTable, excelFileName, 'Sheet', sheetName);
    
            % Display message
            fprintf('Data for %s from file %s has been written to %s\n', subject, fileList(i).name, excelFileName);
        end
        clearvars("matData");
    end
end

%% FUNCTIONS
function [minDice, maxDice] = getMinMaxDice(DiceCellArray)
    numentries = length(DiceCellArray);
    minDice = zeros(1, numentries);
    maxDice = zeros(1, numentries);
    for i = 1:numentries
        minDice(i) = min(DiceCellArray{i});
        maxDice(i) = max(DiceCellArray{i});
    end
end

function varDice = getVarianceDice(DiceCellArray)
    numentries = length(DiceCellArray);
    varDice = zeros(1, numentries);
    for i = 1:numentries
        varDice(i) = var(DiceCellArray{i});
    end
end
