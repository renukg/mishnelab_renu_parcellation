close all;
clear;
clc;

% Directory containing the files
dataDir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\noGSR_thickermask';

% New directory
newdataDir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\noGSR_thickermask_new';

mkdir(newdataDir);
mkdir(fullfile(newdataDir, 'session-1'));
mkdir(fullfile(newdataDir, 'session-2'));
mkdir(fullfile(newdataDir, 'session-3'));

fileList = dir(fullfile(dataDir, '*.nii.gz'));

for i = 1:length(fileList) % processing each .nii file
    fname = fileList(i).name;
    filePath = fullfile(dataDir, fname);

    sub_match = regexp(fname, 'sub-(\w{5})', 'tokens');
    ses_match = regexp(fname, 'ses-(\d+)', 'tokens');
    run_match = regexp(fname, 'run-(\d+)', 'tokens');
    sub_value = sub_match{1}{1};
    ses_value = str2double(ses_match{1}{1});
    run_value = str2double(run_match{1}{1});

    if (ses_value == 1)
        movefile(filePath, fullfile(newdataDir, 'session-1'));
    end
    if (ses_value == 2)
        movefile(filePath, fullfile(newdataDir, 'session-2'));
    end
    if (ses_value == 3)
        movefile(filePath, fullfile(newdataDir, 'session-3'));
    end

    bh=9;
end