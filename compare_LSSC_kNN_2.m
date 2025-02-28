clc;
clear;
close all;

dataID = 4;

if dataID == 1
    lssc_result_dir = 'D:\UCSD_Acads\ProfGal_Research\test_run_norm1_pca0_kNN16_sftune4_hemisphere_1_AllenOrdered';
    kmeans_result_dir = 'D:\UCSD_Acads\ProfGal_Research\test_run_KNN25_hemisphere_replicas100_v3_minclstPix_15_new_AllenOrdered';
    allen_result_dir = 'D:\UCSD_Acads\ProfGal_Research\test_run_AllenMaps_v1';
    dataString = 'Old';
end
if dataID == 2
    lssc_result_dir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\GSR_1voxelmask_results\test_run_norm1_pca0_kNN16_sftune4_hem_AllenOrdered';
    kmeans_result_dir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\GSR_1voxelmask_results\test_run_KNN27_hem_replicas100_minclstPix_15_AllenOrdered';
    allen_result_dir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\GSR_1voxelmask_results\test_run_AllenMaps_v1';
    dataString = 'GSR_1voxelmask';
end
if dataID == 3
    lssc_result_dir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\GSR_thickermask_results\test_run_norm1_pca0_kNN16_sftune4_hem_AllenOrdered';
    kmeans_result_dir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\GSR_thickermask_results\test_run_KNN27_hem_replicas100_minclstPix_15_AllenOrdered';
    allen_result_dir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\GSR_thickermask_results\test_run_AllenMaps_v1';
    dataString = 'GSR_thickermask';
end
if dataID == 4
    lssc_result_dir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\noGSR_thickermask_results\test_run_norm1_pca0_kNN16_sftune4_hem_AllenOrdered';
    kmeans_result_dir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\noGSR_thickermask_results\test_run_KNN27_hem_replicas100_minclstPix_15_AllenOrdered';
    allen_result_dir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\noGSR_thickermask_results\test_run_AllenMaps_v1';
    dataString = 'noGSR_thickermask';
end

saveimagesPath = 'D:\UCSD_Acads\ProfGal_Research\Notes_ppts_images\images_final_allData\Allen_KMeans_LSSC_Dice_Corr_comparison';

%%
matfilepath = fullfile(lssc_result_dir, 'stats_mat_files');
dice_LSSC = load(fullfile(matfilepath, "lssc_avg_dice.mat"));
mean_lssc_dice = mean(dice_LSSC.avg_dice_subwise);
std_lssc_dice = std(dice_LSSC.avg_dice_subwise);


matfilepath1 = fullfile(kmeans_result_dir, 'stats_mat_files');
dice_KNN = load(fullfile(matfilepath1, "knn_avg_dice.mat"));
mean_knn_dice = mean(dice_KNN.avg_dice_subwise);
std_knn_dice = std(dice_KNN.avg_dice_subwise);

% Perform t-test for Dice Coefficients
[~, p_dice] = ttest2(dice_LSSC.avg_dice_subwise, dice_KNN.avg_dice_subwise);

%% Dice comparison
data = [mean_lssc_dice, mean_knn_dice];
errors = [std_lssc_dice, std_knn_dice];
% Create the bar plot
figure;
hold on;
bar_handle = bar(data, 'FaceColor', [0.7, 0.7, 0.7]);  % Set bar color to gray

% Add error bars
errorbar(1:length(data), data, errors, 'k.', 'LineWidth', 1.5);  % 'k.' for black error bars

% Customize x-axis labels
xticks(1:length(data));
xticklabels({'LSSC', 'kNN'});
ylabel('Dice Coefficient');

% Customize the plot appearance
ylim([0, 1]);  % Set y-axis limits
title('Dice comparison'); 
hold off;
title_str = 'Dice comparison';

% Add significance marker if p-value < 0.05
hold on;
if p_dice < 0.05
    y = max(data + errors) + 0.05;  % Position for significance line  
    plot([1, 2], [y, y], 'k-', 'LineWidth', 1.5);  % Horizontal line
    text(1.5, y + 0.02, '*', 'FontSize', 14, 'HorizontalAlignment', 'center');  % Asterisk
end
hold off;

saveas(gcf, fullfile(saveimagesPath,['Data_',dataString, '_', title_str,'.png']));

%% Correlation
matfilepath2 = fullfile(allen_result_dir, 'stats_mat_files');
allen_corr_dat = load(fullfile(matfilepath2, "Allen_corr_means.mat"));

lssc_corr_dat = load(fullfile(matfilepath, "lssc_corr_means.mat"));
knn_corr_dat = load(fullfile(matfilepath1, "KNN_corr_means.mat"));
for cc = 1:2
    if cc == 1
        allen_within = allen_corr_dat.mean_corr_within_subwise;
        lssc_within = lssc_corr_dat.mean_corr_within_subwise;
        knn_within = knn_corr_dat.mean_corr_within_subwise;

        data = [mean(allen_within), mean(lssc_within), mean(knn_within)];
        errors = [std(allen_within), std(lssc_within), std(knn_within)];

        % Perform t-tests between groups
        [~, p_allen_lssc] = ttest2(allen_within, lssc_within);
        [~, p_lssc_knn] = ttest2(lssc_within, knn_within);
        [~, p_allen_knn] = ttest2(allen_within, knn_within);
    else
        allen_cross = allen_corr_dat.mean_corr_across_subwise;
        lssc_cross = lssc_corr_dat.mean_corr_across_subwise;
        knn_cross = knn_corr_dat.mean_corr_across_subwise;

        data = [mean(allen_cross), mean(lssc_cross), mean(knn_cross)];
        errors = [std(allen_cross), std(lssc_cross), std(knn_cross)];

        % Perform t-tests between groups
        [~, p_allen_lssc] = ttest2(allen_cross, lssc_cross);
        [~, p_lssc_knn] = ttest2(lssc_cross, knn_cross);
        [~, p_allen_knn] = ttest2(allen_cross, knn_cross);
    end
    
    % Create the bar plot
    figure;
    hold on;
    bar_handle = bar(data, 'FaceColor', [0.7, 0.7, 0.7]);  % Set bar color to gray
    
    % Add error bars
    errorbar(1:length(data), data, errors, 'k.', 'LineWidth', 1.5);  % 'k.' for black error bars
    
    % Customize x-axis labels
    xticks(1:length(data));
    xticklabels({'Anatomical', 'LSSC', 'kNN'});
    ylabel('Correlation');
    
    % Customize the plot appearance
    ylim([0, 0.8]);  % Set y-axis limits
    if cc == 1
        title('Correlation within'); 
        title_str = 'Correlation within';
    end
    if cc == 2
        title('Correlation across');  
        title_str = 'Correlation across';
    end

    % Add significance markers
    y_offset = 0.02;  % Offset for significance lines

    % Allen vs. LSSC
    if p_allen_lssc < 0.05
        y1 = max(data(1:2) + errors(1:2)) + y_offset;
        plot([1, 2], [y1, y1], 'k-', 'LineWidth', 1.5);
        text(1.5, y1 + 0.01, '*', 'FontSize', 14, 'HorizontalAlignment', 'center');
    end

    % LSSC vs. kNN
    if p_lssc_knn < 0.05
        y2 = max(data(2:3) + errors(2:3)) + y_offset;
        plot([2, 3], [y2, y2], 'k-', 'LineWidth', 1.5);
        text(2.5, y2 + 0.01, '*', 'FontSize', 14, 'HorizontalAlignment', 'center');
    end

    % Allen vs. kNN
    if p_allen_knn < 0.05
        y3 = max(data + errors) + 0.05;  % Higher line for Allen vs. kNN
        plot([1, 3], [y3, y3], 'k-', 'LineWidth', 1.5);
        text(2, y3 + 0.01, '*', 'FontSize', 14, 'HorizontalAlignment', 'center');
    end

    hold off;
    saveas(gcf, fullfile(saveimagesPath,['Data_',dataString, '_', title_str,'.png']));
end