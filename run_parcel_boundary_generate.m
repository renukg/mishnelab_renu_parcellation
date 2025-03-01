clear; 
close all;
clc;

dataID = 1;
for dataID = 1:4

    LSSC = 0;
    KMEANS = 1;
    
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
    
    if LSSC
        file_1 = fullfile(lssc_result_dir, 'run_fmri_sessions\sub_SLC09_ses_1_hemisp_1_lssc_out.mat');
        file_2 = fullfile(lssc_result_dir, 'run_fmri_sessions\sub_SLC09_ses_2_hemisp_1_lssc_out.mat');
    end
    if KMEANS
        file_1 = fullfile(kmeans_result_dir, 'run_knn_sessions\sub_SLC09_ses_1_hemisp_1_knn_out.mat');
        file_2 = fullfile(kmeans_result_dir, 'run_knn_sessions\sub_SLC09_ses_2_hemisp_1_knn_out.mat');
    end
    
    data1 = load(file_1);
    data2 = load(file_2);
    
    mergedA1 = data1.mergedA{1};
    mergedA2 = data2.mergedA{1};
    
    bou1 = getBoundariesFromZerosOnes(mergedA1);
    bou2 = getBoundariesFromZerosOnes(mergedA2);
    bou1 = rot90(bou1, -1);
    bou2 = rot90(bou2, -1);
    
    % Initialize options
    options = struct();
    options.thr_method = 'nrg';   % Use threshold-based contour extraction
    options.maxthr = 0.5;         % Threshold for boundary detection
    options.plot_bck_image = false;
    
    % Combine bou1 and bou2 for background
    NROWS = 32;
    NCOLS = 81;
    Cn = zeros(NROWS, NCOLS);  % Background image (combined boundaries)
    
    % Plot contours for bou1
    figure; % Create a new figure
    [CC1, ~, ~] = plot_contours(mergedA1, Cn, options, false, [], [], 'r'); % Red contours for bou1
    
    % Hold the figure to overlay the next plot
    hold on;
    
    % Plot contours for bou2
    [CC2, ~, ~] = plot_contours(mergedA2, Cn, options, false, [], [], 'b'); % Blue contours for bou2
    
    
    % Create dummy plots for legend
    h1 = plot(NaN, NaN, 'r'); % Dummy red line for Session1
    h2 = plot(NaN, NaN, 'b'); % Dummy blue line for Session2
    
    % Customize the plot
    title_str = '';
    if LSSC
        title_str = 'Sub9 - LSSC';
    end
    if KMEANS
        title_str = 'Sub9 - Kmeans';
    end
    title(title_str);
    legend([h1,h2], {'Session1 (Red)', 'Session2 (Blue)'});
    axis equal;
    hold off;
    
    saveimagesPath = 'D:\UCSD_Acads\ProfGal_Research\Notes_ppts_images\images_final_allData';
    if LSSC
        storepath = fullfile(saveimagesPath, 'parcel_boundary_plots');
        if ~exist(storepath, 'dir') 
           mkdir(storepath); % Create the directory
        end
        saveas(gcf, fullfile(storepath, ['data_', dataString, '_Parcel_boundaries_', title_str, '.png']));
    end
    if KMEANS
        storepath = fullfile(saveimagesPath, 'parcel_boundary_plots');
        if ~exist(storepath, 'dir') 
           mkdir(storepath); % Create the directory
        end
        saveas(gcf, fullfile(storepath, ['data_', dataString, '_Parcel_boundaries_', title_str, '.png']));
    end

clear;
close all;
end
