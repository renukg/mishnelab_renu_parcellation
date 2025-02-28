clc;
clear; 
close all;

ALLEN = 0;
LSSC = 1;
KMEANS = 0;

dataID = 4;

if dataID == 1
    lssc_result_dir = 'D:\UCSD_Acads\ProfGal_Research\test_run_norm1_pca0_kNN16_sftune4_hemisphere_1_AllenOrdered';
    kmeans_result_dir = 'D:\UCSD_Acads\ProfGal_Research\test_run_KNN25_hemisphere_replicas100_v3_minclstPix_15_new_AllenOrdered';
    allen_result_dir = 'D:\UCSD_Acads\ProfGal_Research\test_run_AllenMaps_v1';
    dataString = 'Old';
    dataDir = 'D:\UCSD_Acads\ProfGal_Research\data32\fMRIData\REST';
end
if dataID == 2
    lssc_result_dir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\GSR_1voxelmask_results\test_run_norm1_pca0_kNN16_sftune4_hem_AllenOrdered';
    kmeans_result_dir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\GSR_1voxelmask_results\test_run_KNN27_hem_replicas100_minclstPix_15_AllenOrdered';
    allen_result_dir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\GSR_1voxelmask_results\test_run_AllenMaps_v1';
    dataString = 'GSR_1voxelmask';
    dataDir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\GSR_1voxelmask_new';
end
if dataID == 3
    lssc_result_dir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\GSR_thickermask_results\test_run_norm1_pca0_kNN16_sftune4_hem_AllenOrdered';
    kmeans_result_dir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\GSR_thickermask_results\test_run_KNN27_hem_replicas100_minclstPix_15_AllenOrdered';
    allen_result_dir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\GSR_thickermask_results\test_run_AllenMaps_v1';
    dataString = 'GSR_thickermask';
    dataDir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\GSR_thickermask_new';
end
if dataID == 4
    lssc_result_dir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\noGSR_thickermask_results\test_run_norm1_pca0_kNN16_sftune4_hem_AllenOrdered';
    kmeans_result_dir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\noGSR_thickermask_results\test_run_KNN27_hem_replicas100_minclstPix_15_AllenOrdered';
    allen_result_dir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\noGSR_thickermask_results\test_run_AllenMaps_v1';
    dataString = 'noGSR_thickermask';
    dataDir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\noGSR_thickermask_new';
end

% All sessions in the data
sesList = dir(fullfile(dataDir, 'session-*'));

if LSSC
    %% LSSC parcellation results directory
    processed_dir = fullfile(lssc_result_dir, 'run_fmri_sessions');
    outDir = fullfile(lssc_result_dir, 'maps_plots');
    if ~exist(outDir, 'dir') 
       mkdir(outDir); % Create the directory
    end
        
    % All processed files
    pFileList = dir(fullfile(processed_dir, '*.mat'));
    
    % Group the files subject-wise
    pFileGrp = struct();
    for i = 1:length(pFileList)
        fname = pFileList(i).name;
        fpath = fullfile(processed_dir, fname);
    
        sub_match = regexp(fname, 'sub_(\w{5})', 'tokens');
        ses_match = regexp(fname, 'ses_(\d+)', 'tokens');
        hem_match = regexp(fname, 'hemisp_(\d+)', 'tokens');
        sub_value = sub_match{1}{1};
        ses_value = str2double(ses_match{1}{1});
        hem_value = str2double(hem_match{1}{1});
        
    
        % Add the file to the subject group
        ses_key = sprintf('ses_%d', ses_value);
        hem_key = sprintf('hem_%d', hem_value);
        if ~isfield(pFileGrp, sub_value)
            pFileGrp.(sub_value) = struct();
        end
        pFileGrp.(sub_value).(ses_key).(hem_key) = fpath;
    end
    
    %% Run LSSC plots images generation
    RUN_LSSC_TIMESERIES = 1;
    RUN_LSSC_SCATTER_PLOT = 0;
    RUN_PARCEL_CORR_COLOR_FIG = 0;
    
    if (RUN_LSSC_TIMESERIES)
        curoutDir = fullfile(outDir, 'timeseries_plots');
        if ~exist(curoutDir, 'dir') 
           mkdir(curoutDir); % Create the directory
        end
        for ses = 1:length(sesList) % each session of data
            sesname = sesList(ses).name;
            sesDir = [dataDir, '\', sesname];
            fileList = dir(fullfile(sesDir, '*.nii.gz'));
            
            V = []; % To hold the time series data - 3D matrix
            prev_sub = '';
            prev_ses = '';
            for i = 1:length(fileList) % processing each .nii file
                fname = fileList(i).name;
                filePath = fullfile(sesDir, fname);
            
                sub_match = regexp(fname, 'sub-(\w{5})', 'tokens');
                ses_match = regexp(fname, 'ses-(\d+)', 'tokens');
                run_match = regexp(fname, 'run-(\d+)', 'tokens');
                sub_value = sub_match{1}{1};
                ses_value = str2double(ses_match{1}{1});
                run_value = str2double(run_match{1}{1});
                
                V0 = niftiread(filePath);
                V0 = squeeze(V0);
            
                if ((i == 1 || strcmp(prev_sub, sub_value)) && i ~= length(fileList))    % runs belong to same subject
                    if isempty(V)
                        V = V0;
                    else
                        V = cat(3, V, V0); % concatenate runs
                    end 
                else   % All runs concatenated
                    if (i == length(fileList)) % if last file
                        prev_sub = sub_value;
                        prev_ses = ses_value;
                    end
    
                    MIDLINE1 = 29;
                    MIDLINE2 = 32;
                    % Separating the brain region into two halves along the midline
                    % Extract the top and bottom halves
                    V_top = V(1:MIDLINE1-1, :, :);        % Top half
                    V_bottom = V(MIDLINE2+1:end, :, :);  % Bottom half
                    
                    % Flip the bottom half vertically (along the row axis)
                    V_bottom_flipped = flip(V_bottom, 1);  % Flip along the 1st dimension (rows)
    
                    num_time_frames = size(V_top, 3);
                    
                    % Pre-allocate a new matrix to hold the interleaved data
                    ROW_PAD = 4; % To avoid out of index error in LSSC algorithm
                    V_hem = zeros(size(V_top, 1)+ROW_PAD, size(V_top, 2), num_time_frames);
    
                    % Run for both hemispheres SEPARATELY
                    for hemisp = 1:2
                        for t = 1:num_time_frames
                            if (hemisp == 1)
                                V_hem(1:MIDLINE1-1, :, t) = V_top(:, :, t);
                            else
                                V_hem(1:MIDLINE1-1, :, t) = V_bottom_flipped(:, :, t);
                            end
                        end
                        
                        % Retrieve LSSC parcel map from pre-saved data
                        sub_k = prev_sub;
                        ses_k = sprintf('ses_%d', prev_ses);
                        hem_k = sprintf('hem_%d', hemisp);
                        cur_lssc_data = load(pFileGrp.(sub_k).(ses_k).(hem_k));
                        lssc_map_hem = cur_lssc_data.labels{1};
                        lssc_map_hem = rot90(lssc_map_hem, 1);
    
                        % Create brain mask for the hemispheres
                        maxV = max(V_hem, [], 3);
                        brain_mask = maxV ~= 0;
                        [R,C] = size(brain_mask);
        
                        % Flatten the spatial dimensions
                        V_flat = reshape(V_hem, [], size(V_hem, 3));
                        
                        % Extract relevant pixels
                        allregionspix = find(brain_mask);
                        dFoF_masked = V_flat(allregionspix, :);
    
                        lssc_map_hem = lssc_map_hem.*brain_mask;
                        clusterIDs = unique(lssc_map_hem(:)); % cluster IDs
                        clusterIDs(clusterIDs == 0) = [];
                        num_clusters = length(clusterIDs);
                        
                        % plot time series
                        parcel_response = zeros(num_clusters, num_time_frames);
                        % Iterate over each cluster ID
                        for id = 1:num_clusters
                            clusterID = clusterIDs(id);
                            
                            % Find the mask for the current cluster
                            clusterMask = (lssc_map_hem == clusterID);
                            
                            % cluster pixels
                            clus_pix_ids = find(clusterMask);
            
                            % Extract the time series data for the pixels in this cluster
                            clusterTimeSeries = V_hem(repmat(clusterMask, [1, 1, num_time_frames]));
                            
                            % Reshape the time series data into a 2D matrix (pixels x time)
                            clusterTimeSeries = reshape(clusterTimeSeries, [], num_time_frames);
                            
                            parcel_response(id, :) = mean(clusterTimeSeries);
                        end
                        % Define an offset for stacking the time series
                        offset = 80;  % Adjust based on data scale
                        
                        % Create a figure
                        title_str = sprintf("LSSC Parcel Response: Sub %s, Sess %d, Hem %d", prev_sub, prev_ses, hemisp);
                        figure;
                        hold on;
                        
                        % Plot each time series with an offset
                        for j = 1:num_clusters
                            plot(parcel_response(j, :) + (j-1) * offset, 'DisplayName', sprintf('p%d', clusterIDs(j)));
                        end
                        yticks((0:num_clusters-1) * offset);
                        yticklabels(arrayfun(@(x) sprintf('p%d', clusterIDs(x)), 1:num_clusters, 'UniformOutput', false));
                        set(gca, 'YDir', 'reverse');  % Invert y-axis
                        
                        xlabel('Time');
                        title(title_str);
                        %legend('show', 'Location', 'northeastoutside');  % Optional legend outside the plot
                        hold off;
                        saveas(gcf, fullfile(curoutDir,['LSSC_Parcel_Response_',prev_sub,...
                            '_sess_',num2str(prev_ses), '_hem_', num2str(hemisp),'.png']));
                        close;
    
                        % Pairwise correlation
                        pairwise_Correlation = zeros(num_clusters);
                        for j1 = 1:num_clusters
                            for j2 = j1:num_clusters
                                clust1 = clusterIDs(j1);
                                clust2 = clusterIDs(j2);
                            
                                % Find the masks
                                clustMask1 = (lssc_map_hem == clust1);
                                clustMask2 = (lssc_map_hem == clust2);
                                
                                % cluster pixels
                                clus_pix_ids1 = find(clustMask1);
                                clus_pix_ids2 = find(clustMask2);
                
                                % Extract the time series data for the pixels in this cluster
                                clusterTimeSeries1 = V_hem(repmat(clustMask1, [1, 1, num_time_frames]));
                                clusterTimeSeries2 = V_hem(repmat(clustMask2, [1, 1, num_time_frames]));
                                
                                % Reshape the time series data into a 2D matrix (pixels x time)
                                clusterTimeSeries1 = reshape(clusterTimeSeries1, [], num_time_frames);
                                clusterTimeSeries2 = reshape(clusterTimeSeries2, [], num_time_frames);
    
                                cross_corr_ = corr(clusterTimeSeries1', clusterTimeSeries2');
                                pairwise_Correlation(j1, j2) = mean(cross_corr_(:));
                                pairwise_Correlation(j2, j1) = pairwise_Correlation(j1, j2);
                            end
                        end
    
                        title_str = sprintf("LSSC Temp correlation: Sub %s, Sess %d, Hem %d", prev_sub, prev_ses, hemisp);
                        figure;
                        imagesc(pairwise_Correlation);
                        colorbar;
                        title(title_str);
                        saveas(gcf, fullfile(curoutDir,['LSSC_TempCorr_pairwise_',prev_sub,...
                            '_sess_',num2str(prev_ses), '_hem_', num2str(hemisp),'.png']));
                        close;
    
                    end
        
                    V = V0; % initialize to the new subject
                end
                prev_sub = sub_value;
                prev_ses = ses_value;
            end
        end
    end

    if RUN_LSSC_SCATTER_PLOT
        curoutDir = fullfile(outDir, 'scatter_plots');
        rFileList = dir(fullfile(resultDir, '*.mat'));
        for i = 1:length(rFileList)
            fname = rFileList(i).name;
            fpath = fullfile(resultDir, fname);
            
            sub_match = regexp(fname, 'sub_(\w{5})', 'tokens');
            ses_match = regexp(fname, 'session_(\d+)', 'tokens');
            hem_match = regexp(fname, 'hemisp_(\d+)', 'tokens');
            sub_value = sub_match{1}{1};
            ses_value = str2double(ses_match{1}{1});
            hemisp = str2double(hem_match{1}{1});

            cur_data_forcorr = load(fpath);
            ses_k = sprintf('ses_%d', ses_value);
            hem_k = sprintf('hem_%d', hemisp);
            
            NROWS = cur_data_forcorr.NROWS;
            NCOLS = cur_data_forcorr.NCOLS;
            mergedA = cur_data_forcorr.mergedA;
            labels = zeros(NROWS,NCOLS);
            for i2= 1:size(mergedA,2)
                labels(mergedA(:,i2)>0) = i2;
            end
            cluster_sizes = sum(mergedA, 1);
            within_clust_corr = cur_data_forcorr.clusterwise_within_corr;
            figure;
            scatter(cluster_sizes, within_clust_corr, 20, 'b', 'filled', 'MarkerEdgeColor', 'k'); 
            xlabel('Cluster Size');
            ylabel('Within-Cluster Correlation');
            title('Cluster Size vs. Within-Cluster Correlation');
            saveas(gcf, fullfile(curoutDir,['LSSC_ClustSize_vs_Corr_',sub_value,...
                            '_sess_',num2str(ses_value), '_hem_', num2str(hemisp),'.png']));
            close;
        end
    end

    % Plot where each parcel is colored by its average correlation
    if RUN_PARCEL_CORR_COLOR_FIG
        curoutDir = fullfile(outDir, 'parcel_corr_color_plots');
        if ~exist(curoutDir, 'dir') 
            mkdir(curoutDir); % Create the directory
        end
        rFileList = dir(fullfile(resultDir, '*.mat'));
        for i = 1:length(rFileList)
            fname = rFileList(i).name;
            fpath = fullfile(resultDir, fname);
            
            sub_match = regexp(fname, 'sub_(\w{5})', 'tokens');
            ses_match = regexp(fname, 'session_(\d+)', 'tokens');
            hem_match = regexp(fname, 'hemisp_(\d+)', 'tokens');
            sub_value = sub_match{1}{1};
            ses_value = str2double(ses_match{1}{1});
            hemisp = str2double(hem_match{1}{1});

            cur_data_forcorr = load(fpath);
            ses_k = sprintf('ses_%d', ses_value);
            hem_k = sprintf('hem_%d', hemisp);
            
            NROWS = cur_data_forcorr.NROWS;
            NCOLS = cur_data_forcorr.NCOLS;
            mergedA = cur_data_forcorr.mergedA;
            labels = zeros(NROWS,NCOLS);
            for i2= 1:size(mergedA,2)
                labels(mergedA(:,i2)>0) = i2;
            end
            cluster_sizes = sum(mergedA, 1);
            within_clust_corr = cur_data_forcorr.clusterwise_within_corr;
            dispCorrColor = zeros(size(labels));
            parcels = unique(labels);
            parcels(parcels == 0) = [];
            for i3 = 1:length(parcels)
                dispCorrColor(labels == i3) = within_clust_corr(i3);
            end

            figure;
            imagesc(dispCorrColor);
            colormap(jet);
            colorbar;
            clim([min(within_clust_corr) max(within_clust_corr)]); % Set color limits
            axis off; % Hide axis ticks (optional)
            title('Within-Cluster Correlation Map');
            saveas(gcf, fullfile(curoutDir,['LSSC_ParcelCorr_color_',sub_value,...
                            '_sess_',num2str(ses_value), '_hem_', num2str(hemisp),'.png']));
            close;
        end
    end
end

%%
if KMEANS
    processed_dir = fullfile(kmeans_result_dir, 'run_knn_sessions');
    outDir = fullfile(kmeans_result_dir, 'maps_plots');
    if ~exist(outDir, 'dir') 
       mkdir(outDir); % Create the directory
    end

    RUN_KNN_SCATTER_PLOT = 0;
    RUN_PARCEL_CORR_COLOR_FIG = 1;

    if RUN_KNN_SCATTER_PLOT
        curoutDir = fullfile(outDir, 'scatter_plots');
        pFileList = dir(fullfile(processed_dir, '*.mat'));
        for i = 1:length(pFileList)
            fname = pFileList(i).name;
            fpath = fullfile(processed_dir, fname);
            
            sub_match = regexp(fname, 'sub_(\w{5})', 'tokens');
            ses_match = regexp(fname, 'ses_(\d+)', 'tokens');
            hem_match = regexp(fname, 'hemisp_(\d+)', 'tokens');
            sub_value = sub_match{1}{1};
            ses_value = str2double(ses_match{1}{1});
            hemisp = str2double(hem_match{1}{1});

            cur_data_forcorr = load(fpath);
            ses_k = sprintf('ses_%d', ses_value);
            hem_k = sprintf('hem_%d', hemisp);
            
            labels = cur_data_forcorr.labels{1};
            NROWS = size(labels, 1);
            NCOLS = size(labels, 2);
            mergedA = cur_data_forcorr.mergedA{1};
            cluster_sizes = sum(mergedA, 1);
            within_clust_corr = cur_data_forcorr.clusterwise_within_corr;
            figure;
            scatter(cluster_sizes, within_clust_corr, 20, 'b', 'filled', 'MarkerEdgeColor', 'k'); 
            xlabel('Cluster Size');
            ylabel('Within-Cluster Correlation');
            title('Cluster Size vs. Within-Cluster Correlation');
            saveas(gcf, fullfile(curoutDir,['Kmeans_ClustSize_vs_Corr_',sub_value,...
                            '_sess_',num2str(ses_value), '_hem_', num2str(hemisp),'.png']));
            close;
        end
    end

    % Plot where each parcel is colored by its average correlation
    if RUN_PARCEL_CORR_COLOR_FIG
        curoutDir = fullfile(outDir, 'parcel_corr_color_plots');
        if ~exist(curoutDir, 'dir') 
            mkdir(curoutDir); % Create the directory
        end
        pFileList = dir(fullfile(processed_dir, '*.mat'));
        for i = 1:length(pFileList)
            fname = pFileList(i).name;
            fpath = fullfile(processed_dir, fname);
            
            sub_match = regexp(fname, 'sub_(\w{5})', 'tokens');
            ses_match = regexp(fname, 'ses_(\d+)', 'tokens');
            hem_match = regexp(fname, 'hemisp_(\d+)', 'tokens');
            sub_value = sub_match{1}{1};
            ses_value = str2double(ses_match{1}{1});
            hemisp = str2double(hem_match{1}{1});

            cur_data_forcorr = load(fpath);
            ses_k = sprintf('ses_%d', ses_value);
            hem_k = sprintf('hem_%d', hemisp);
            
            labels = cur_data_forcorr.labels{1};
            NROWS = size(labels, 1);
            NCOLS = size(labels, 2);
            mergedA = cur_data_forcorr.mergedA{1};
            cluster_sizes = sum(mergedA, 1);
            within_clust_corr = cur_data_forcorr.clusterwise_within_corr;
            dispCorrColor = zeros(size(labels));
            parcels = unique(labels);
            parcels(parcels == 0) = [];
            for i3 = 1:length(parcels)
                dispCorrColor(labels == i3) = within_clust_corr(i3);
            end
            %dispCorrColor = dispCorrColor.*10;
            figure;
            imagesc(dispCorrColor);
            colormap(jet);
            colorbar;
            clim([min(within_clust_corr) max(within_clust_corr)]); % Set color limits
            axis off; 
            title('Within-Cluster Correlation Map');
            saveas(gcf, fullfile(curoutDir,['Kmeans_ParcelCorr_color_',sub_value,...
                            '_sess_',num2str(ses_value), '_hem_', num2str(hemisp),'.png']));
            close;
        end
    end
end
%%
if ALLEN
    processed_dir = fullfile(allen_result_dir, 'run_allen_sessions');
    outDir = fullfile(allen_result_dir, 'maps_plots');
    if ~exist(outDir, 'dir') 
       mkdir(outDir); % Create the directory
    end
    RUN_ALLEN_SCATTER_PLOT = 1;

    if RUN_ALLEN_SCATTER_PLOT
        curoutDir = fullfile(outDir, 'scatter_plots');
        pFileList = dir(fullfile(processed_dir, '*.mat'));
        for i = 1:length(pFileList)
            fname = pFileList(i).name;
            fpath = fullfile(processed_dir, fname);
            
            sub_match = regexp(fname, 'sub_(\w{5})', 'tokens');
            ses_match = regexp(fname, 'ses_(\d+)', 'tokens');
            hem_match = regexp(fname, 'hemisp_(\d+)', 'tokens');
            sub_value = sub_match{1}{1};
            ses_value = str2double(ses_match{1}{1});
            hemisp = str2double(hem_match{1}{1});

            cur_data_forcorr = load(fpath);
            ses_k = sprintf('ses_%d', ses_value);
            hem_k = sprintf('hem_%d', hemisp);
            
            labels = cur_data_forcorr.labels;
            NROWS = size(labels, 1);
            NCOLS = size(labels, 2);
            mergedA = cur_data_forcorr.mergedA;
            cluster_sizes = sum(mergedA, 1);
            within_clust_corr = cur_data_forcorr.clusterwise_within_corr;
            figure;
            scatter(cluster_sizes, within_clust_corr, 20, 'b', 'filled', 'MarkerEdgeColor', 'k'); 
            xlabel('Cluster Size');
            ylabel('Within-Cluster Correlation');
            title('Cluster Size vs. Within-Cluster Correlation');
            saveas(gcf, fullfile(curoutDir,['Allen_ClustSize_vs_Corr_',sub_value,...
                            '_sess_',num2str(ses_value), '_hem_', num2str(hemisp),'.png']));
            close;
        end
    end
end