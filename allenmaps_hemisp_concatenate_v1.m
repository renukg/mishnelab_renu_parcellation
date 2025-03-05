clc;
clear; 
close all;
%% allenmaps_v1.m is for hemispheres. allenmaps_v0.m is for full brain images
% Plots time series, pairwise correlations, temporal correlation - within and across

% Maps directory
mapDir = 'D:\UCSD_Acads\ProfGal_Research\Allen maps';

dataID = 2;

if dataID == 1
    dataDir = 'D:\UCSD_Acads\ProfGal_Research\data32\fMRIData\REST';
    resultDir = 'D:\UCSD_Acads\ProfGal_Research\test_run_AllenMaps_hem_concatenate_v1';
end
if dataID == 2
    dataDir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\GSR_1voxelmask_new';
    resultDir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\GSR_1voxelmask_results\test_run_AllenMaps_hem_concatenate_v1';
end
if dataID == 3
    dataDir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\GSR_thickermask_new';
    resultDir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\GSR_thickermask_results\test_run_AllenMaps_hem_concatenate_v1';
end
if dataID == 4
    dataDir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\noGSR_thickermask_new';
    resultDir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\noGSR_thickermask_results\test_run_AllenMaps_hem_concatenate_v1';
end

AllenfilePath = fullfile(mapDir, '2D_calcium_atlas.nii');

allen_map = niftiread(AllenfilePath);

% All sessions in the directory
sesList = dir(fullfile(dataDir, 'session-*'));

AllenCorrStorePath = fullfile(resultDir, 'run_allen_sessions');
RUN_ALLEN_PROCESSING = 0;
RUN_ALLEN_REPORTING = 1;

if (RUN_ALLEN_PROCESSING)
    for ses = 1:length(sesList) % processing each session
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

                % Resize the Allen map to the size of the brain image using nearest-neighbor interpolation
                allenMapResized = imresize(allen_map, [size(V, 1), size(V, 2)], 'nearest');

                MIDLINE1 = 29;
                MIDLINE2 = 32;
                % Separating the brain region into two halves along the midline
                % Extract the top and bottom halves
                V_top = V(1:MIDLINE1-1, :, :);        % Top half
                V_bottom = V(MIDLINE2+1:end, :, :);  % Bottom half

                allenMap_top = allenMapResized(1:MIDLINE1-1, :);
                allenMap_bottom = allenMapResized(MIDLINE2+1:end, :, :);
                
                % Flip the bottom half vertically (along the row axis)
                V_bottom_flipped = flip(V_bottom, 1);  % Flip along the 1st dimension (rows)
                allenMap_bottom_flipped = flip(allenMap_bottom, 1);

                num_time_frames = size(V_top, 3);
                
                % Pre-allocate a new matrix to hold the interleaved data
                ROW_PAD = 4; % To avoid out of index error in LSSC algorithm
                V_hem = zeros(size(V_top, 1)+ROW_PAD, size(V_top, 2), 2*num_time_frames);
                allenMap_hem = zeros(size(V_top, 1)+ROW_PAD, size(V_top, 2));
                
                %Interleave the top and bottom halves in time
                for t = 1:num_time_frames
                    V_hem(1:MIDLINE1-1, :, 2*t-1) = V_top(:, :, t);  % Odd indices: V_top
                    V_hem(1:MIDLINE1-1, :, 2*t) = V_bottom_flipped(:, :, t);  % Even indices: V_bottom_flipped
                end

                allenMap_hem(1:MIDLINE1-1, :) = allenMap_top;  % use the top half for now. @renu - revisit

                % Create brain mask for the hemispheres
                maxV = max(V_hem, [], 3);
                brain_mask = maxV ~= 0;
                [R,C] = size(brain_mask);

                % Flatten the spatial dimensions
                V_flat = reshape(V_hem, [], size(V_hem, 3));
                
                % Extract relevant pixels
                allregionspix = find(brain_mask);
                dFoF_masked = V_flat(allregionspix, :);
                
                % Plot allen maps in brain region (for ppt)
                allenMap_hem = allenMap_hem.*brain_mask;
                %figure;imagesc(label2rgb(rot90(allenMap_hem, -1))); % for brain region with clusters

                clusterIDs = unique(allenMap_hem(:)); % cluster IDs
                clusterIDs(clusterIDs == 0) = [];
                num_clusters = length(clusterIDs);
                
                % plot time series
                parcel_response = zeros(num_clusters, 2*num_time_frames);
                % Iterate over each cluster ID
                for id = 1:num_clusters
                    clusterID = clusterIDs(id);
                    
                    % Find the mask for the current cluster
                    clusterMask = (allenMap_hem == clusterID);
                    
                    % cluster pixels
                    clus_pix_ids = find(clusterMask);
    
                    % Extract the time series data for the pixels in this cluster
                    clusterTimeSeries = V_hem(repmat(clusterMask, [1, 1, 2*num_time_frames]));
                    
                    % Reshape the time series data into a 2D matrix (pixels x time)
                    clusterTimeSeries = reshape(clusterTimeSeries, [], 2*num_time_frames);
                    
                    parcel_response(id, :) = mean(clusterTimeSeries);
                end
                % Define an offset for stacking the time series
                offset = 80;  % Adjust based on data scale
                
                % Create a figure
                title_str = sprintf("Allen Parcel Response: Sub %s, Sess %d", prev_sub, prev_ses);
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
                hold off;
                saveas(gcf, fullfile(resultDir,['Allen_Parcel_Response_',prev_sub,...
                    '_sess_',num2str(prev_ses), '.png']));
                close;

                % Pairwise correlation
                
                for j1 = 1:num_clusters
                    for j2 = j1:num_clusters
                        clust1 = clusterIDs(j1);
                        clust2 = clusterIDs(j2);
                    
                        % Find the masks
                        clustMask1 = (allenMap_hem == clust1);
                        clustMask2 = (allenMap_hem == clust2);
                        
                        % cluster pixels
                        clus_pix_ids1 = find(clustMask1);
                        clus_pix_ids2 = find(clustMask2);
        
                        % Extract the time series data for the pixels in this cluster
                        clusterTimeSeries1 = V_hem(repmat(clustMask1, [1, 1, 2*num_time_frames]));
                        clusterTimeSeries2 = V_hem(repmat(clustMask2, [1, 1, 2*num_time_frames]));
                        
                        % Reshape the time series data into a 2D matrix (pixels x time)
                        clusterTimeSeries1 = reshape(clusterTimeSeries1, [], 2*num_time_frames);
                        clusterTimeSeries2 = reshape(clusterTimeSeries2, [], 2*num_time_frames);

                        cross_corr_ = corr(clusterTimeSeries1', clusterTimeSeries2');
                        pairwise_Correlation(j1, j2) = mean(cross_corr_(:));
                        pairwise_Correlation(j2, j1) = pairwise_Correlation(j1, j2);
                    end
                end
                title_str = sprintf("Allen Temp correlation: Sub %s, Sess %d", prev_sub, prev_ses);
                figure;
                imagesc(pairwise_Correlation);
                colorbar;
                title(title_str);
                saveas(gcf, fullfile(resultDir,['Allen_TempCorr_pairwise_',prev_sub,...
                    '_sess_',num2str(prev_ses), '.png']));
                close;
                
                % TemporalCorr
                clusterwise_within_corr = [];
                clusterwise_across_corr = [];
                % compute temporal correlation between pixels
                cluster_centres = zeros(num_clusters, 2);
                [nR, nC] = size(brain_mask);

                mergedA = zeros(nR*nC, num_clusters);
                for i1 = 1:num_clusters
                    t_id = find(allenMap_hem == clusterIDs(i1));
                    mergedA(t_id, i1) = 1;
                end

                % within parcel 
                clusterwise_within_corr = zeros(1, num_clusters);
                %parfor cl = 1:num_clusters  % @renu check the dFoF_masked for parfor 
                for cl = 1:num_clusters
                    cl_pix_ids = find(mergedA(:, cl)); % pixel ids as per original image
                    [~, cl_pix_pos_mask] = ismember(cl_pix_ids, allregionspix); % position of corresponding pixel time series in maksed data - PixxTime_dff
                    pix_tseries = dFoF_masked(cl_pix_pos_mask, :);
                    corr_within_mat = corr(pix_tseries');            

                    % Compute average correlation value
                    upT = triu(corr_within_mat, 1);
                    nzupT = upT(upT ~= 0);
                    clusterwise_within_corr(cl) = mean(nzupT);

                    [p_r, p_c] = ind2sub([nR, nC], cl_pix_ids);
                    cluster_centres(cl, 1) = mean(p_r);
                    cluster_centres(cl, 2) = mean(p_c);
                end
                
                % across parcels
                NUM_NEAREST = 2;
                cl_distances = pdist2(cluster_centres, cluster_centres);
                clusterwise_across_corr = zeros(1, num_clusters);
                for cl = 1:num_clusters
                    dists = cl_distances(cl, :);
                    [~, sorted_ids] = sort(dists);
                    nearest_nbrs = sorted_ids(2:NUM_NEAREST+1); % 1st one is always itself (zero dist)
                    
                    % Time series' of pixels in the current parcel
                    cl_pix_ids = find(mergedA(:, cl)); % pixel ids as per original image
                    [~, cl_pix_pos_mask] = ismember(cl_pix_ids, allregionspix); % position of corresponding pixel time series in maksed data - dFoF_masked
                    pix_tseries = dFoF_masked(cl_pix_pos_mask, :);
                    
                    cumul_corr = 0;
                    cumul_pix = 0;
                    for nn = 1:NUM_NEAREST
                        cl_pix_ids_nn = find(mergedA(:, nearest_nbrs(nn))); % pixel ids as per original image
                        [~, cl_pix_pos_mask_nn] = ismember(cl_pix_ids_nn, allregionspix); % position of corresponding pixel time series in maksed data - dFoF_masked
                        pix_tseries_nn = dFoF_masked(cl_pix_pos_mask_nn, :);

                        cross_corr_ = corr(pix_tseries', pix_tseries_nn');
                        cumul_corr = cumul_corr + sum(cross_corr_(:));
                        cumul_pix = cumul_pix + size(cross_corr_, 1)*size(cross_corr_, 2);
                    end
                    clusterwise_across_corr(cl) = cumul_corr / cumul_pix;
                end

                % save outputs
                outfname = ['sub_' prev_sub '_ses_' num2str(prev_ses) '_allen_out.mat'];
                storepath = fullfile(resultDir, 'run_allen_sessions');
                if ~exist(storepath, 'dir') 
                   mkdir(storepath); % Create the directory
                end
                fulloutpath = fullfile(storepath, outfname);
                results = struct();
                results.filename = outfname;
                results.labels = allenMap_hem;
                results.mergedA = mergedA;
                results.clusterwise_within_corr = clusterwise_within_corr;
                results.clusterwise_across_corr = clusterwise_across_corr;
                save(fulloutpath, '-struct', 'results');
    
                V = V0; % initialize to the new subject
            end
            prev_sub = sub_value;
            prev_ses = ses_value;
        end
    end
end

%%
if (RUN_ALLEN_REPORTING)
     % All processed files
    pFileList = dir(fullfile(AllenCorrStorePath, '*.mat'));

    % Group the files subject-wise
    pFileGrp = struct();
    for i = 1:length(pFileList)
        fname = pFileList(i).name;
        fpath = fullfile(AllenCorrStorePath, fname);

        sub_match = regexp(fname, 'sub_(\w{5})', 'tokens');
        ses_match = regexp(fname, 'ses_(\d+)', 'tokens');
        sub_value = sub_match{1}{1};
        ses_value = str2double(ses_match{1}{1});

        % Add the file to the subject group
        ses_key = sprintf('ses_%d', ses_value);
        if ~isfield(pFileGrp, sub_value)
            pFileGrp.(sub_value) = struct();
        end
        pFileGrp.(sub_value).(ses_key) = fpath;
    end

    subNames = fieldnames(pFileGrp);
    NUM_SESS = 3;
    mean_corr_within_subwise = zeros(1, length(subNames));
    mean_corr_across_subwise = zeros(1, length(subNames));
    var_corr_within_subwise = zeros(1, length(subNames));
    var_corr_across_subwise = zeros(1, length(subNames));

    min_corr_within_subwise = zeros(1, length(subNames));
    min_corr_across_subwise = zeros(1, length(subNames));
    max_corr_within_subwise = zeros(1, length(subNames));
    max_corr_across_subwise = zeros(1, length(subNames));
    displayCorr = zeros(length(subNames), 2);

    for s = 1:length(subNames)
        sub = subNames{s};
        sesNames = fieldnames(pFileGrp.(sub));
        corr_within_parcel = [];
        corr_across_parcel = [];
        for s1 = 1:length(sesNames)
            ses = sesNames{s1};
            dat_parcel = load(pFileGrp.(sub).(ses));
            corr_within_parcel = [corr_within_parcel, dat_parcel.clusterwise_within_corr];
            corr_across_parcel = [corr_across_parcel, dat_parcel.clusterwise_across_corr];
        end
        mean_corr_within_subwise(s) = mean(corr_within_parcel);
        mean_corr_across_subwise(s) = mean(corr_across_parcel);
        var_corr_within_subwise(s) = var(corr_within_parcel);
        var_corr_across_subwise(s) = var(corr_across_parcel);
        min_corr_within_subwise(s) = min(corr_within_parcel);
        min_corr_across_subwise(s) = min(corr_across_parcel);
        max_corr_within_subwise(s) = max(corr_within_parcel);
        max_corr_across_subwise(s) = max(corr_across_parcel);
        displayCorr(s, 1) = mean_corr_within_subwise(s);
        displayCorr(s, 2) = mean_corr_across_subwise(s);
    end
    subLabels = {'SLC01', 'SLC03', 'SLC04', 'SLC05', 'SLC06', 'SLC07', 'SLC08', 'SLC09', 'SLC10'};

    % Create figure
    figure;
    hold on;
    
    % Plotting box-like shapes using min, max, mean, and variance
    for i = 1:length(subLabels)
        % Calculate lower and upper bounds of variance
        lowerBound = mean_corr_within_subwise(i) - sqrt(var_corr_within_subwise(i));
        upperBound = mean_corr_within_subwise(i) + sqrt(var_corr_within_subwise(i));
        
        % Draw the box (representing variance range)
        patch([i-0.2, i+0.2, i+0.2, i-0.2], [lowerBound, lowerBound, upperBound, upperBound], 'b', 'FaceAlpha', 0.3);
    
        % Draw the line from min to max (whiskers)
        plot([i, i], [min_corr_within_subwise(i), max_corr_within_subwise(i)], 'k-', 'LineWidth', 1.5); 
        
        % Draw the mean as a red line inside the box
        plot([i-0.2, i+0.2], [mean_corr_within_subwise(i), mean_corr_within_subwise(i)], 'r-', 'LineWidth', 2);
    end
    
    % Set labels and title
    xlabel('Subjects');
    ylabel('Values');
    title('Temporal correlation within parcel');
    
    % Set x-ticks to match group indices and label them as subjects
    xticks(1:length(subLabels));
    xticklabels({'SLC01', 'SLC03', 'SLC04', 'SLC05', 'SLC06', 'SLC07', 'SLC08', 'SLC09', 'SLC10'});
    ylim([0 1]);
    hold off;

    figure;
    hold on;
    
    % Plotting box-like shapes using min, max, mean, and variance
    for i = 1:length(subLabels)
        % Calculate lower and upper bounds of variance
        lowerBound = mean_corr_across_subwise(i) - sqrt(var_corr_across_subwise(i));
        upperBound = mean_corr_across_subwise(i) + sqrt(var_corr_across_subwise(i));
        
        % Draw the box (representing variance range)
        patch([i-0.2, i+0.2, i+0.2, i-0.2], [lowerBound, lowerBound, upperBound, upperBound], 'b', 'FaceAlpha', 0.3);
    
        % Draw the line from min to max (whiskers)
        plot([i, i], [min_corr_across_subwise(i), max_corr_across_subwise(i)], 'k-', 'LineWidth', 1.5); 
        
        % Draw the mean as a red line inside the box
        plot([i-0.2, i+0.2], [mean_corr_across_subwise(i), mean_corr_across_subwise(i)], 'r-', 'LineWidth', 2);
    end

    storepath = fullfile(resultDir, 'stats_mat_files');
    if ~exist(storepath, 'dir') 
       mkdir(storepath); % Create the directory
    end
    save(fullfile(storepath, "Allen_corr_means.mat"), 'mean_corr_within_subwise', 'mean_corr_across_subwise');
    
    % Set labels and title
    xlabel('Subjects');
    ylabel('Values');
    title('Temporal correlation across parcels');
    
    % Set x-ticks to match group indices and label them as subjects
    xticks(1:length(subLabels));
    xticklabels({'SLC01', 'SLC03', 'SLC04', 'SLC05', 'SLC06', 'SLC07', 'SLC08', 'SLC09', 'SLC10'});
    ylim([0 1]);
    hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
