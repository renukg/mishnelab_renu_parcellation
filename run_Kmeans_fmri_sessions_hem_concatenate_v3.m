clc;
clear;
close all;

dataID = 1;

if dataID == 1
    dataDir = 'D:\UCSD_Acads\ProfGal_Research\data32\fMRIData\REST';
    resultDir = 'D:\UCSD_Acads\ProfGal_Research\test_run_KNN27_hem_concatenate_replic100_minclstPix_15_AllenOrdered';
end
if dataID == 2
    dataDir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\GSR_1voxelmask_new';
    resultDir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\GSR_1voxelmask_results\test_run_KNN27_hem_concatenate_replic100_minclstPix_15_AllenOrdered';
end
if dataID == 3
    dataDir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\GSR_thickermask_new';
    resultDir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\GSR_thickermask_results\test_run_KNN27_hem_concatenate_replic100_minclstPix_15_AllenOrdered';
end
if dataID == 4
    dataDir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\noGSR_thickermask_new';
    resultDir = 'D:\UCSD_Acads\ProfGal_Research\data_fMRI_GSR_processed\noGSR_thickermask_results\test_run_KNN27_hem_concatenate_replic100_minclstPix_15_AllenOrdered';
end

% All sessions in the directory
sesList = dir(fullfile(dataDir, 'session-*'));

% Allen map
mapDir = 'D:\UCSD_Acads\ProfGal_Research\Allen maps';
AllenFilePath = fullfile(mapDir, '2D_calcium_atlas.nii'); 

RUN_KNN = 0;
RUN_DICE_SIMILARITY = 1;
RUN_TEMPORAL_CORR = 1;
cfg.thrcluster=[0.9]; % @renu - handle this well 

if (RUN_KNN)
    % Run KNN after concatenating all runs of a session for each subject
    for ses = 1:length(sesList) % processing each session
        sesname = sesList(ses).name;
        sesDir = [dataDir, '\', sesname];
        fileList = dir(fullfile(sesDir, '*.nii.gz'));
        
        V = [];
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
            else      % New subject - perform LSSC to V matrix and initialize to new subject's data
                if (i == length(fileList)) % if last file
                    prev_sub = sub_value;
                    prev_ses = ses_value;
                end
                % Run for the current data (V)
                num_time_samples = size(V, 3);
                
                MIDLINE1 = 29;
                MIDLINE2 = 32;
                % Separating the brain region into two halves along the
                % midline
                % Extract the top and bottom halves
                V_top = V(1:MIDLINE1-1, :, :);        % Top half
                V_bottom = V(MIDLINE2+1:end, :, :);  % Bottom half
                
                % Flip the bottom half vertically (along the row axis)
                V_bottom_flipped = flip(V_bottom, 1);  % Flip along the 1st dimension (rows)
                
                num_time_frames = size(V_top, 3);
                
                % Pre-allocate a new matrix to hold the interleaved data
                ROW_PAD = 4; % To avoid out of index error in LSSC algorithm
                %V_combined_h = zeros(size(V_top, 1)+ROW_PAD, size(V_top, 2), num_time_frames);
                V_combined_h = zeros(size(V_top, 1)+ROW_PAD, size(V_top, 2), 2*num_time_frames);

                %Interleave the top and bottom halves in time
                for t = 1:num_time_frames
                    V_combined_h(1:MIDLINE1-1, :, 2*t-1) = V_top(:, :, t);  % Odd indices: V_top
                    V_combined_h(1:MIDLINE1-1, :, 2*t) = V_bottom_flipped(:, :, t);  % Even indices: V_bottom_flipped
                end

                % Create brain mask for the hemispheres
                maxV = max(V_combined_h, [], 3);
                brain_mask = maxV ~= 0;
                [R,C] = size(brain_mask);

                % Flatten the spatial dimensions
                V_flat = reshape(V_combined_h, [], size(V_combined_h, 3));
                
                % Extract relevant pixels
                allregionspix = find(brain_mask);
                dFoF_masked = V_flat(allregionspix, :);
                nPIX = size(dFoF_masked, 1); % number of pixels in brain region
                
                % KNN Prameters
                N_KNN_CLUSTERS = 27;
                cfg.ComputeTemporalCorr = true;
                
                % Apply k-means clustering with Euclidean distance
                % metric
                %[idx, clust_centers] = kmeans(dFoF_masked, N_KNN_CLUSTERS, 'MaxIter', 1000, 'Replicates', 5);
                [PixIdx, clust_centers] = kmeans(dFoF_masked, N_KNN_CLUSTERS, 'MaxIter', 1000, 'Replicates', 100);
                % PixIdx stores cluster assignment for each pixel

                % Identify disconnected clusters and assign different cluster ids
                cur_labels = zeros(R, C);
                cur_labels(allregionspix) = PixIdx;
                cur_clusters = unique(PixIdx);
                cur_labels_new = zeros(size(cur_labels));
                AssignID = 1;
                for c1 = 1:length(cur_clusters)
                    clus = cur_clusters(c1);
                    cur_cluster_mask = (cur_labels == clus);
                    CC = bwconncomp(cur_cluster_mask);
                    num_connected_comps = CC.NumObjects;
                    if (num_connected_comps == 1)
                        cur_pix_grp = CC.PixelIdxList{1};
                        cur_labels_new(cur_pix_grp) = AssignID;
                        AssignID = AssignID + 1;
                        continue;
                    end
                    PixGroups = CC.PixelIdxList;
                    for c2 = 1:num_connected_comps
                        cur_pix_grp = PixGroups{c2};
                        cur_labels_new(cur_pix_grp) = AssignID;
                        AssignID = AssignID + 1;
                    end
                end
                

                % Reassign clusters with less pixels to nearby clusters
                % find size of each cluster
                PixIdx = cur_labels_new(allregionspix);
                clusters = unique(cur_labels_new);
                clusters(clusters == 0) = [];
                num_clust = length(clusters);
                clust_sizes = zeros(num_clust, 1); % no. of pixels per cluster
                for ic = 1:num_clust
                    clust_sizes(ic) = sum(PixIdx == (ic));
                end

                % Display each cluster before reassigning
                % figure;
                % clusters = unique(cur_labels_new); 
                % clusters(clusters == 0) = []; % Remove zero id
                % 
                % brain_mask = (cur_labels_new > 0);
                % 
                % % Loop through each cluster and display it
                % for i = 1:length(clusters)
                %     cluster_id = clusters(i);
                %     cluster_mask = (cur_labels_new == cluster_id);
                %     vis_image = zeros(size(cur_labels_new)); 
                %     vis_image(brain_mask) = 1; 
                %     vis_image(cluster_mask) = 5; 
                % 
                %     % Display the cluster
                %     imagesc(label2rgb(vis_image)); 
                %     colormap('gray'); 
                %     title(['Cluster ' num2str(cluster_id)]); 
                %     axis off; 
                % 
                %     % Pause 
                %     pause(1); 
                % end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                min_clust_size = 15; % min pixels in a cluster
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Reject pixels
                reject_clust = find(clust_sizes < min_clust_size);
                reject_ids = ismember(PixIdx, reject_clust);
                PixIdx(reject_ids) = 0; % Mark rejected pixels as unassigned (0)
                reject_pixels = allregionspix(reject_ids);

                % Reassign rejected pixels
                Pix2Resolve = reject_pixels;
                prev_num_pix2resolve = 0;
                cnt_iters = 0;
                while (~isempty(Pix2Resolve))
                    for ir = 1:length(Pix2Resolve)
                        cur_pix = Pix2Resolve(ir);
                        cur_pix_id = find(allregionspix == cur_pix);
                        pixel_data = dFoF_masked(cur_pix_id, :); % Extract data for the rejected pixel
                        
                        % Convert current pixel index to subscripts (row, column)
                        [f, g] = ind2sub([R, C], cur_pix);
                        
                        % Create a grid of neighboring pixel subscripts
                        [ff, gg] = meshgrid(f + [-1:1], g + [-1:1]);

                        % Ensure indices are within bounds
                        if any(ff(:) < 1 | ff(:) > R) || any(gg(:) < 1 | gg(:) > C)
                            error('Indices in ff or gg are out of range.');
                        end
                        
                        % Convert neighboring subscripts to linear indices
                        nbr_pix = sub2ind([R, C], ff(:), gg(:));
                        
                        % Remove indices that correspond to non-region pixels
                        nbr_pix = setdiff(nbr_pix, find(brain_mask == false));
                        %nbr_pix_ids = find(allregionspix == nbr_pix);
                        [isNbr, nbr_pix_ids] = ismember(nbr_pix, allregionspix);

                        % Check if all neighboring pixels are missing or have more than 4 missing pixels
                        LS_PIX = 4;
                        if (prev_num_pix2resolve == length(Pix2Resolve))
                            LS_PIX = 8;
                        end
                        if sum(PixIdx(nbr_pix_ids) == 0) > LS_PIX  
                            continue;
                        end

                        % Find max correlation among neighbours that
                        % belong to a cluster
                        maxcorr = -Inf;
                        max_cor_clust = 0;
                        for jr = 1:length(nbr_pix)
                            if PixIdx(nbr_pix_ids(jr)) ~= 0
                                cluster_id = PixIdx(nbr_pix_ids(jr));
                                cluster_data = dFoF_masked(PixIdx == cluster_id, :); % Time series of all pixels in the cluster
                                mean_cluster_data = mean(cluster_data, 1); % Mean time series of the cluster
                                cur_corr = corr(pixel_data', mean_cluster_data'); % Correlation
                                if maxcorr < cur_corr
                                    maxcorr = cur_corr;
                                    max_cor_clust = cluster_id;
                                end
                            end
                        end
                        % Assign pixel to the ROI with the highest correlation
                        PixIdx(cur_pix_id) = max_cor_clust;
                    end
                    prev_num_pix2resolve = length(Pix2Resolve);
                    reject_ids = find(PixIdx == 0);
                    Pix2Resolve = allregionspix(reject_ids);
                    %prev_num_pix2resolve
                    cnt_iters = cnt_iters + 1;
                    % if cnt_iters > 5
                    %     stap=9;
                    % end
                end

                % Renaming clusters ids to be in a sequence as 1,2,3...
                idx_cpy = PixIdx;
                PixIdx = zeros(size(idx_cpy));
                clust_ids = unique(idx_cpy);
                for ic1 = 1:length(clust_ids)
                    tp_pix = find(idx_cpy == clust_ids(ic1));
                    PixIdx(tp_pix) = ic1;
                end

                % Proceed
                labels = zeros(R, C);
                labels(allregionspix) = PixIdx;

                % Reorder cluster IDs as per Allen cluster IDs
                labels_before = labels;
                hemisp = 1;
                labels = ReorderLikeAllen(rot90(labels, -1), AllenFilePath, hemisp);
                labels = rot90(labels, 1);

                clust_idx = unique(labels);
                clust_idx(clust_idx==0) = [];
                mergedA = zeros(R*C, length(clust_idx));
                for i1 = 1:length(clust_idx)
                    t_id = find(labels == clust_idx(i1));
                    mergedA(t_id, clust_idx(i1)) = 1;
                end

                title_str = ['sub_' prev_sub '_session_' num2str(prev_ses)];
                imagesc(label2rgb(labels));
                saveas(gcf, fullfile(resultDir,['KNN_Labels_',title_str,'_nkNN_', num2str(N_KNN_CLUSTERS),'.png']));
                
                % TemporalCorr
                clusterwise_within_corr = [];
                clusterwise_across_corr = [];
                if cfg.ComputeTemporalCorr
                    % compute temporal correlation between pixels
                    num_clusters = size(mergedA, 2);
                    cluster_centres = zeros(num_clusters, 2);
                    [nR, nC] = size(brain_mask);
    
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
                end

                % save outputs
                labels_all = cell(length(cfg.thrcluster), 1); % currently only 1
                labels_all{1} = labels;
                mergedA_all = cell(length(cfg.thrcluster), 1);
                mergedA_all{1} = mergedA;
                outfname = ['sub_' prev_sub '_ses_' num2str(prev_ses) '_knn_out.mat'];
                storepath = fullfile(resultDir, 'run_knn_sessions');
                if ~exist(storepath, 'dir') 
                   mkdir(storepath); % Create the directory
                end
                fulloutpath = fullfile(storepath, outfname);
                results = struct();
                results.filename = outfname;
                results.labels = labels_all;
                results.mergedA = mergedA_all;
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
if (RUN_DICE_SIMILARITY)
    processed_dir = fullfile(resultDir, 'run_knn_sessions');
    
    % All processed files
    pFileList = dir(fullfile(processed_dir, '*.mat'));

    % Group the files subject-wise
    pFileGrp = struct();
    for i = 1:length(pFileList)
        fname = pFileList(i).name;
        fpath = fullfile(processed_dir, fname);

        sub_match = regexp(fname, 'sub_(\w{5})', 'tokens');
        ses_match = regexp(fname, 'ses_(\d+)', 'tokens');
        sub_value = sub_match{1}{1};
        ses_value = str2double(ses_match{1}{1});

        % Add the file to the subject group
        ses_key = sprintf('ses_%d', ses_value);
        if ~isfield(pFileGrp, sub_value)
            pFileGrp.(sub_value) = struct();
        end
        pFileGrp.(sub_value).(ses_key)= fpath;
    end

    % Compute Dice similarity between all the session-pairs for each subject
    subNames = fieldnames(pFileGrp);
    NUM_SES_PAIRS = 3;
    displayDice = zeros(length(subNames), NUM_SES_PAIRS);
    for i = 1:length(subNames)
        sub = subNames{i};

        dice_pairwise = cell(1, NUM_SES_PAIRS);
        pairs = cell(1, NUM_SES_PAIRS);
        for j = 1:NUM_SES_PAIRS % 3 session pairs
            r1 = j;
            r2 = rem(j, NUM_SES_PAIRS) + 1;
            ke1 = sprintf('ses_%d', r1);
            ke2 = sprintf('ses_%d', r2);

            LSSC_out_pair = cell(1, 2);
            LSSC_out_pair{1} = load(pFileGrp.(sub).(ke1));
            LSSC_out_pair{2} = load(pFileGrp.(sub).(ke2));
            
            % Compute dice coefficient for each cluster threshold value
            
            dice_values = zeros(1, length(cfg.thrcluster));
            for thr_id = 1:length(cfg.thrcluster)
                % Pairing similar cluster labels 
                pairing = pairComponents(LSSC_out_pair{1}.mergedA{thr_id}, LSSC_out_pair{2}.mergedA{thr_id});
                p_1 = pairing.p1;
                p_2 = pairing.p2;
            
                labels_1 = LSSC_out_pair{1}.labels{thr_id};
                labels_2 = LSSC_out_pair{2}.labels{thr_id};
            
                % Re-labeling labels for 2nd one
                for k1=1:size(labels_2, 1)
                    for k2=1:size(labels_2, 2)
                        if (labels_2(k1, k2) ~= 0)
                            tp = find(p_2 == labels_2(k1, k2));
                            if (isempty(tp))
                                labels_2(k1, k2) = 0;
                            else
                                labels_2(k1, k2) = p_1(tp);
                            end
                            
                        end
                    end
                end
            
                dice_similarity = multiclass_dice_coefficient(labels_1, labels_2);
                dice_values(thr_id) = dice_similarity;
            end
            dice_pairwise{j} = dice_values;
            pairs{j} = sprintf('ses_%d-%d', r1, r2);
            displayDice(i, j) = dice_values;
        end
        pFileGrp.(sub).dice = dice_pairwise;
        pFileGrp.(sub).pairs = pairs;
    end
    
    data_dice = displayDice;
    %data_dice(:, 4) = mean(displayDice, 2);
    subLabels = {'SLC01', 'SLC03', 'SLC04', 'SLC05', 'SLC06', 'SLC07', 'SLC08', 'SLC09', 'SLC10'};
    groupLabels = {'Session 1-2', 'Session 2-3', 'Session 3-1'};
    figure;
    %bar(data_dice);
    bar(displayDice);
    set(gca, 'XTickLabel', subLabels);
    title('Comparison between sessions');
    xlabel('Subjects');
    ylabel('Dice value');
    legend('Session 1-2', 'Session 2-3', 'Session 3-1');
    ylim([0 1]);

    figure;
    bar(mean(displayDice, 2));
    set(gca, 'XTickLabel', subLabels);
    title('Comparison between sessions');
    xlabel('Subjects');
    ylabel('Dice value');
    legend('Average');
    ylim([0 1]);

    avg_dice_subwise = mean(displayDice, 2);
    storepath = fullfile(resultDir, 'stats_mat_files');
    if ~exist(storepath, 'dir') 
       mkdir(storepath); % Create the directory
    end
    save(fullfile(storepath, "knn_avg_dice.mat"), 'avg_dice_subwise');
end
%%
if (RUN_TEMPORAL_CORR)
     % All processed files
    processed_dir = fullfile(resultDir, 'run_knn_sessions');
    pFileList = dir(fullfile(processed_dir, '*.mat'));

    % Group the files subject-wise
    pFileGrp = struct();
    for i = 1:length(pFileList)
        fname = pFileList(i).name;
        fpath = fullfile(processed_dir, fname);

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
    %groupLabels = {'Within\_parcel', 'Across\_parcel'};
    % figure;
    % boxplot(displayCorr);
    % %bar(displayCorr);
    % set(gca, 'XTickLabel', subLabels);
    % title('Average Temporal correlation between pixels');
    % xlabel('Subjects');
    % ylabel('Temporal correlation');
    % legend('Within\_parcel', 'Across\_parcel');
    % ylim([0 1]);

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
    save(fullfile(storepath, "KNN_corr_means.mat"), 'mean_corr_within_subwise', 'mean_corr_across_subwise');
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%