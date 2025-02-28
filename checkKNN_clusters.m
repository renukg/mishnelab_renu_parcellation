clc;
clear;
close all;

% Directory containing the files
dataDir = 'D:\UCSD_Acads\ProfGal_Research\data32\fMRIData\REST';

% Results diectory
resultDir = 'D:\UCSD_Acads\ProfGal_Research\test_run_KNN25_hemisphere_replicas100_v3_minclstPix_15_new';
processed_dir = fullfile(resultDir, 'run_knn_sessions');
    
% All processed files
pFileList = dir(fullfile(processed_dir, '*.mat'));

DISPLAY_CLUSTERS = 0;

if DISPLAY_CLUSTERS == 1
    figure;
end
for fi = 1:length(pFileList)
    filename = pFileList(fi).name;
    %filename = 'sub_SLC06_ses_2_hemisp_1_knn_out.mat';
    cur_data = load(fullfile(processed_dir, filename));
    cur_labels = cur_data.labels{1};

    clusters = unique(cur_labels); 
    clusters(clusters == 0) = []; % Remove zero id

    brain_mask = (cur_labels > 0);
    
    if DISPLAY_CLUSTERS == 1
        % Loop through each cluster and display it
        for i = 1:length(clusters)
            cluster_id = clusters(i);
            cluster_mask = (cur_labels == cluster_id);
            vis_image = zeros(size(cur_labels)); 
            vis_image(brain_mask) = 1; 
            vis_image(cluster_mask) = 5; 
            
            % Display the cluster
            imagesc(label2rgb(vis_image)); 
            colormap('gray'); 
            title(['Cluster ' num2str(cluster_id)]); 
            axis off; 
            
            % Pause 
            pause(1); 
        end
    end

    for c1 = 1:length(clusters)
        clus = clusters(c1);
        cur_cluster_mask = (cur_labels == clus);
        CC = bwconncomp(cur_cluster_mask);
        num_connected_comps = CC.NumObjects;
        if (num_connected_comps == 1)
            continue;
        end
        clus
    end
    nj=10;

end