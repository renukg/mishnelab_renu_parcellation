function ReorderedOutput = ReorderLikeAllen(ReorderInput, AllenFilePath, hemisp)
%function ReorderedOutput = ReorderLikeAllen(MatFilePath, AllenFilePath)
    % input_ = load(MatFilePath);
    % ReorderInput = input_.labels{1};

    allen_map = niftiread(AllenFilePath);
    % Resize the Allen map to the size of the brain image using nearest-neighbor interpolation
    allenMapResized = imresize(allen_map, [60, 81], 'nearest');
    MIDLINE1 = 29;
    MIDLINE2 = 32;
    allenMap_top = allenMapResized(1:MIDLINE1-1, :);
    allenMap_bottom = allenMapResized(MIDLINE2+1:end, :, :);
    allenMap_bottom_flipped = flip(allenMap_bottom, 1);

    ROW_PAD = 4;
    allenMap_hem = zeros(28+ROW_PAD, 81);
    %hemisp = 1;
    if (hemisp == 1)
        allenMap_hem(1:MIDLINE1-1, :) = allenMap_top;
    else
        allenMap_hem(1:MIDLINE1-1, :) = allenMap_bottom_flipped;
    end

    allenMap_hem = rot90(allenMap_hem, -1);
    allen_clusters = unique(allenMap_hem);
    allen_clusters(allen_clusters == 0) = [];
    allen_clusters = sort(allen_clusters);

    inp_clusters = unique(ReorderInput);
    inp_clusters(inp_clusters == 0) = [];
    
    Group_clust_AllenBased = struct(); 
    for i0 = 1:length(allen_clusters)
       Group_clust_AllenBased.(sprintf('A%d', allen_clusters(i0))) = [];
    end
    
    MapClust_To_AllenIDs = zeros(size(inp_clusters));
    for i=1:length(inp_clusters)
       cur_clust_mask = (ReorderInput == inp_clusters(i));
       cur_clusters_Allen = allenMap_hem(cur_clust_mask);
       cur_clusters_Allen(cur_clusters_Allen == 0) = []; % remove cluster '0'
       assign_AllenID = mode(cur_clusters_Allen);
       sfield = sprintf('A%d', assign_AllenID);
       Group_clust_AllenBased.(sfield) = [Group_clust_AllenBased.(sfield), inp_clusters(i)];
       MapClust_To_AllenIDs(i) = assign_AllenID;
    end
    
    ReorderedOutput = zeros(size(ReorderInput));
    Afields = fieldnames(Group_clust_AllenBased);
    NewID = 1;
    for ai = 1:length(Afields)
        cur_clusters = Group_clust_AllenBased.(Afields{ai});
        for c1 = 1:length(cur_clusters)
            pmask = (ReorderInput == cur_clusters(c1));
            ReorderedOutput(pmask) = NewID;
            NewID = NewID + 1;
        end
    end   
    % figure;
    % imagesc(label2rgb(allenMap_hem));
    % figure;
    % imagesc(label2rgb(ReorderInput));
    % figure;
    % imagesc(label2rgb(ReorderedOutput));
end