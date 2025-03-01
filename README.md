# Validation of LSSC parcellation algorithm for fMRI data

Data and results: https://drive.google.com/drive/folders/1VuCujrljncQ8kfwf0e_B2W9Na8W1Nj48?usp=drive_link
## Please follow the instructions below to run parcellation on the FMRI data:
<table border="1">
  <tr>
    <th>Parcellation Technique</th>
    <th>Code file</th>
    <th>Instructions</th>
  </tr>
  <tr>
    <td>LSSC</td>
    <td>run_LSSC_fmri_sessions_hemisphere.m</td>
    <td> 
      <ul>
        <li>Update the directories. </li>
        <li>dataID = 1/2/3/4 - for 4 data sets. </li>
        <li>RUN_LSSC = 0/1 - to run and generate LSSC parcellation images, mat files.</li>
        <li>RUN_DICE_SIMILARITY = 0/1 - Dice computation subject-wise. Run only after LSSC.</li>
        <li>RUN_TEMPORAL_CORR = 0/1 - Within and across temporal correlations subject-wise. Run only after LSSC.</li>
        <li>Look for "cfg" structure and "runROI_meso_nlm_new_v2.m" file to modify LSSC parameters.</li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>Kmeans</td>
    <td>run_Kmeans_fmri_sessions_hemisphere_v3.m</td>
    <td>
      <ul>
        <li>Update the directories. </li>
        <li>dataID = 1/2/3/4 - for 4 data sets. </li>
        <li>RUN_KNN = 0/1 - to run and generate Kmeans parcellation images, mat files.</li>
        <li>RUN_DICE_SIMILARITY = 0/1 - Dice computation subject-wise. Run only after Kmeans.</li>
        <li>RUN_TEMPORAL_CORR = 0/1 - Within and across temporal correlations subject-wise. Run only after Kmeans.</li>
        <li>Paramters to modify: N_KNN_CLUSTERS=27, min_clust_size=15</li>
      </ul>
    </td>
  </tr>
  <tr>
    <td>Allen</td>
    <td>allenmaps_v1.m</td>
    <td>
      <ul>
        <li>Update the data and Allen atlas (2D_calcium_atlas.nii) directories. </li>
        <li>dataID = 1/2/3/4 - for 4 data sets. </li>
        <li>RUN_ALLEN_PROCESSING = 0/1 - to process data, generate time series and pairwise correlation plots, and generate mat files for correlation</li>
        <li>RUN_ALLEN_REPORTING = 0/1 - Within and across temporal correlations subject-wise. Run only after ALLEN_PROCESSING.</li>
      </ul>
    </td>
  </tr>
</table>

## To perform dice, within/across correlation comparison study across LSCC, Kmeans, and Allen (anatomical):
Code file: compare_LSSC_kNN_2.m <br>
Instructions:
<ul>
  <li>Update the directories. </li>
  <li>Images will be saved in the "saveimagesPath". </li>
  <li>dataID = 1/2/3/4 - for 4 data sets. </li>
  <li>Three images will be generated for each data set - Dice comparison, Withn correlation comparison, Across correlation comparison</li>
</ul>
