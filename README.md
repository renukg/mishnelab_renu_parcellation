# Validation of LSSC parcellation algorithm for fMRI data

<table border="1">
  <tr>
    <th>Parcellation Technique</th>
    <th>Code file</th>
    <th>Instructions</th>
  </tr>
  <tr>
    <td>LSSC</td>
    <td>run_LSSC_fmri_sessions_hemisphere.m</td>
    <td> Update the directories <br>
      dataID = 1/2/3/4 - for 4 data sets <br>
      RUN_LSSC = 0/1 - to run and generate parcellation images, mat files.<br>
      RUN_DICE_SIMILARITY = 0/1 - Dice computation subject-wise. Run only after LSSC.<br>
      RUN_TEMPORAL_CORR = 0/1 - Within and across temporal correlations subject-wise. Run only after LSSC.
    </td>
  </tr>
  <tr>
    <td>Kmeans</td>
    <td>Value 5</td>
    <td>Value 6</td>
  </tr>
  <tr>
    <td>Allen</td>
    <td>Value 5</td>
    <td>Value 6</td>
  </tr>
</table>
