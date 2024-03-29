## Introduction
Codes and datasets for the <br>"A Dorsomedial Prefrontal Cortex-based Dynamic Functional Connectivity Model of Rumination", <br> Nat. Comm., 2023
https://www.nature.com/articles/s41467-023-39142-9

## Dependencies
All the analyses are running based on <br>
https://github.com/cocoanlab/CocoanCore <br>
https://github.com/canlab/CanlabCore <br> 
https://github.com/spm/spm12

## Descriptions

### data
"Datasets.mat": <br>
    Each struct "Study1", "Study2", and ""Study3" holds DCC variance data from 
    all DMN seeds, and also RRS subscales. <br>
    Struct "Supp" holds DCC variance data from "pre" and "post" movie run based on the dMPFC seed. <br>
    It also includes RRS subscales, "alertness" and "self-relevance" ratings which were reported in <br>
    Supplementary Figure 6. <br><br>

"summaries.mat": <br>
    Summary statistiscs and  information comprising model <br>
    and parcellation information. <br><br>

"cluster_Fan_Net_r280.mat" & "Fan_et_al_atlas_r280.nii": <br>
    Parcellation Information. <br>
    (also could found in https://github.com/cocoanlab/cocoanCORE/Canonical_brains/Brainnetome) <br><br>

"dMPFCmodel.mat" & "dmpfc_based_connectivity_model.nii": <br>
    Final (full) model in mat & nifti format. <br><br>

"negative_weights.nii" & "positive_weights.nii": <br>
    Final (full) model divided in positive and negative weights for visualization. <br><br>

"dmpfc_based_connectivity_model_refined.nii": <br>
    Refined model with 21 weights. <br><br>
    
### codes
"makefigs.m": <br>
    All the codes in step for generating main figures.  <br>
    Includes virtual lesion analysis. <br><br>
"functions": <br>
    Folder containing collateral codes for making figures.
    



