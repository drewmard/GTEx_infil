# Age, Sex, and Genetics Influence the Abundance of Infiltrating Immune Cells in Human Tissues


## Robust estimation of immune cell types in bulk RNA-seq profiles
Directory: ./scripts/SynMix_Sims/


#### In silico simulation of synthetic mixes:
GenerateSyntheticMixture-ALL.sh

#### Deconvolution of TPMs from in silico mixes
cibersort_simulations.R \
xcell_sim_gen.R \
xcell_ciber_simComp.R



## Evaluating infiltration across human tissues by using deconvolution
Directory: ./scripts/GTEx_Deconv/


#### CIBERSORT: 
running_cibersort.R \
cibersort_out.R

#### xCell: 
xCell_generation.R
xCell_process.R

#### Heatmaps of immune content (Supp Fig. 1):
heatmaps.R

#### Pairwise correlations of infiltration phenotypes (Supp Fig. 2):
corr_heatmaps.R

#### Case examples of differences between deconvolution algorithms (Supp Fig. 3):
xcell_cibersort_CaseExample.R

#### Hierarchical clustering (Fig 1b, Supp Fig 4-5):
< insert script here >

#### Boxplots of immune content (Fig 1c, Supp Fig 6-9):
infiltration_boxplots.R

#### t-SNE of immune content (Supp Fig 10):
< insert script here >



## Evaluating infiltration across human tissues by using deconvolution
Directory: ./scripts/<insert directory>/






