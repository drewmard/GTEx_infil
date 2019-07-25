# Age, Sex, and Genetics Influence the Abundance of Infiltrating Immune Cells in Human Tissues


## Robust estimation of immune cell types in bulk RNA-seq profiles
Directory: ./scripts/SynMix_Sims/


#### In silico simulation of synthetic mixes:
GenerateSyntheticMixture-ALL.sh \
(TPMs were subsequently generated by an in-house pipeline developed by Akanksha Verma, which has been described in the methods section)

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

#### Filtering of infiltration phenotypes (212 to 73):
filter.R

#### Infiltration signatures in bulk gene expression values
< insert script here >



## Identification and characterization of extreme infiltrating immune cell patterns
Directory: ./scripts/HotCold_Cluster/


#### Algorithm:
< munna's script here >

#### Widespread vs tissue-specific extreme inflammation patterns
< insert script here >



## Association of age and sex with immune infiltration
Directory: ./scripts/AgeSex/

#### Statistical analysis, Table 1a,b, and Fig 2a,b:
AgeSex_Analysis.R



## Association of genetic variants with infiltrating immune cells
Directory: ./scripts/GeneticAnalysis/


#### Generate PCs:
./GTEx_Genetic_PCA/tissue_subset_wrapper.sh \
(which runs ./GTEx_Genetic_PCA/PCA_calc.sh)

#### Prepare data for external GWAS software:
GWAS_preprocess.R

#### Run GWAS
GWAS_local.sh \
(or GWAS.sh in parallel cloud computing environment, followed up by running Merge_Chr_GWAS_wrapper.sh)

#### Combine related p-values using Empirical Brown's method, and identify significant results
Empirical_Brown_pval_wrapper.sh \
(runs Empirical_Brown_pval.R in parallel cloud computing environment)

#### Significant results analysis
significant_results.R

#### QQ plots
qqplots.R

#### Genotype - phenotype plots
extract_SNP.sh (extract SNP) \
genotypebyphenotype_plots.R (create plots)

#### Gene expression - phenotype plots
expressionbyphenotype_plots.R


## Downstream analysis of genetic results
Directory: ./scripts/GeneticAnalysis_2/


#### eQTL/ieQTL overlap (Method 1):
< insert script here >

#### eQTL/ieQTL overlap (Method 2):
< insert script here >

#### GeneMania:
< insert script here >

#### Pleiotropic effects:
< insert script here >



