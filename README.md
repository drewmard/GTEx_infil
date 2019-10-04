# Age, Sex, and Genetics Influence the Abundance of Infiltrating Immune Cells in Human Tissues


## Robust estimation of immune cell types in bulk RNA-seq profiles
Directory: ./scripts/SynMix_Sims/


#### In silico simulation of synthetic mixes:
GenerateSyntheticMixture-ALL.sh \
(TPMs were subsequently generated by an in-house pipeline developed by Akanksha Verma, which has been described in the methods section)

#### Deconvolution of TPMs from in silico mixes
cibersort_simulations.R # using cibersort \ 
xcell_sim_gen.R # using xcell \

#### Comparison of CIB-Rel, CIB-Abs, and xCell
xcell_ciber_simComp.R



## Evaluating infiltration across human tissues by using deconvolution
Directory: ./scripts/GTEx_Deconv/


#### CIBERSORT: 
running_cibersort.R \
cibersort_out.R

#### xCell: 
xCell_generation.R \
xCell_process.R

#### Pairwise correlations of infiltration phenotypes (Supp Fig. 3):
corr_heatmaps.R

#### Case examples of differences between deconvolution algorithms (Supp Fig. 1):
xcell_cibersort_CaseExample.R

#### Hierarchical clustering (Fig 2a, Supp Fig 2):
clustering_heatmap.R

#### Boxplots of immune content (Fig 2b, Supp Fig 4-7):
infiltration_boxplots.R \
celltype_boxplots.R

#### t-SNE of immune content (Supp Fig 8):
immune_content_clusters_tsne.R

#### Filtering of infiltration phenotypes (190):
filter.R

#### Infiltration signatures in bulk gene expression values
generate_exprPCA_matrix.R \
calculate_exprPCA_infil_cor.R



## Identification and characterization of extreme infiltrating immune cell patterns
Directory: ./scripts/HotCold_Cluster/


#### Algorithm:
< munna's script here >

#### Widespread vs tissue-specific extreme inflammation patterns
hot_tissueSpecific.R



## Association of age and sex with immune infiltration
Directory: ./scripts/AgeSex/

#### Statistical analysis:
AgeSex_Analysis.R

#### Create age/sex plots of most significant associations (Fig 3b-e)
AgeSex_Plots2.R

#### Heatmap summary of age/sex association results (Fig 3f-g)
age_sex_heatmap.R

#### Breast t-sne, males vs females (Supplementary Fig 11)
breast_content_clusters_tsne.R


## Association of genetic variants with infiltrating immune cells
Directory: ./scripts/GeneticAnalysis/


#### Generate PCs:
./GTEx_Genetic_PCA/PCA_calc_all.sh

#### Prepare data for external GWAS software:
GWAS_preprocess.R

#### Run GWAS
GWAS_local.sh (or GWAS.sh in parallel cloud computing environment) \
Merge_Chr_GWAS_wrapper.sh (to merge different chromosomes together)

#### Combine related p-values using Empirical Brown's method, and identify significant results
Empirical_Brown_pval_wrapper.sh \
(runs Empirical_Brown_pval.R in parallel computing environment)

#### Analyzing significant results in genetic analysis.
significant_results.R

#### Script to extract a SNP for genotype - phenotype plots
extract_SNP.sh

#### Figure panel: qqplots, genotype-phenotype plots & expression-phenotype plots
figure_panel.R \
separate_qq.R

## Downstream analysis of genetic results
Directory: ./scripts/GeneticAnalysis_2/

#### Analysis of related genetic variants also associated with thyroiditis using GeneHancer regions & phenoscanner
GWAS_search.R

#### ieQTL enrichment:
eQTL_enrichment_method2.R \
eqtl_enrich_plots.R

#### GeneMania:
eQTL_network_gen.R (use eQTL_network_gen.R for *ieGene.txt output for input into GeneMania server)

#### Pleiotropic effects:
pleiotropy.R

#### Venn diagrams (Fig 3a, Fig 4a-b)
venndiagram.R


