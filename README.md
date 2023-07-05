# NLC_differentiation

This GitHub repository contains all the scripts and main files that were used for the paper .... analyses and figures. Additonal data and processed files can be available upon reasonable request to the corresponding authors.

# CellChat
Scripts made by B. Guidoux D'Halluin, A. Camus

# Complex Heatmap
- Process_AxialSkeletonWeek7_and_draw_heatmap.R: script made to prepare Zhou et al. dataset for heatmap generation
- Process_NotochordWeek8_and_draw_heatmap.R: script made to prepare fetal notochord dataset from this study for heatmap generation
- Process_iPS_ESC-NLC_and_draw_heatmap.R: script made to prepare iPS/ESC-NLC integrated dataset from this study for heatmap generation
  
Scripts made by B. Guidoux D'Halluin, A. Camus

# Correlation
- prepare_data_and_perform_correlation.R: script used to perform correlation analyses between the notochord clusters from iPS-NLC and ESC-NLC dataset.

Script made by V. Tam

# Enrichment Map
- TOP100_DEG_AxialSkeletonWeek7.csv: gene list of the top 100 differentially expressed genes (DEG) from AxialSkeletonWeek7 dataset used for comparison
- TOP100_DEG_NotochordWeek8.csv: gene list of the top 100 DEG from NotochordWeek8 dataset used for comparison
- TOP100_DEG_iPS_ESC-NLC.csv: gene list of the top 100 DEG from iPS/ESC-NLC integrated dataset used for comparison
- EnrichmentMap_integrated_invitro_vs_invivo_all_integrated_scripts.Rmd: script used to apply hypergeometric comparison test for Enrichment Map creation

Scripts made by D. Yin

# Gene lists
- TOP50_DEG_ESC-NLC.xlsx: gene list of the top 50 DEG used for dataset annotation
- TOP50_DEG_NotochordWeek8.csv: gene list of the top 50 DEG used for dataset annotation
- TOP50_DEG_iPS-NLC.xlsx: gene list of the top 50 DEG used for dataset annotation
- TOP50_DEG_iPSESC-NLC.xlsx: gene list of the top 50 DEG used for dataset annotation
  
Generated by J. Warin, A. Camus
# Integration
- AxialSkeletonWeek7.R: integration of the sample 3 and 4 of week 7 embryo from Zhou et al. 2023 (https://doi.org/10.1002/advs.202206296)
- AxialSkeletonWeek7_NotochordWeek8_Harmony_integration.R: selection of the samples to be integrated between fetal notochord datasets (Zhou et al. 2023 and this study)
- iPS-NLC_ESC-NLC_Harmony_integration.R: integration of the iPS-NLC and ESC-NLC datasets using Harmony
  
Scripts made by D. Yin
# Matrisome
- Matrisome_heatmap.R: script to generate matrisome analysis in NotochordWeek8 dataset and compare it to iPS/ESC-NLC dataset
- Module_score.R: script to calculate the module scorev for matrisome gene set and notochord gene set
  
Script made by V. Tam
# Preprocess
Scripts made by V. Tam, F. Riemers, J. Warin, A. Camus
# SCENIC
Scripts made by J. Warin, L. David and A. Camus
# Velocity
Scripts made by X. Lin
# VennDiagram
Scripts made by J. Warin and A. Camus
# WGCNA
Based on Meistermann et al. 2021 (https://doi.org/10.1016/j.stem.2021.04.027), adapted by J. Warin with the help of L. David, A. Camus
# qPCR and IF process
Scripts made by A. Humeau
