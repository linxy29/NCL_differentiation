# NLC_differentiation

This GitHub repository contains all the scripts and main files that were used for the paper .... analyses and figures. Additonal data and processed files can be available upon reasonable request to the corresponding authors.

# CellChat
- CellChat_NOTO_SB_iPS_NLC_trimean.R: Creation and analyse of a cellchat object for the NOTO SB condition with trimean method
- CellChat_NOTO_iPS_NLC_trimean.R: Creation and analyse of a cellchat object for the NOTO condition with trimean method
- CellChat_NT_iPS_NLC_trimean.R: Creation and analyse of a cellchat object for the NT condition with trimean method
- CellChat_iPS_NLC_trimean.R: Creation and analyse of a cellchat object for iPS-NLC dataset with trimean method
- CellChat_netVisual_chord_bug_correction.R: Retrieved June 30, 2023, from https://github.com/sqjin/CellChat/issues/266
- CellChat_tables_LR.R: Obtain a table of all the receptor ligands of the signalling pathways: COLLAGEN, HH, SEMA3 and SPP1.The tables come from the entire CellChat database, NT, NOTO and NOTO SB conditions.
- CellChat_visualisation_LR_interaction_chord_diagram.R: Visualisation of ligands emitted by the notochord towards other clusters in the signalling pathways : COLLAGEN, HH, SEMA3 and SPP1, under NT, NOTO and NOTO SB conditions.
- CellChat_visualisation_interaction_chord_diagram.R: Visualisation of interactions inside signalling pathways : COLLAGEN, HH, SEMA3 and SPP1, under NT, NOTO and NOTO SB conditions
  
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
- Preprocess_ESC-NLCdataset_and_NotochordWeek8.R: script made to process ESC-NLC and NotochordWeek8 dataset. Contains filtering, normalization, CellCycle regression and differential gene analysis (V. Tam)
- Primary analysis_iPS-NLC.R: script made to process iPS-NLC dataset. Contains filtering, normalization, and differential gene analysis (F. Riemers, J. Warin, A. Camus)
- functions01.R: list of functions used in "Primary analysis_iPS-NLC." for automatisation of anlayses
- iPS-NLC_demultiplex.R: demultiplexing of iPS-NLC dataset before primary analyses using Seurat. Contains the tag assignment, multiplet removal(F. Riemers)
  
Scripts made by V. Tam, F. Riemers, J. Warin, A. Camus
# SCENIC
SCENIC analysis was performed as recommended by the authors (Aibar et al. 2017, https://doi.org/10.1038%2Fnmeth.4463)

- SCENIC_regulons_module_brown_annotatedNO.xlsx: xls file containing all the regulons associated with the WGCNA module specific to the notochord cluster
- allTFs_hg38.txt: transcription factor list based on human genome
  
Scripts made by J. Warin, L. David and A. Camus
# Velocity
Velocity was performed using scVelo and veloAE packages

- iPS-NLC_clusterAXP_VelocityGenes.png: panel of the TOP100 genes driving the velocity arrows arising from AXP cluster
- iPS-NLC_clusterLPXP_VelocityGenes.png: panel of the TOP100 genes driving the velocity arrows arising from LPXP cluster
- iPS-NLC_clusterNO_VelocityGenes.png: panel of the TOP100 genes driving the velocity arrows arising from NO cluster
  
Scripts made by X. Lin
# VennDiagram
- VennDiagram_script: script used to compare gene lists and draw VennDiagram figure
- genelist_human_AxialSkeletonWeek7_NotochordWeek8_vivo.txt: addition of all positive DEG from TBXT+ cell cluster from Zhou et al. 2023, and from notochord cluster in NotochordWeek8 dataset (this study)
- genelist_iPS_ESC-NLC_vitro.txt: all positive DEG from notochord cluster in iPS/ESC-NLC dataset
- genelist_mouse_Tamplin_Peck_Wymeersch.txt: list of mouse notochord genes extracted from Tamplin et al. 2008 (https://doi.org/10.1186/1471-2164-9-511), Tamplin et al. 2011 (https://doi.org/10.1016/j.ydbio.2011.10.002), Peck et al. 2017 (https://doi.org/10.1038/s41598-017-10692-5), Wymeersch et al. 2019 (https://doi.org/10.1242/dev.168161) 
- intersection_human_mouse_vivo.csv: list of genes related to notochord that are shared between mouse and fetal human
- intersection_human_vivo_vitro.csv: list of genes related to notochord that are shared between fetal human an in vitro notochordal cells
- intersection_to_all_vitro_vivo_human_mouse.csv: list of genes related to notochord that are shared between mouse, fetal human and in vitro notochordal cells
- intersection_vitro_mouse.csv: list of genes related to notochord that are shared between mouse and in vitro notochordal cells
- unique_AxialSkeletonWeek7_NotochordWeek8_vivo.csv: list of genes unique to human fetal notochordal cells
- unique_iPS_ESC-NLC.csv: list of genes unique to in vitro notochordal-like cells
- unique_mouse_Tamplin_Peck_Wymeersch.csv: list of genes unique to mouse notochordal cells
  
Scripts made by J. Warin and A. Camus
# WGCNA
TO DO: add annotation and update name of gene lists for modules described on paper

Adapted from Meistermann et al. 2021 (https://doi.org/10.1016/j.stem.2021.04.027), by J. Warin with the help of L. David, A. Camus
# qPCR and IF process
Scripts made by A. Humeau
