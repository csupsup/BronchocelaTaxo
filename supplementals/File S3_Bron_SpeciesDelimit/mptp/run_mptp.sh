## Single-locus Species Delimitation using mPTP
## Date: 08 May 2025
## Author: C.E. Supsup

## Multi-rate Poisson Tree Process (mPTP; Kapli et al. 2017)
## Input - tree (in newick format) from RAxML, PhyML or IQTree analysis
## Analysis implemented in bash terminal (see installation here: https://github.com/Pas-Kapli/mptp)

## 1. Find minimum branch length
mptp --tree_file Bron_partition.nex.contree \
	 --minbr_auto Bron_90ND2_02May2025.fasta \
	 --output_file Bron_delimitation \
	 --outgroup Ceratophora_stoddartii_WHT_1682,Cophotis_dumbara_GQ_502785,Aphaniotis_fusca_THC_57874,Gonocephalus_chamaeleontinus_LSUHC_3788,Gonocephalus_semperi_KU_303821_Mindoro_PH,Gonocephalus_semperi_KU_305748_Cebu_PH

## 2. Run with ML
mptp --tree_file Bron_partition.nex.contree \
	 --output_file Bron_delimitation_ml \
	 --ml --multi \
	 --minbr  0.0011636090 \
	 --outgroup Ceratophora_stoddartii_WHT_1682,Cophotis_dumbara_GQ_502785,Aphaniotis_fusca_THC_57874,Gonocephalus_chamaeleontinus_LSUHC_3788,Gonocephalus_semperi_KU_303821_Mindoro_PH,Gonocephalus_semperi_KU_305748_Cebu_PH

## 3. Run with MCMC
mptp --tree_file Bron_partition.nex.contree \
	 --output_file Bron_delimitation_mcmc \
	 --multi \
	 --minbr 0.0011636090 \
	 --mcmc 10000000 \
	 --mcmc_runs 2 \
	 --mcmc_burnin 1000000 \
	 --outgroup Ceratophora_stoddartii_WHT_1682,Cophotis_dumbara_GQ_502785,Aphaniotis_fusca_THC_57874,Gonocephalus_chamaeleontinus_LSUHC_3788,Gonocephalus_semperi_KU_303821_Mindoro_PH,Gonocephalus_semperi_KU_305748_Cebu_PH