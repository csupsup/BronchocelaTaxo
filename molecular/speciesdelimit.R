#Species delimitation for single-locus using multi-rate Poisson Tree Process (mPTP; Kapli et al. 2017)
#Input - tree (in newick format) from RAxML, PhyML or IQTree analysis
#Analysis implemented in bash terminal (see installation here: https://github.com/Pas-Kapli/mptp)

#1 - perform analysis using maximum likelihood
mptp --tree_file Bronchocela_partition.nex.contree \
	 --minbr_auto Bronchocela_108_trimmed_30April2023.fasta \
	 --output_file Bronchocela_delimitation \
	 --outgroup Ceratophora_stoddartii_WHT_1682,Cophotis_dumbara_GQ_502785,Aphaniotis_fusca_THC_57874,Aphaniotis_fuscus_JAM_1141,Gonocephalus_chamaeleontinus_LSUHC_3788,Gonocephalus_semperi_KU_305748_Cebu_PH,Gonocephalus_semperi_KU_303821_Mindoro_PH,Gonocephalus_semperi_KU_303819_Mindoro_PH,Gonocephalus_semperi_KU_303820_Mindoro_PH

mptp --tree_file Bronchocela_partition.nex.contree \
	 --output_file Bronchocela_delimitation_muti \
	 --ml \
	 --multi \ 
	 --minbr 0.0011674981 \
	 --outgroup Ceratophora_stoddartii_WHT_1682,Cophotis_dumbara_GQ_502785,Aphaniotis_fusca_THC_57874,Aphaniotis_fuscus_JAM_1141,Gonocephalus_chamaeleontinus_LSUHC_3788,Gonocephalus_semperi_KU_305748_Cebu_PH,Gonocephalus_semperi_KU_303821_Mindoro_PH,Gonocephalus_semperi_KU_303819_Mindoro_PH,Gonocephalus_semperi_KU_303820_Mindoro_PH

#2 - perform analysis using MCMC

mptp --tree_file Bronchocela_partition.nex.contree \
     --output_file Bronchocela_delimitation_muti \
     --multi \
     --mcmc 10000000 \
     --mcmc_sample 10000 \
     --mcmc_runs 2 \
     --mcmc_burnin 1000000 \
     --outgroup Ceratophora_stoddartii_WHT_1682,Cophotis_dumbara_GQ_502785,Aphaniotis_fusca_THC_57874,Aphaniotis_fuscus_JAM_1141,Gonocephalus_chamaeleontinus_LSUHC_3788,Gonocephalus_semperi_KU_305748_Cebu_PH,Gonocephalus_semperi_KU_303821_Mindoro_PH,Gonocephalus_semperi_KU_303819_Mindoro_PH,Gonocephalus_semperi_KU_303820_Mindoro_PH 
