#Single-locus Species Delimitation using mPTP
#Date: 15 March 2024
#Author: C.E. Supsup

#Multi-rate Poisson Tree Process (mPTP; Kapli et al. 2017)
#Input - tree (in newick format) from RAxML, PhyML or IQTree analysis
#Analysis implemented in bash terminal (see installation here: https://github.com/Pas-Kapli/mptp)

#1.find minimum branch length
mptp --tree_file Bronchocela_partition.nex.contree \
	 --minbr_auto Bronchocela_88ND2.fasta \
	 --output_file Bronchocela_delimitation \
	 --outgroup Ceratophora_stoddartii_WHT_1682,Cophotis_dumbara_GQ_502785,Aphaniotis_fusca_THC_57874,Gonocephalus_chamaeleontinus_LSUHC_3788,Gonocephalus_semperi_KU_303821_Mindoro_PH,Gonocephalus_semperi_KU_305748_Cebu_PH

#2.run with ML
mptp --tree_file Bronchocela_partition.nex.contree \
	 --output_file Bronchocela_delimitation_ml \
	 --ml --multi \
	 --minbr 0.0011677408 \
	 --outgroup Ceratophora_stoddartii_WHT_1682,Cophotis_dumbara_GQ_502785,Aphaniotis_fusca_THC_57874,Gonocephalus_chamaeleontinus_LSUHC_3788,Gonocephalus_semperi_KU_303821_Mindoro_PH,Gonocephalus_semperi_KU_305748_Cebu_PH

#3.run with MCMC
mptp --tree_file Bronchocela_partition.nex.contree \
	 --output_file Bronchocela_delimitation_mcmc \
	 --multi \
	 --minbr 0.0011674981 \
	 --mcmc 10000000 \
	 --mcmc_runs 2 \
	 --mcmc_burnin 1000000 \
	 --outgroup Ceratophora_stoddartii_WHT_1682,Cophotis_dumbara_GQ_502785,Aphaniotis_fusca_THC_57874,Gonocephalus_chamaeleontinus_LSUHC_3788,Gonocephalus_semperi_KU_303821_Mindoro_PH,Gonocephalus_semperi_KU_305748_Cebu_PH