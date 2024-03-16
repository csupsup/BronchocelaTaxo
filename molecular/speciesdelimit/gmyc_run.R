#Single-locus Species Delimitation using GMYC
#Date: 15 March 2024
#Author: C.E. Supsup

#Generalized Mixed Yule Coalescent (GMYC; Pons et al. 2006)
#Input - ultrametric tree (from Bayesian Inference) in nexus format
#Analysis implemented in R

#1.summarize tree and conver to .nex using treeannotator function from BEAST
cd /Applications/BEAST\ 2.7.4/bin 

treeannotator -burnin 1000 -heights mean Bronchocela_88ND2_BI.trees Bronchocela_88ND2_BI.nex

#2.install required packages in R
install.packages(c("ape", "paran", "rncl"))
install.packages("splits", repos = "http://R-Forge.R-project.org")

#load packages
library(ape)
library(paran)
library(rncl)
library(splits)

#3.read data in R
bronchocela.bi.tree <- read_nexus_phylo("Bronchocela_88ND2_BI.nex")

#3.run GMYC analysis
bronchocela_gmyc <- gmyc(bronchocela.bi.tree, method = "single")

#4.get summary and plot
summary(bronchocela_gmyc)
plot(bronchocela_gmyc)

#5.inspect delimited species
spec.list(bronchocela_gmyc)

#6.plot tree with support values
bronchocela_gmyc_support <- gmyc.support(bronchocela_gmyc) 				# estimate support values
is.na(bronchocela_gmyc_support[bronchocela_gmyc_support == 0]) <- TRUE 	# only show values for affected nodes
plot(bronchocela.bi.tree, cex=.6, no.margin=TRUE)         				 		# plot tree
nodelabels(round(bronchocela_gmyc_support, 2), cex=.7)    				# add support values on tree

#7.export list of delimited species
sp.delim <- spec.list(bronchocela_gmyc)
write.csv(sp.delim, "sp.delim.csv")
