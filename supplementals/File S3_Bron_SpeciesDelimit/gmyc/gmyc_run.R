## Single-locus Species Delimitation using GMYC
## Date: 08 May 2025
## Author: C.E. Supsup

## Generalized Mixed Yule Coalescent (GMYC; Pons et al. 2006)
## Input - ultrametric tree (from Bayesian Inference) in nexus format
## Analysis implemented in R

## 1.Summarize tree and convert to .nex using treeannotator function from BEAST
cd /Applications/BEAST\ 2.7.7/bin 

treeannotator -burnin 10 -height mean Bron_90ND2.trees Bron_90ND2_BITree.nex

## 2. Install required packages in R
install.packages(c("ape", "paran", "rncl"))
install.packages("splits", repos = "http://R-Forge.R-project.org")

## Load packages
library(ape)
library(paran)
library(rncl)
library(splits)

## 3. Read data in R
bron.bi.tree <- read_nexus_phylo("Bron_90ND2_BITree.nex")

## 4. Run GMYC analysis
bron_gmyc <- gmyc(bron.bi.tree, method = "multiple")

## 5. Get summary and plot
summary(bron_gmyc)
plot(bron_gmyc)

## 6. Inspect delimited species
spec.list(bron_gmyc)

## 7. Plot tree with support values
bron_gmyc_support <- gmyc.support(bron_gmyc) 				
is.na(bron_gmyc_support[bron_gmyc_support == 0]) <- TRUE 	
plot(bron.bi.tree, cex=.6, no.margin=TRUE)         				
nodelabels(round(bron_gmyc_support, 2), cex=.7)    				

## 8. Export list of delimited species
sp.delim <- spec.list(bron_gmyc)
write.csv(sp.delim, "sp.delim_multi.csv")
