---
title: "Bronchocela Morphology - Univariate Analysis"
author: "Christian Supsup"
date: "24 July 2023"
output:
  html_document:
    df_print: paged
  pdf_document: default
indent: no
---
```{r setup, include=FALSE, message=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
```
**1. Read data**
```{r read data, results='hide'}
library(ggplot2)
library(GroupStruct)
library(ggstatsplot)
library(dplyr)
library(tidyverse)
library(ggbiplot)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(readr)
library(PupillometryR)
library(corrplot)
library(broom)
library(rstatix)
library(psych)

broncho.data <- read.csv("Bronchocela_23July2023_univar.csv")
broncho.mens <- subset(broncho.data, select = c(1, 3:16)) #get mensural data
broncho.mens.trans <- allom(broncho.mens, type = "population1") #perform allometry correction
broncho.mer <- subset(broncho.data, select = c(1, 17:26)) #get meristic data
broncho.data.mod <- cbind(broncho.mens.trans, broncho.mer[2:11])
```
```{r, veiw data}
head(broncho.data.mod)

#summary stats
#tapply(broncho.data$SVL, broncho.data$Pop, summary)
#use psych package
#svl.sum <- describeBy(broncho.data$SVL, group = broncho.data$Pop, mat = TRUE)

#loop through all characters
sum.res <- data.frame()
for(col in broncho.data[,3:28]){
  res <- describeBy(col, group = broncho.data$Pop, mat = TRUE)
  sum.res <- rbind(sum.res,res)}

#rename item names
m.name <- colnames(broncho.data[3:ncol(broncho.data)])
pop <- length(unique(broncho.data$Pop)) # pop is the number of population (i.e., Pop column)
m.name <-  as.data.frame(rep(c(m.name), each = paste(pop))) 
colnames(m.name)[1] <- "char"
sum.res["item"] <- m.name

head(sum.res)

#export summary stats
write.csv(sum.res, "sum.res.csv")
```

**2. Perform normality test**
```{r norm test}
#shapiro-wilk test for mensural data
mens.norm.res <- broncho.mens.trans %>%
  group_by(Pop) %>%
  summarise_all(.funs = funs(statistic = shapiro.test(.)$statistic, 
                             p.value = shapiro.test(.)$p.value))
  head(mens.norm.res)

#shapiro-wilk test for meristic data
broncho.mer.sub <- subset(broncho.mer, select = c(1:7,10,11)) #subset data to exclude DorSPB & DorSPU - to analyze separately.
mer.norm.res <- broncho.mer.sub %>%
  group_by(Pop) %>%
  summarise_all(.funs = funs(statistic = shapiro.test(.)$statistic, 
                             p.value = shapiro.test(.)$p.value))
  head(mer.norm.res)

#combine mensural and meristic data results
norm.st.pval.comb <- cbind(mens.norm.res[,1:15], mer.norm.res[,2:9], 
                              mens.norm.res[,17:29], mer.norm.res[10:17])

head(norm.st.pval.comb)

#export normality test results
write.csv(norm.st.pval.comb, "Supplementary_File_S3.csv")

#test DorSPB & DorSPU using qqplot
broncho.mer.sub2 <- subset(broncho.mer, select = c(1,8:9))

broncho.mer.sub2 %>%
  ggplot(aes(sample = DorSPB)) +
  geom_qq() + geom_qq_line() +
  facet_wrap(~Pop, scales = "free_y")

broncho.mer.sub2 %>%
  ggplot(aes(sample = DorSPU)) +
  geom_qq() + geom_qq_line() +
  facet_wrap(~Pop, scales = "free_y")

#test for tail
broncho.tail <- read.csv('Bronchocela_23July2023_tail.csv')
broncho.tail <- broncho.tail[3:8]

#perform allometry correction using SVL
broncho.tail.trans <- allom(broncho.tail, type = "population1") 

#get summary
tail.sum <- describeBy(broncho.tail$TL, group = broncho.tail$Pop, mat = TRUE)
head(tail.sum)

#export summary stat
write.csv(tail.sum, "tail.sum.csv")

#perform normality test
tail.norm <- broncho.tail %>%
  group_by(Pop) %>%
  summarise_all(.funs = funs(statistic = shapiro.test(.)$statistic, 
                             p.value = shapiro.test(.)$p.value))
head(tail.norm)

#export normality test result for tail
write.csv(tail.norm, "tail.norm.csv")
```
