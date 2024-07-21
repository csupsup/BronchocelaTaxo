---
title: "Bronchocela Morphology - Univariate Analysis"
author: "Christian Supsup"
date: "29 March 2024"
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
```{r read data, results='hide', message=FALSE}
library(ggplot2)
library(GroupStruct)
library(ggstatsplot)
library(Rmisc)
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
library(car)
library(Hmisc)
library(ComplexHeatmap)
library(ggstatsplot)

broncho.data.raw <- read.csv("Bronchocela_21March2024_univar.csv")
broncho.data <- subset(broncho.data.raw, select = c(6,12,14:37))
names(broncho.data)[1] <- "Pop"

head(broncho.data)
```

**2. Get data summary**
```{r data summary}
#summary stats
#tapply(broncho.data$SVL, broncho.data$Pop, summary)
#use psych package
#svl.sum <- describeBy(broncho.data$SVL, group = broncho.data$Pop, mat = TRUE)

#loop through all characters
sum.res <- data.frame()
for(col in broncho.data[,3:26]){
  res <- describeBy(col, group = broncho.data$Pop, mat = TRUE)
  sum.res <- rbind(sum.res,res)}

#rename item names
m.name <- colnames(broncho.data[3:ncol(broncho.data)])
pop <- length(unique(broncho.data$Pop)) # pop is the number of population (i.e., "Pop" column)
m.name <-  as.data.frame(rep(c(m.name), each = paste(pop))) 
colnames(m.name)[1] <- "char"
sum.res["item"] <- m.name

head(sum.res)

#export summary stats
write.csv(sum.res, "pooled_data_summary_stats.csv")
```

**3. Test sexual dimorphism with combined datasets**
```{r test.pooled.dimorph, results = "hide"}
#read data
broncho.ancova.raw <- read.csv("Bronchocela_21March2024_univar.csv")
broncho.ancova.sub <- subset(broncho.ancova.raw, select = c(6, 10,12:36))
names(broncho.ancova.sub)[1] <- "Pop"
#count samples and get equal sample for each sex
broncho.ancova.sub %>% dplyr::count(Sex)
set.seed(122588)
broncho.ancova.sub <- broncho.ancova.sub %>% group_by(Sex) %>% slice_sample(n=99)
broncho.ancova.sub <- broncho.ancova.sub [4:18]
write.csv(broncho.ancova.sub, "broncho.ancova.sub.csv")
broncho.ancova.sub <- read.csv("broncho.ancova.sub.csv")
broncho.ancova.sub <- allom(broncho.ancova.sub[2:16], type = "population1")
```
```{r check norm}
#check normality
test.norm <- broncho.ancova.sub %>%
  group_by(Sex) %>%
  summarise_all(.funs = funs(statistic = shapiro.test(.)$statistic, 
                             p.value = shapiro.test(.)$p.value))
test.norm

write.csv(test.norm, "norm.test.res.bysex.csv")

comb.test.norm <- do.call(rbind, lapply(broncho.ancova.sub[2:15], 
                                        function(x) shapiro.test(x)[c("statistic", "p.value")]))

comb.test.norm

write.csv(comb.test.norm, "norm.test.res.comb.csv")

#drop normal characters
broncho.ntnorm <- broncho.ancova.sub[, !(colnames(broncho.ancova.sub) %in% c("SVL", "ED", "OD", "HW", "HD", 
  "HL", "HdL", "FbL", "HbL", "SnsW"))]
#run paired wilcoxon test
FtL.wt <- broncho.ntnorm %>% wilcox_test (FtL ~ Sex, paired = TRUE) %>%
                   adjust_pvalue(method = "bonferroni") %>%
                   add_significance("p.adj")
#loop wilcox test
results <- NULL
results <- as.data.frame(results)
for(col in 2:ncol(broncho.ntnorm))
  {results[col-1,1] <- names(broncho.ntnorm)[col]
  wt.res <- broncho.ntnorm %>% wilcox_test(formula(paste0(names(broncho.ntnorm)[col],"~","Sex")), paired=TRUE) %>%
            adjust_pvalue(method = "bonferroni") %>%
            add_significance("p.adj")
    results[col-1,2] <- wt.res$statistic[1]
    results[col-1,3] <- wt.res$p.adj[1]
    results[col-1,4] <- wt.res$p.adj.signif[1]}
colnames(results)[2] <- "statistic"
colnames(results)[3] <- "p.adj"
colnames(results)[4] <- "p.adj.signif"

View(results)
write.csv(results, "wilcox.test.results.csv")

#run t-test
#drop non-normal characters
broncho.norm <- broncho.ancova.sub[, !(colnames(broncho.ancova.sub) %in% c("FtL", "FnsH", "SnsH", "FnsW"))]
#single run
broncho.norm %>% t_test (ED ~ Sex, paired = FALSE)  %>%
            adjust_pvalue(method = "bonferroni") %>%
            add_significance("p.adj")
#loop paired t-test
results <- NULL
results <- as.data.frame(results)
for(col in 2:ncol(broncho.norm))
  {results[col-1,1] <- names(broncho.norm)[col]
  ttest.res <- broncho.norm %>% t_test(formula(paste0(names(broncho.norm)[col],"~","Sex")), paired=FALSE) %>%
            adjust_pvalue(method = "bonferroni") %>%
            add_significance("p.adj")
    results[col-1,2] <- ttest.res$statistic[1]
    results[col-1,3] <- ttest.res$df[1]
    results[col-1,4] <- ttest.res$p.adj[1]
    results[col-1,5] <- ttest.res$p.adj.signif[1]}
colnames(results)[2] <- "statistic"
colnames(results)[3] <- "df"
colnames(results)[4] <- "p.adj"
colnames(results)[5] <- "p.adj.signif"

View(results)
write.csv(results, "t.test.results.csv")

#create plots for characters with significant differences
plot_color <- c("#bdd7e7", "#2171b5")
my_comparisons = list( c("f", "m"))
#show boxplot for SVL
ggplot(broncho.ancova.sub, aes(x = Sex, y = SVL, fill = Sex)) + geom_boxplot() +
      theme_bw() +
      theme(axis.text=element_text(size=15, color = "black"), 
            axis.title=element_text(size=15, face="bold")) +
      ylab("Snout-Vent Length (mm)") + 
      xlab("") + 
      scale_fill_manual(name = "", labels = c("Female", "Male"), values = plot_color) + 
      stat_compare_means(comparisons = my_comparisons, label.y = c(2.10), 
                         method = "t.test", label = "p.signif")
```

**4. Test sexual dimporphism for each population**
```{r test.pop.dimorph, results = "hide"}
#subset data
broncho.data.sub <- subset(broncho.data.raw, select = c(6, 13:37))
colnames(broncho.data.sub)[1] <- "Pop"
broncho.mens <- subset(broncho.data, select = c(1, 3:16)) #get mensural data
broncho.mens.trans <- allom(broncho.mens, type = "population1") #perform allometry correction for each population
broncho.mens.trans <- cbind(broncho.data.sub$Sex, broncho.mens.trans)
colnames(broncho.mens.trans)[1] <- "Sex"
```
```{r check dimorph}
#check sexual dimorphism with two-way ANOVA using SVL
broncho.data.2aov <- broncho.mens.trans
broncho.data.2aov$Pop <- factor(broncho.data.2aov$Pop)
broncho.data.2aov$Sex <- factor(broncho.data.2aov$Sex)
broncho.data.2aov <- broncho.data.2aov[-(which(broncho.data.2aov$Pop %in% "Lubang")),]
broncho.data.2aov <- broncho.data.2aov[-(which(broncho.data.2aov$Pop %in% "Romblon")),]
write.csv(broncho.data.2aov, "broncho_sex_data.csv")

model <- aov(SVL ~ Sex + Pop + Sex:Pop, data = broncho.data.2aov, type = "III")
summary(model)

#plot error bars
broncho.data.2aov <- broncho.mens.trans
broncho.data.2aov$Pop <- factor(broncho.data.2aov$Pop)
broncho.data.2aov$Sex <- factor(broncho.data.2aov$Sex)

svl.sum <- summarySE(broncho.data.2aov, measurevar="SVL", groupvars=c("Sex","Pop"))

svl.error <- ggplot(svl.sum, aes(x=Pop, y=SVL, colour=Sex, group = Sex)) + 
              theme_bw() +
              geom_errorbar(aes(ymin=SVL-se, ymax=SVL+se), width=.1) +
              geom_line(linewidth = 1) +
              geom_point(size = 3) +
              scale_color_manual(name = "",labels = c("Female", "Male"), 
                  values = c("#bdd7e7", "#2171b5")) +
              ylab("Snout-Vent Length (mm)") +
              xlab("") +
              scale_shape_manual(name = "", labels = c("Female", "Male"), values = c(15,16)) +
              theme(axis.text = element_text(size = 15, color = "black"), 
                      axis.title = element_text(size = 15, face ="bold"), 
                      legend.text = element_text(size = 15, face = "bold"),
                      plot.title = element_text(size = 17, face = "bold"),
                      legend.position = c(0.90, 0.95),
                      legend.background = element_blank())

#export result
png("svl_error.png", res=300, width = 11, height = 6, unit="in") 
svl.error
dev.off()

pdf("svl_error.pdf", width = 11, height = 6) 
svl.error
dev.off()

knitr::include_graphics("svl_error.png")
```

**4. Test character differences with female data**
```{r female.char.diff, results='hide'}
#create directory for outputs
dir.create("tk_female")
dir.create("dn_female")
#get female data
broncho.data.f <- broncho.data.sub[broncho.data.sub$Sex=="f",]

#data summary
sum.res <- data.frame()
for(col in broncho.data.f[,3:26]){
  res <- describeBy(col, group = broncho.data.f$Pop, mat = TRUE)
  sum.res <- rbind(sum.res,res)}
#rename item names
m.name <- colnames(broncho.data.f[3:ncol(broncho.data.f)])
pop <- length(unique(broncho.data.f$Pop)) # pop is the number of population (i.e., "Pop" column)
m.name <-  as.data.frame(rep(c(m.name), each = paste(pop))) 
colnames(m.name)[1] <- "char"
sum.res["item"] <- m.name

head(sum.res)

#export summary stats
write.csv(sum.res, "female_data_summary_stats.csv")

broncho.f.trans <- allom(subset(broncho.data.f, select = c(1, 3:16)), type = "population1") #perform allometry correction for each population
broncho.f.norm <- cbind(broncho.f.trans[1:13], subset(broncho.data.f, select = (24))) #characters only with normal residuals for ANOVA
#data for plot
broncho.female <- cbind(broncho.f.trans, subset(broncho.data.f, select = c(17:26)))
#check data distribution
hist.data.frame(broncho.f.norm [2:13])

#perform anova on a single character
aov.mod <- aov(SVL ~ Pop, data = broncho.f.norm)
summary(aov.mod)
#check residual
par(mfrow = c(2,2))
plot(aov.mod)
#perform shapiro-wilk test on residuals
aov_residuals <- residuals(object = aov.mod)
shapiro.test(x = aov_residuals)

#loop through all characters
results <- NULL
results <- as.data.frame(results)
for(col in 2:ncol(broncho.f.norm ))
  {results[col-1,1] <- names(broncho.f.norm )[col]
    aov.res <- broncho.f.norm  %>% anova_test(formula(paste0(names(broncho.f.norm )[col],"~","Pop"))) %>%
              adjust_pvalue(method = "bonferroni") %>%
              add_significance("p.adj")
    results[col-1,2] <- aov.res$F[1]
    results[col-1,3] <- aov.res$p[1]
    results[col-1,4] <- aov.res$p.adj[1]
    results[col-1,5] <- aov.res$p.adj.signif[1]}
colnames(results)[2] <- "F"
colnames(results)[3] <- "pvalue"
colnames(results)[4] <- "p.adj"
colnames(results)[5] <- "p.adj.signif"

write.csv(results, "aov_female.csv")

#using r base function
aov.res <- apply(broncho.f.norm [,2:ncol(broncho.f.norm )], 2, function(x) aov(x ~ Pop, data = broncho.f.norm ))
#combine anova results
aov.res.df <- do.call(rbind, lapply(aov.res, broom::tidy))
#tukey post-hoc test
tukey.res <- sapply(aov.res, function(x) TukeyHSD(x, "Pop", ordered = TRUE))
#combine tukey results
tukey.res.df <- as.data.frame(do.call(rbind, Map(cbind, Name = names(tukey.res), tukey.res)))

#perform tukey hsd manually
SVL.tk <- aov(SVL ~ Pop, data = broncho.f.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(SVL.tk, "tk_female/tk.SVL.csv")
ED.tk <- aov(ED ~ Pop, data = broncho.f.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(ED.tk, "tk_female/tk.ED.csv")
OD.tk <- aov(OD ~ Pop, data = broncho.f.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(OD.tk, "tk_female/tk.OD.csv")
HW.tk <- aov(HW ~ Pop, data = broncho.f.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(HW.tk, "tk_female/tk.HW.csv")
HD.tk <- aov(HD ~ Pop, data = broncho.f.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(HD.tk, "tk_female/tk.HD.csv")
HL.tk <- aov(HL ~ Pop, data = broncho.f.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(HL.tk, "tk_female/tk.HL.csv")
HdL.tk <- aov(HdL ~ Pop, data = broncho.f.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(HdL.tk, "tk_female/tk.HdL.csv")
FbL.tk <- aov(FbL ~ Pop, data = broncho.f.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(FbL.tk, "tk_female/tk.FbL.csv")
FtL.tk <- aov(FtL ~ Pop, data = broncho.f.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(FtL.tk, "tk_female/tk.FtL.csv")
HbL.tk <- aov(HbL ~ Pop, data = broncho.f.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(HbL.tk, "tk_female/tk.HbL.csv")
FnsH.tk <- aov(FnsH ~ Pop, data = broncho.f.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(FnsH.tk, "tk_female/tk.FnsH.csv")
SnsH.tk <- aov(SnsH ~ Pop, data = broncho.f.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(SnsH.tk, "tk_female/tk.SnsH.csv")
SlT4.tk <- aov(SlT4 ~ Pop, data = broncho.f.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(SlT4.tk, "tk_female/tk.SlT4.csv")

#perform kruskal-wallis on non-normal characters
broncho.f.non.norm <- cbind(subset(broncho.f.trans, select = c(1, 14:15)), subset(broncho.data.f, select = c(17:23, 25, 26)))

#loop through all characters
results <- NULL
results <- as.data.frame(results)
for(col in 2:ncol(broncho.f.non.norm))
  { results[col-1,1] <- names(broncho.f.non.norm)[col]
    kw.res <- broncho.f.non.norm %>% kruskal_test(formula(paste0(names(broncho.f.non.norm)[col],"~","Pop"))) %>%
              adjust_pvalue(method = "bonferroni") %>%
              add_significance("p.adj")
    results[col-1,2] <- kw.res$statistic[1]
    results[col-1,3] <- kw.res$p[1]
    results[col-1,4] <- kw.res$p.adj.signif[1]}
colnames(results)[2] <- "F"
colnames(results)[3] <- "pvalue"
colnames(results)[3] <- "pvalue.signif"

write.csv(results, 'kw_female.csv')

#perform dunn test
FnsW <- broncho.f.non.norm %>% dunn_test(FnsW ~ Pop, p.adjust.method = "bonferroni")
write.csv(FnsW, "dn_female/dn.FnsW.csv")
SnsW <- broncho.f.non.norm %>% dunn_test(SnsW ~ Pop, p.adjust.method = "bonferroni")
write.csv(SnsW, "dn_female/dn.SnsW.csv")
MidSR <- broncho.f.non.norm %>% dunn_test(MidSR ~ Pop, p.adjust.method = "bonferroni")
write.csv(MidSR, "dn_female/dn.MidSR.csv")
VerSC <- broncho.f.non.norm %>% dunn_test(VerSC ~ Pop, p.adjust.method = "bonferroni")
write.csv(VerSC, "dn_female/dn.VerSC.csv")
PoPV <- broncho.f.non.norm %>% dunn_test(PoPV ~ Pop, p.adjust.method = "bonferroni")
write.csv(PoPV, "dn_female/dn.PoPV.csv")
UpPV <- broncho.f.non.norm %>% dunn_test(UpPV ~ Pop, p.adjust.method = "bonferroni")
write.csv(UpPV, "dn_female/dn.UpPV.csv")
SL <- broncho.f.non.norm %>% dunn_test(SL ~ Pop, p.adjust.method = "bonferroni")
write.csv(SL, "dn_female/dn.SL.csv")
IL <- broncho.f.non.norm %>% dunn_test(IL ~ Pop, p.adjust.method = "bonferroni")
write.csv(IL, "dn_female/dn.IL.csv")
SlF3 <- broncho.f.non.norm %>% dunn_test(SlF3 ~ Pop, p.adjust.method = "bonferroni")
write.csv(SlF3, "dn_female/dn.SlF3.csv")
NS <- broncho.f.non.norm %>% dunn_test(NS ~ Pop, p.adjust.method = "bonferroni")
write.csv(NS, "dn_female/dn.NS.csv")
CS <- broncho.f.non.norm %>% dunn_test(CS ~ Pop, p.adjust.method = "bonferroni")
write.csv(CS, "dn_female/dn.CS.csv")

#get p-values from tukey test
tk.tests <- list.files("tk_female/", pattern = 'tk.', full.names = TRUE)
tk.pval <- as.data.frame(lapply(tk.tests, function(x) read.csv(x)$p.adj))
names(tk.pval) <- tools::file_path_sans_ext(basename(tk.tests))
#fix column names
names(tk.pval) <- sub("tk.", "", names(tk.pval))

#get p-values from dunn test
dn.tests <- list.files("dn_female/", pattern = 'dn.', full.names = TRUE)
dn.pval <- as.data.frame(lapply(dn.tests, function(x) read.csv(x)$p.adj))
names(dn.pval) <- tools::file_path_sans_ext(basename(dn.tests))
#fix column names
names(dn.pval) <- sub("dn.", "", names(dn.pval))
#get row labels
#tk.pval$Comparison <- paste(ED.tk$group1,"-",ED.tk$group2)
#get row labels
row.labs <- paste(ED.tk$group1,"-",ED.tk$group2)

#combine tukey and dunn test results
phoc.pval <- cbind(tk.pval, dn.pval)
phoc.pval <- phoc.pval %>%
              mutate_all(as.numeric)
phoc.pval <- cbind(row.labs, phoc.pval)
female.phoc.pval <- phoc.pval %>% tibble::column_to_rownames('row.labs')
phcol <- colnames(broncho.data.sub[3:26])
female.phoc.pval <- female.phoc.pval[phcol]

write.csv(female.phoc.pval, "phoc_pval_female.csv")

#tail stats
broncho.tail <- read.csv('Bronchocela_21March2024_tail.csv')
broncho.tail.f <- broncho.tail[broncho.tail$Sex=="f",]
broncho.tail.f <- subset(broncho.tail.f, select = c(6,12:20))
names(broncho.tail.f)[1] <- "Pop"

#perform allometry correction using SVL
broncho.tail.f.trans <- allom(broncho.tail.f, type = "population1") 
#get summary
tail.sum <- describeBy(broncho.tail.f$TL, group = broncho.tail.f$Pop, mat = TRUE)
head(tail.sum)
#export summary stat
write.csv(tail.sum, "tail_female.sum.csv")

#perform anova on a single character
aov.mod <- aov(TL ~ Pop, data = broncho.tail.f)
summary(aov.mod)
#check residual
par(mfrow = c(2,2))
plot(aov.mod)
#perform shapiro-wilk test on residuals
aov_residuals <- residuals(object = aov.mod)
shapiro.test(x = aov_residuals)

TL.tk <- aov(TL ~ Pop, data = broncho.tail.f ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(TL.tk, "tk.TL_female.csv")
```
```{r, female pheatmap}
#plot heatmap
f.phoc.map <- ComplexHeatmap::pheatmap(female.phoc.pval, cluster_rows = F, cluster_cols = F,
                         display_numbers = T, #cutree_rows = 4, 
                         color = c("#6baed6", "#eff3ff"), breaks = c(0.000,0.04,1.0), 
                         border_color = "black", number_color = "black",
                         fontsize_number = 12, fontsize = 15, fontface = "bold", column_names_side = c("top"),
                         angle_col = c("0"), heatmap_legend_param = list(title = "p-value", labels_gp = gpar(fontsize = 15), 
                         title_gp = gpar(fontsize = 16, fontface = "bold")), cellwidth = 50, cellheight = 17)

#export result
png("female.phoc.pheatmap.png", res=300, width = 21, height = 9, unit="in")
f.phoc.map
dev.off()

pdf("female.phoc.pheatmap.pdf", width = 21, height = 9)
f.phoc.map
dev.off()

knitr::include_graphics("female.phoc.pheatmap.png")
```

**5. Test character differences with male data**
```{r male.char.diff, results='hide'}
#create directory for outputs
dir.create("tk_male")
dir.create("dn_male")
#get female data
broncho.data.sub <- subset(broncho.data.raw, select = c(6, 13:37))
colnames(broncho.data.sub)[1] <- "Pop"
broncho.data.m <- broncho.data.sub[broncho.data.sub$Sex=="m",]

#data summary
sum.res <- data.frame()
for(col in broncho.data.m[,3:26]){
  res <- describeBy(col, group = broncho.data.m$Pop, mat = TRUE)
  sum.res <- rbind(sum.res,res)}
#rename item names
m.name <- colnames(broncho.data.m[3:ncol(broncho.data.m)])
pop <- length(unique(broncho.data.m$Pop)) # pop is the number of population (i.e., "Pop" column)
m.name <-  as.data.frame(rep(c(m.name), each = paste(pop))) 
colnames(m.name)[1] <- "char"
sum.res["item"] <- m.name

head(sum.res)

#export summary stats
write.csv(sum.res, "male_data_summary_stats.csv")

broncho.m.trans <- allom(subset(broncho.data.m, select = c(1, 3:16)), type = "population1") #perform allometry correction for each population
broncho.m.norm <- subset(broncho.m.trans, select = c(1, 3:13)) #characters only with normal residuals for ANOVA
#data for plot
broncho.male <- cbind(broncho.m.trans, subset(broncho.data.m, select = c(17:26)))
#check data distribution
broncho.m.norm <- cbind(broncho.m.norm, subset(broncho.data.m, select = c(17,23,24)))

hist.data.frame(broncho.m.norm [2:12])

#perform anova on a single character
aov.mod <- aov(ED ~ Pop, data = broncho.m.norm)
summary(aov.mod)
#check residual
par(mfrow = c(2,2))
plot(aov.mod)
#perform shapiro-wilk test on residuals
aov_residuals <- residuals(object = aov.mod)
shapiro.test(x = aov_residuals)

#loop through all characters
results <- NULL
results <- as.data.frame(results)
for(col in 2:ncol(broncho.m.norm ))
  {results[col-1,1] <- names(broncho.m.norm )[col]
    aov.res <- broncho.m.norm  %>% anova_test(formula(paste0(names(broncho.m.norm )[col],"~","Pop"))) %>%
              adjust_pvalue(method = "bonferroni") %>%
              add_significance("p.adj")
    results[col-1,2] <- aov.res$F[1]
    results[col-1,3] <- aov.res$p[1]
    results[col-1,4] <- aov.res$p.adj[1]
    results[col-1,5] <- aov.res$p.adj.signif[1]}
colnames(results)[2] <- "F"
colnames(results)[3] <- "pvalue"
colnames(results)[4] <- "p.adj"
colnames(results)[5] <- "p.adj.signif"

write.csv(results, "aov_male.csv")

#using r base function
aov.res <- apply(broncho.m.norm [,2:ncol(broncho.m.norm )], 2, function(x) aov(x ~ Pop, data = broncho.m.norm ))
#combine anova results
aov.res.df <- do.call(rbind, lapply(aov.res, broom::tidy))
#tukey post-hoc test
tukey.res <- sapply(aov.res, function(x) TukeyHSD(x, "Pop", ordered = TRUE))
#combine tukey results
tukey.res.df <- as.data.frame(do.call(rbind, Map(cbind, Name = names(tukey.res), tukey.res)))

#perform tukey hsd manually
ED.tk <- aov(ED ~ Pop, data = broncho.m.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(ED.tk, "tk_male/tk.ED.csv")
OD.tk <- aov(OD ~ Pop, data = broncho.m.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(OD.tk, "tk_male/tk.OD.csv")
HW.tk <- aov(HW ~ Pop, data = broncho.m.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(HW.tk, "tk_male/tk.HW.csv")
HD.tk <- aov(HD ~ Pop, data = broncho.m.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(HD.tk, "tk_male/tk.HD.csv")
HL.tk <- aov(HL ~ Pop, data = broncho.m.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(HL.tk, "tk_male/tk.HL.csv")
HdL.tk <- aov(HdL ~ Pop, data = broncho.m.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(HdL.tk, "tk_male/tk.HdL.csv")
FbL.tk <- aov(FbL ~ Pop, data = broncho.m.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(FbL.tk, "tk_male/tk.FbL.csv")
FtL.tk <- aov(FtL ~ Pop, data = broncho.m.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(FtL.tk, "tk_male/tk.FtL.csv")
HbL.tk <- aov(HbL ~ Pop, data = broncho.m.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(HbL.tk, "tk_male/tk.HbL.csv")
FnsH.tk <- aov(FnsH ~ Pop, data = broncho.m.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(FnsH.tk, "tk_male/tk.FnsH.csv")
SnsH.tk <- aov(SnsH ~ Pop, data = broncho.m.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(SnsH.tk, "tk_male/tk.SnsH.csv")
MidSR.tk <- aov(MidSR ~ Pop, data = broncho.m.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(MidSR.tk, "tk_male/tk.MidSR.csv")
SlF3.tk <- aov(SlF3 ~ Pop, data = broncho.m.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(SlF3.tk, "tk_male/tk.SlF3.csv")
SlT4.tk <- aov(SlT4 ~ Pop, data = broncho.m.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(SlT4.tk, "tk_male/tk.SlT4.csv")

#perform kruskal-wallis on non-normal characters
broncho.m.non.norm <- cbind(subset(broncho.m.trans, select = c(1:2, 14,15)), subset(broncho.data.m, select = c(18:22, 25:26)))

#loop through all characters
results <- NULL
results <- as.data.frame(results)
for(col in 2:ncol(broncho.m.non.norm))
  { results[col-1,1] <- names(broncho.m.non.norm)[col]
    kw.res <- broncho.m.non.norm %>% kruskal_test(formula(paste0(names(broncho.m.non.norm)[col],"~","Pop"))) %>%
              adjust_pvalue(method = "bonferroni") %>%
              add_significance("p.adj")
    results[col-1,2] <- kw.res$statistic[1]
    results[col-1,3] <- kw.res$p[1]
    results[col-1,4] <- kw.res$p.adj.signif[1]}
colnames(results)[2] <- "F"
colnames(results)[3] <- "pvalue"
colnames(results)[3] <- "pvalue.signif"

write.csv(results, 'kw_male.csv')

#perform dunn test
SVL <- broncho.m.non.norm %>% dunn_test(SVL ~ Pop, p.adjust.method = "bonferroni")
write.csv(SVL, "dn_male/dn.SVL.csv")
FnsW <- broncho.m.non.norm %>% dunn_test(FnsW ~ Pop, p.adjust.method = "bonferroni")
write.csv(FnsW, "dn_male/dn.FnsW.csv")
SnsW <- broncho.m.non.norm %>% dunn_test(SnsW ~ Pop, p.adjust.method = "bonferroni")
write.csv(SnsW, "dn_male/dn.SnsW.csv")
VerSC <- broncho.m.non.norm %>% dunn_test(VerSC ~ Pop, p.adjust.method = "bonferroni")
write.csv(VerSC, "dn_male/dn.VerSC.csv")
PoPV <- broncho.m.non.norm %>% dunn_test(PoPV ~ Pop, p.adjust.method = "bonferroni")
write.csv(PoPV, "dn_male/dn.PoPV.csv")
UpPV <- broncho.m.non.norm %>% dunn_test(UpPV ~ Pop, p.adjust.method = "bonferroni")
write.csv(UpPV, "dn_male/dn.UpPV.csv")
SL <- broncho.m.non.norm %>% dunn_test(SL ~ Pop, p.adjust.method = "bonferroni")
write.csv(SL, "dn_male/dn.SL.csv")
IL <- broncho.m.non.norm %>% dunn_test(IL ~ Pop, p.adjust.method = "bonferroni")
write.csv(IL, "dn_male/dn.IL.csv")
NS <- broncho.m.non.norm %>% dunn_test(NS ~ Pop, p.adjust.method = "bonferroni")
write.csv(NS, "dn_male/dn.NS.csv")
CS <- broncho.m.non.norm %>% dunn_test(CS ~ Pop, p.adjust.method = "bonferroni")
write.csv(CS, "dn_male/dn.CS.csv")

#get p-values from tukey test
tk.tests <- list.files("tk_male/", pattern = 'tk.', full.names = TRUE)
tk.pval <- as.data.frame(lapply(tk.tests, function(x) read.csv(x)$p.adj))
names(tk.pval) <- tools::file_path_sans_ext(basename(tk.tests))
#fix column names
names(tk.pval) <- sub("tk.", "", names(tk.pval))

#get p-values from dunn test
dn.tests <- list.files("dn_male/", pattern = 'dn.', full.names = TRUE)
dn.pval <- as.data.frame(lapply(dn.tests, function(x) read.csv(x)$p.adj))
names(dn.pval) <- tools::file_path_sans_ext(basename(dn.tests))
#fix column names
names(dn.pval) <- sub("dn.", "", names(dn.pval))
#get row labels
#tk.pval$Comparison <- paste(ED.tk$group1,"-",ED.tk$group2)
#get row labels
row.labs <- paste(ED.tk$group1,"-",ED.tk$group2)

#combine tukey and dunn test results
phoc.pval <- cbind(tk.pval, dn.pval)
phoc.pval <- phoc.pval %>%
              mutate_all(as.numeric)
phoc.pval <- cbind(row.labs, phoc.pval)
male.phoc.pval <- phoc.pval %>% tibble::column_to_rownames('row.labs')
phcol <- colnames(broncho.data.sub[3:26])
male.phoc.pval <- male.phoc.pval[phcol]

write.csv(male.phoc.pval, "phoc_pval_male.csv")

#tail stats
broncho.tail <- read.csv('Bronchocela_21March2024_tail.csv')
broncho.tail.m <- broncho.tail[broncho.tail$Sex=="m",]
broncho.tail.m <- subset(broncho.tail.m, select = c(6,12:20))
names(broncho.tail.m)[1] <- "Pop"

#perform allometry correction using SVL
broncho.tail.m.trans <- allom(broncho.tail.m, type = "population1") 
#get summary
tail.sum <- describeBy(broncho.tail.m$TL, group = broncho.tail.m$Pop, mat = TRUE)
head(tail.sum)
#export summary stat
write.csv(tail.sum, "tail_male.sum.csv")

#perform anova on a single character
aov.mod <- aov(TL ~ Pop, data = broncho.tail.m.trans)
summary(aov.mod)
#check residual
par(mfrow = c(2,2))
plot(aov.mod)
#perform shapiro-wilk test on residuals
aov_residuals <- residuals(object = aov.mod)
shapiro.test(x = aov_residuals)

#kw test
broncho.tail.m.trans %>% kruskal_test(TL ~ Pop)

TL <- broncho.tail.m.trans %>% dunn_test(TL ~ Pop, p.adjust.method = "bonferroni")
write.csv(TL, "dn.TL_male.csv")
```
```{r, male heatmap}
#plot pheatmap
m.phoc.map <- ComplexHeatmap::pheatmap(male.phoc.pval, cluster_rows = F, cluster_cols = F,
                         display_numbers = T, #cutree_rows = 4, 
                         color = c("#6baed6", "#eff3ff"), breaks = c(0.000,0.04,1.0), 
                         border_color = "black", number_color = "black",
                         fontsize_number = 12, fontsize = 15, fontface = "bold", column_names_side = c("top"),
                         angle_col = c("0"), heatmap_legend_param = list(title = "p-value", labels_gp = gpar(fontsize = 15), 
                         title_gp = gpar(fontsize = 16, fontface = "bold")), cellwidth = 50, cellheight = 17)

#export result
png("male.phoc.pheatmap.png", res=300, width = 20, height = 6, unit="in")
m.phoc.map
dev.off()

pdf("male.phoc.pheatmap.pdf", width = 20, height = 6)
m.phoc.map
dev.off()

knitr::include_graphics("male.phoc.pheatmap.png")
```

**6. Test character differences with combined dataset**
```{r combined.data, results='hide'}
#create directory for outputs
dir.create("tk_pooled")
dir.create("dn_pooled")
#get female data
broncho.data.sub <- subset(broncho.data.raw, select = c(6, 14:37))
colnames(broncho.data.sub)[1] <- "Pop"
broncho.pld.trans <- allom(subset(broncho.data.sub, select = c(1:15)), type = "population1") #perform allometry correction for each population
broncho.pld.norm <- subset(broncho.pld.trans, select = c(1, 3:11, 14,15)) #characters only with normal residuals for ANOVA
#data for plot
broncho.pld <- cbind(broncho.pld.trans, subset(broncho.data.sub, select = c(16:25)))
#check data distribution
broncho.pld.norm <- cbind(broncho.pld.norm, subset(broncho.data.sub, select = c(23)))

hist.data.frame(broncho.pld.norm [2:13])

#perform anova on a single character
aov.mod <- aov(ED ~ Pop, data = broncho.pld.norm)
summary(aov.mod)
#check residual
par(mfrow = c(2,2))
plot(aov.mod)
#perform shapiro-wilk test on residuals
aov_residuals <- residuals(object = aov.mod)
shapiro.test(x = aov_residuals)

#loop through all characters
results <- NULL
results <- as.data.frame(results)
for(col in 2:ncol(broncho.pld.norm ))
  {results[col-1,1] <- names(broncho.pld.norm )[col]
    aov.res <- broncho.pld.norm  %>% anova_test(formula(paste0(names(broncho.pld.norm )[col],"~","Pop"))) %>%
              adjust_pvalue(method = "bonferroni") %>%
              add_significance("p.adj")
    results[col-1,2] <- aov.res$F[1]
    results[col-1,3] <- aov.res$p[1]
    results[col-1,4] <- aov.res$p.adj[1]
    results[col-1,5] <- aov.res$p.adj.signif[1]}
colnames(results)[2] <- "F"
colnames(results)[3] <- "pvalue"
colnames(results)[4] <- "p.adj"
colnames(results)[5] <- "p.adj.signif"

write.csv(results, "aov_pooled.csv")

#using r base function
aov.res <- apply(broncho.pld.norm [,2:ncol(broncho.pld.norm )], 2, function(x) aov(x ~ Pop, data = broncho.pld.norm ))
#combine anova results
aov.res.df <- do.call(rbind, lapply(aov.res, broom::tidy))
#tukey post-hoc test
tukey.res <- sapply(aov.res, function(x) TukeyHSD(x, "Pop", ordered = TRUE))
#combine tukey results
tukey.res.df <- as.data.frame(do.call(rbind, Map(cbind, Name = names(tukey.res), tukey.res)))

#perform tukey hsd manually
ED.tk <- aov(ED ~ Pop, data = broncho.pld.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(ED.tk, "tk_pooled/tk.ED.csv")
OD.tk <- aov(OD ~ Pop, data = broncho.pld.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(OD.tk, "tk_pooled/tk.OD.csv")
HW.tk <- aov(HW ~ Pop, data = broncho.pld.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(HW.tk, "tk_pooled/tk.HW.csv")
HD.tk <- aov(HD ~ Pop, data = broncho.pld.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(HD.tk, "tk_pooled/tk.HD.csv")
HL.tk <- aov(HL ~ Pop, data = broncho.pld.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(HL.tk, "tk_pooled/tk.HL.csv")
HdL.tk <- aov(HdL ~ Pop, data = broncho.pld.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(HdL.tk, "tk_pooled/tk.HdL.csv")
FbL.tk <- aov(FbL ~ Pop, data = broncho.pld.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(FbL.tk, "tk_pooled/tk.FbL.csv")
FtL.tk <- aov(FtL ~ Pop, data = broncho.pld.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(FtL.tk, "tk_pooled/tk.FtL.csv")
HbL.tk <- aov(HbL ~ Pop, data = broncho.pld.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(HbL.tk, "tk_pooled/tk.HbL.csv")
FnsW.tk <- aov(FnsW ~ Pop, data = broncho.pld.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(FnsW.tk, "tk_pooled/tk.FnsW.csv")
SnsW.tk <- aov(SnsW ~ Pop, data = broncho.pld.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(SnsW.tk, "tk_pooled/tk.SnsW.csv")
SlT4.tk <- aov(SlT4 ~ Pop, data = broncho.pld.norm ) %>% tukey_hsd(p.adjust.method = "bonferroni")
write.csv(SlT4.tk, "tk_pooled/tk.SlT4.csv")

#perform kruskal-wallis on non-normal characters
broncho.pld.non.norm <- cbind(subset(broncho.pld.trans, select = c(1:2, 12, 13)), 
  subset(broncho.data.sub, select = c(16:22, 24:25)))

#loop through all characters
results <- NULL
results <- as.data.frame(results)
for(col in 2:ncol(broncho.pld.non.norm ))
  { results[col-1,1] <- names(broncho.pld.non.norm )[col]
    kw.res <- broncho.pld.non.norm  %>% kruskal_test(formula(paste0(names(broncho.pld.non.norm )[col],"~","Pop"))) %>%
              adjust_pvalue(method = "bonferroni") %>%
              add_significance("p.adj")
    results[col-1,2] <- kw.res$statistic[1]
    results[col-1,3] <- kw.res$p[1]
    results[col-1,4] <- kw.res$p.adj.signif[1]}
colnames(results)[2] <- "F"
colnames(results)[3] <- "pvalue"
colnames(results)[3] <- "pvalue.signif"

write.csv(results, 'kw_pooled.csv')

#perform dunn test
SVL <- broncho.pld.non.norm %>% dunn_test(SVL ~ Pop, p.adjust.method = "bonferroni")
write.csv(SVL, "dn_pooled/dn.SVL.csv")
FnsH <- broncho.pld.non.norm %>% dunn_test(FnsH ~ Pop, p.adjust.method = "bonferroni")
write.csv(FnsH, "dn_pooled/dn.FnsH.csv")
SnsH <- broncho.pld.non.norm %>% dunn_test(SnsH ~ Pop, p.adjust.method = "bonferroni")
write.csv(SnsH, "dn_pooled/dn.SnsH.csv")
MidSR <- broncho.pld.non.norm %>% dunn_test(MidSR ~ Pop, p.adjust.method = "bonferroni")
write.csv(MidSR, "dn_pooled/dn.MidSR.csv")
VerSC <- broncho.pld.non.norm %>% dunn_test(VerSC ~ Pop, p.adjust.method = "bonferroni")
write.csv(VerSC, "dn_pooled/dn.VerSC.csv")
PoPV <- broncho.pld.non.norm %>% dunn_test(PoPV ~ Pop, p.adjust.method = "bonferroni")
write.csv(PoPV, "dn_pooled/dn.PoPV.csv")
UpPV <- broncho.pld.non.norm %>% dunn_test(UpPV ~ Pop, p.adjust.method = "bonferroni")
write.csv(UpPV, "dn_pooled/dn.UpPV.csv")
SL <- broncho.pld.non.norm %>% dunn_test(SL ~ Pop, p.adjust.method = "bonferroni")
write.csv(SL, "dn_pooled/dn.SL.csv")
IL <- broncho.pld.non.norm %>% dunn_test(IL ~ Pop, p.adjust.method = "bonferroni")
write.csv(IL, "dn_pooled/dn.IL.csv")
SlF3 <- broncho.pld.non.norm %>% dunn_test(SlF3 ~ Pop, p.adjust.method = "bonferroni")
write.csv(SlF3, "dn_pooled/dn.SlF3.csv")
NS <- broncho.pld.non.norm %>% dunn_test(NS ~ Pop, p.adjust.method = "bonferroni")
write.csv(NS, "dn_pooled/dn.NS.csv")
CS <- broncho.pld.non.norm %>% dunn_test(CS ~ Pop, p.adjust.method = "bonferroni")
write.csv(CS, "dn_pooled/dn.CS.csv")

#get p-values from tukey test
tk.tests <- list.files("tk_pooled/", pattern = 'tk.', full.names = TRUE)
tk.pval <- as.data.frame(lapply(tk.tests, function(x) read.csv(x)$p.adj))
names(tk.pval) <- tools::file_path_sans_ext(basename(tk.tests))
#fix column names
names(tk.pval) <- sub("tk.", "", names(tk.pval))

#get p-values from dunn test
dn.tests <- list.files("dn_pooled/", pattern = 'dn.', full.names = TRUE)
dn.pval <- as.data.frame(lapply(dn.tests, function(x) read.csv(x)$p.adj))
names(dn.pval) <- tools::file_path_sans_ext(basename(dn.tests))
#fix column names
names(dn.pval) <- sub("dn.", "", names(dn.pval))
#get row labels
#tk.pval$Comparison <- paste(ED.tk$group1,"-",ED.tk$group2)
#get row labels
row.labs <- paste(ED.tk$group1,"-",ED.tk$group2)

#combine tukey and dunn test results
phoc.pval <- cbind(tk.pval, dn.pval)
phoc.pval <- phoc.pval %>%
              mutate_all(as.numeric)
phoc.pval <- cbind(row.labs, phoc.pval)
pld.phoc.pval <- phoc.pval %>% tibble::column_to_rownames('row.labs')
phcol <- colnames(broncho.data.sub[2:25])
pld.phoc.pval <- pld.phoc.pval[phcol]

write.csv(pld.phoc.pval, "phoc_pval_pooled.csv")

#tail stats
broncho.tail <- read.csv('Bronchocela_21March2024_tail.csv')
broncho.tail.pld <- subset(broncho.tail, select = c(6,12:20))
names(broncho.tail.pld)[1] <- "Pop"

#perform allometry correction using SVL
broncho.tail.pld.trans <- allom(broncho.tail.pld, type = "population1")
#get summary
tail.sum <- describeBy(broncho.tail$TL, group = broncho.tail$Region, mat = TRUE)
head(tail.sum)
#export summary stat
write.csv(tail.sum, "tail_pooled.sum.csv")

#perform anova on a single character
aov.mod <- aov(TL ~ Pop, data = broncho.tail.pld.trans)
summary(aov.mod)
#check residual
par(mfrow = c(2,2))
plot(aov.mod)
#perform shapiro-wilk test on residuals
aov_residuals <- residuals(object = aov.mod)
shapiro.test(x = aov_residuals)

#kw test
broncho.tail.pld.trans %>% kruskal_test(TL ~ Pop)

TL <- broncho.tail.pld.trans %>% dunn_test(TL ~ Pop, p.adjust.method = "bonferroni")
write.csv(TL, "dn.TL_pooled.csv")
```
```{r, pooled heatmap}
#plot pheatmap
m.phoc.map <- ComplexHeatmap::pheatmap(pld.phoc.pval, cluster_rows = F, cluster_cols = F,
                         display_numbers = T, #cutree_rows = 4, 
                         color = c("#6baed6", "#eff3ff"), breaks = c(0.000,0.04,1.0), 
                         border_color = "black", number_color = "black",
                         fontsize_number = 12, fontsize = 15, fontface = "bold", column_names_side = c("top"),
                         angle_col = c("0"), heatmap_legend_param = list(title = "p-value", labels_gp = gpar(fontsize = 15), 
                         title_gp = gpar(fontsize = 16, fontface = "bold")), cellwidth = 50, cellheight = 17)

#export result
png("pooled.phoc.pheatmap.png", res=300, width = 21, height = 9, unit="in")
m.phoc.map
dev.off()

pdf("pooled.phoc.pheatmap.pdf", width = 21, height = 9)
m.phoc.map
dev.off()

knitr::include_graphics("pooled.phoc.pheatmap.png")
```

**7. Plot selected significant characters**
```{r, plot.diag.char}
broncho.pld$Pop <- factor(broncho.pld$Pop, levels = c("Luzon", "West Visayas", "Romblon", "Lubang", 
  "Mindoro", "Palawan", "Mindanao", "Moluccas", "Borneo"))
broncho.male$Pop <- factor(broncho.male$Pop, levels = c("Luzon", "West Visayas", 
  "Mindoro", "Palawan", "Mindanao", "Moluccas", "Borneo"))
broncho.female$Pop <- factor(broncho.female$Pop, levels = c("Luzon", "West Visayas", "Romblon", "Lubang", 
  "Mindoro", "Palawan", "Mindanao", "Moluccas", "Borneo"))

pop.lab <- c("Luzon", "West Visayas", "Romblon", "Lubang", 
  "Mindoro", "Palawan", "Mindanao", "Moluccas", "Borneo")

plot_color <- c("#7f404a", "#cc1489", "#cc9666", "#408073", "#cc3d3d", 
  "#8c994d", "#7b3dcc", "#d970cb", "#5b5b5b")

plot_color.m <- c("#7f404a", "#cc1489", "#cc3d3d", 
  "#8c994d", "#7b3dcc", "#d970cb", "#5b5b5b")

#plot mensural data
svl.pld <- ggplot(broncho.pld,aes(x=Pop,y=SVL, fill = Pop, color = Pop))+
         theme_bw() +
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = SVL), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab("SVL") + xlab("") + ggtitle("Combined") + 
         scale_color_manual(values = plot_color) + scale_fill_manual(values = plot_color) + 
         theme(axis.text=element_text(size=15, color = "black"), 
                    axis.title=element_text(size=15, face="bold"), 
                    legend.text=element_text(size=15, face = "bold"),
                    plot.title = element_text(size = 17, face = "bold")) + 
         scale_y_continuous(labels = scales::label_number(accuracy = 0.1), limits = c(1.8, 2.10)) +
         coord_flip()
svl.male <- ggplot(broncho.male,aes(x=Pop,y=SVL, fill = Pop, color = Pop))+
         theme_bw() +
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = SVL), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab("SVL") + xlab("") + ggtitle("Male") + 
         scale_color_manual(values = plot_color.m) + scale_fill_manual(values = plot_color.m) + 
         theme(axis.text=element_text(size=15, color = "black"), 
                    axis.title=element_text(size=15, face="bold"), 
                    legend.text=element_text(size=15, face = "bold"),
                    plot.title = element_text(size = 17, face = "bold")) + 
         scale_y_continuous(labels = scales::label_number(accuracy = 0.1), limits = c(1.8, 2.10)) +
         coord_flip()
svl.female <- ggplot(broncho.female,aes(x=Pop,y=SVL, fill = Pop, color = Pop))+
         theme_bw() +
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = SVL), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab("SVL") + xlab("") + ggtitle("Female") + 
         scale_color_manual(values = plot_color) + scale_fill_manual(values = plot_color) + 
         theme(axis.text=element_text(size=15, color = "black"), 
                    axis.title=element_text(size=15, face="bold"), 
                    legend.text=element_text(size=15, face = "bold"),
                    plot.title = element_text(size = 17, face = "bold")) + 
         scale_y_continuous(labels = scales::label_number(accuracy = 0.1), limits = c(1.8, 2.10)) +
         coord_flip()
hbl.pld <- ggplot(broncho.pld,aes(x=Pop,y=HbL, fill = Pop, color = Pop))+
         theme_bw() +
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = HbL), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab("HbL") + xlab("") + ggtitle("") + 
         scale_color_manual(values = plot_color) + scale_fill_manual(values = plot_color) + 
         theme(axis.text=element_text(size=15, color = "black"), 
                    axis.title=element_text(size=15, face="bold"), 
                    legend.text=element_text(size=15, face = "bold"),
                    plot.title = element_text(size = 17, face = "bold")) + 
         scale_y_continuous(labels = scales::label_number(accuracy = 0.1), limits = c(1.81, 1.97)) +
         coord_flip()
hbl.male <- ggplot(broncho.male,aes(x=Pop,y=HbL, fill = Pop, color = Pop))+
         theme_bw() +
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = HbL), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab("HbL") + xlab("") + ggtitle("") + 
         scale_color_manual(values = plot_color.m) + scale_fill_manual(values = plot_color.m) + 
         theme(axis.text=element_text(size=15, color = "black"), 
                    axis.title=element_text(size=15, face="bold"), 
                    legend.text=element_text(size=15, face = "bold"),
                    plot.title = element_text(size = 17, face = "bold")) + 
         scale_y_continuous(labels = scales::label_number(accuracy = 0.1), limits = c(1.81, 1.97)) +
         coord_flip()
hbl.female <- ggplot(broncho.female,aes(x=Pop,y=HbL, fill = Pop, color = Pop))+
         theme_bw() +
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = HbL), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab("HbL") + xlab("") + ggtitle("") + 
         scale_color_manual(values = plot_color) + scale_fill_manual(values = plot_color) + 
         theme(axis.text=element_text(size=15, color = "black"), 
                    axis.title=element_text(size=15, face="bold"), 
                    legend.text=element_text(size=15, face = "bold"),
                    plot.title = element_text(size = 17, face = "bold")) + 
         scale_y_continuous(labels = scales::label_number(accuracy = 0.1), limits = c(1.81, 1.97)) +
         coord_flip()
snsh.pld <- ggplot(broncho.pld,aes(x=Pop,y=SnsH, fill = Pop, color = Pop))+
         theme_bw() +
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = SnsH), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab("SnsH") + xlab("") + ggtitle("") + 
         scale_color_manual(values = plot_color) + scale_fill_manual(values = plot_color) + 
         theme(axis.text=element_text(size=15, color = "black"), 
                    axis.title=element_text(size=15, face="bold"), 
                    legend.text=element_text(size=15, face = "bold"),
                    plot.title = element_text(size = 17, face = "bold")) + 
         scale_y_continuous(labels = scales::label_number(accuracy = 0.1), limits = c(-0.5, 1.0)) +
         coord_flip()
snsh.male <- ggplot(broncho.male,aes(x=Pop,y=SnsH, fill = Pop, color = Pop))+
         theme_bw() +
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = SnsH), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab("SnsH") + xlab("") + ggtitle("") + 
         scale_color_manual(values = plot_color.m) + scale_fill_manual(values = plot_color.m) + 
         theme(axis.text=element_text(size=15, color = "black"), 
                    axis.title=element_text(size=15, face="bold"), 
                    legend.text=element_text(size=15, face = "bold"),
                    plot.title = element_text(size = 17, face = "bold")) + 
         scale_y_continuous(labels = scales::label_number(accuracy = 0.1), limits = c(-0.5, 1.0)) +
         coord_flip()
snsh.female <- ggplot(broncho.female,aes(x=Pop,y=SnsH, fill = Pop, color = Pop))+
         theme_bw() +
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = SnsH), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab("SnsH") + xlab("") + ggtitle("") + 
         scale_color_manual(values = plot_color) + scale_fill_manual(values = plot_color) + 
         theme(axis.text=element_text(size=15, color = "black"), 
                    axis.title=element_text(size=15, face="bold"), 
                    legend.text=element_text(size=15, face = "bold"),
                    plot.title = element_text(size = 17, face = "bold")) + 
         scale_y_continuous(labels = scales::label_number(accuracy = 0.1), limits = c(-0.5, 1.0)) +
         coord_flip()

bx.mens <- ggarrange(svl.pld,svl.male,svl.female,hbl.pld,hbl.male,hbl.female, snsh.pld, snsh.male, snsh.female,
            nrow = 3, ncol=3, legend = "none", common.legend = TRUE) 

#export result
png("mens_boxplot.png", res=300, width = 15, height = 12, unit="in")
bx.mens
dev.off()

pdf("mens_boxplot.pdf", width = 15, height = 12)
bx.mens
dev.off()

knitr::include_graphics("mens_boxplot.png")

#plot meristic data
versc.pld <- ggplot(broncho.pld,aes(x=Pop,y=VerSC, fill = Pop, color = Pop))+
         theme_bw() +
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = VerSC), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab("VerSC") + xlab("") + ggtitle("Combined") + 
         scale_color_manual(values = plot_color) + scale_fill_manual(values = plot_color) + 
         theme(axis.text=element_text(size=15, color = "black"), 
                    axis.title=element_text(size=15, face="bold"), 
                    legend.text=element_text(size=15, face = "bold"),
                    plot.title = element_text(size = 17, face = "bold")) + 
         scale_y_continuous(labels = scales::label_number(accuracy = 1), limits = c(25, 45)) +
         coord_flip()
versc.male <- ggplot(broncho.male,aes(x=Pop,y=VerSC, fill = Pop, color = Pop))+
         theme_bw() +
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = VerSC), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab("VerSC") + xlab("") + ggtitle("Male") + 
         scale_color_manual(values = plot_color.m) + scale_fill_manual(values = plot_color.m) + 
         theme(axis.text=element_text(size=15, color = "black"), 
                    axis.title=element_text(size=15, face="bold"), 
                    legend.text=element_text(size=15, face = "bold"),
                    plot.title = element_text(size = 17, face = "bold")) + 
         scale_y_continuous(labels = scales::label_number(accuracy = 1), limits = c(25, 45)) +
         coord_flip()
versc.female <- ggplot(broncho.female,aes(x=Pop,y=VerSC, fill = Pop, color = Pop))+
         theme_bw() +
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = VerSC), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab("VerSC") + xlab("") + ggtitle("Female") + 
         scale_color_manual(values = plot_color) + scale_fill_manual(values = plot_color) + 
         theme(axis.text=element_text(size=15, color = "black"), 
                    axis.title=element_text(size=15, face="bold"), 
                    legend.text=element_text(size=15, face = "bold"),
                    plot.title = element_text(size = 17, face = "bold")) + 
         scale_y_continuous(labels = scales::label_number(accuracy = 1), limits = c(25, 45)) +
         coord_flip()
popv.pld <- ggplot(broncho.pld,aes(x=Pop,y=PoPV, fill = Pop, color = Pop))+
         theme_bw() +
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = PoPV), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab("PopV") + xlab("") + ggtitle("") + 
         scale_color_manual(values = plot_color) + scale_fill_manual(values = plot_color) + 
         theme(axis.text=element_text(size=15, color = "black"), 
                    axis.title=element_text(size=15, face="bold"), 
                    legend.text=element_text(size=15, face = "bold"),
                    plot.title = element_text(size = 17, face = "bold")) + 
         scale_y_continuous(labels = scales::label_number(accuracy = 1), limits = c(2, 10)) +
         coord_flip()
popv.male <- ggplot(broncho.male,aes(x=Pop,y=PoPV, fill = Pop, color = Pop))+
         theme_bw() +
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = PoPV), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab("PopV") + xlab("") + ggtitle("") + 
         scale_color_manual(values = plot_color.m) + scale_fill_manual(values = plot_color.m) + 
         theme(axis.text=element_text(size=15, color = "black"), 
                    axis.title=element_text(size=15, face="bold"), 
                    legend.text=element_text(size=15, face = "bold"),
                    plot.title = element_text(size = 17, face = "bold")) + 
         scale_y_continuous(labels = scales::label_number(accuracy = 1), limits = c(2, 10)) +
         coord_flip()
popv.female <- ggplot(broncho.female,aes(x=Pop,y=PoPV, fill = Pop, color = Pop))+
         theme_bw() +
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = PoPV), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab("PopV") + xlab("") + ggtitle("") + 
         scale_color_manual(values = plot_color) + scale_fill_manual(values = plot_color) + 
         theme(axis.text=element_text(size=15, color = "black"), 
                    axis.title=element_text(size=15, face="bold"), 
                    legend.text=element_text(size=15, face = "bold"),
                    plot.title = element_text(size = 17, face = "bold")) + 
         scale_y_continuous(labels = scales::label_number(accuracy = 1), limits = c(2, 10)) +
         coord_flip()
slf3.pld <- ggplot(broncho.pld,aes(x=Pop,y=SlF3, fill = Pop, color = Pop))+
         theme_bw() +
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = SlF3), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab("SlF3") + xlab("") + ggtitle("") + 
         scale_color_manual(values = plot_color) + scale_fill_manual(values = plot_color) + 
         theme(axis.text=element_text(size=15, color = "black"), 
                    axis.title=element_text(size=15, face="bold"), 
                    legend.text=element_text(size=15, face = "bold"),
                    plot.title = element_text(size = 17, face = "bold")) + 
         scale_y_continuous(labels = scales::label_number(accuracy = 1), limits = c(20, 35)) +
         coord_flip()
slf3.male <- ggplot(broncho.male,aes(x=Pop,y=SlF3, fill = Pop, color = Pop))+
         theme_bw() +
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = SlF3), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab("SlF3") + xlab("") + ggtitle("") + 
         scale_color_manual(values = plot_color.m) + scale_fill_manual(values = plot_color.m) + 
         theme(axis.text=element_text(size=15, color = "black"), 
                    axis.title=element_text(size=15, face="bold"), 
                    legend.text=element_text(size=15, face = "bold"),
                    plot.title = element_text(size = 17, face = "bold")) + 
         scale_y_continuous(labels = scales::label_number(accuracy = 1), limits = c(20, 35)) +
         coord_flip()
slf3.female <- ggplot(broncho.female,aes(x=Pop,y=SlF3, fill = Pop, color = Pop))+
         theme_bw() +
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = SlF3), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab("SlF3") + xlab("") + ggtitle("") + 
         scale_color_manual(values = plot_color) + scale_fill_manual(values = plot_color) + 
         theme(axis.text=element_text(size=15, color = "black"), 
                    axis.title=element_text(size=15, face="bold"), 
                    legend.text=element_text(size=15, face = "bold"),
                    plot.title = element_text(size = 17, face = "bold")) + 
         scale_y_continuous(labels = scales::label_number(accuracy = 1), limits = c(20, 35)) +
         coord_flip()

bx.mer <- ggarrange(versc.pld,versc.male,versc.female,popv.pld,popv.male,popv.female, slf3.pld, slf3.male, slf3.female,
            nrow = 3, ncol=3, legend = "none", common.legend = TRUE) 

#export result
png("mer_boxplot.png", res=300, width = 15, height = 12, unit="in")
bx.mer
dev.off()

pdf("mer_boxplot.pdf", width = 15, height = 12)
bx.mer
dev.off()

knitr::include_graphics("mer_boxplot.png") 
```
