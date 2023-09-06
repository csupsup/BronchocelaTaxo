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
**3. Perform ANOVA and Kurskal-Wallis on mensural characters**
```{r mens.anova.kw}
#get characters with normal distribution based on shapiro-wilk test
broncho.norm <- subset(broncho.mens.trans, select = c(1, 3:7, 9))

#perform anova on a single character
aov.ED <- broncho.mens.trans %>% anova_test(ED ~ Pop)

#loop through all characters
results <- NULL
results <- as.data.frame(results)
for(col in 2:ncol(broncho.norm))
  {results[col-1,1] <- names(broncho.norm)[col]
    aov.res <- broncho.norm %>% anova_test(formula(paste0(names(broncho.norm)[col],"~","Pop")))
    results[col-1,2] <- aov.res$p[1]
    results[col-1,3] <- aov.res$F[1]}
colnames(results)[2] <- "pvalue"
colnames(results)[3] <- "F"

#write results
write.csv(results, "aov_results.csv")

#perform tukey hsd
ED.tk <- aov(ED ~ Pop, data = broncho.mens.trans) %>% tukey_hsd()
write.csv(ED.tk, "ED.tk.csv")
OD.tk <- aov(OD ~ Pop, data = broncho.mens.trans) %>% tukey_hsd()
write.csv(OD.tk, "OD.tk.csv")
HW.tk <- aov(HW ~ Pop, data = broncho.mens.trans) %>% tukey_hsd()
write.csv(HW.tk, "HW.tk.csv")
HD.tk <- aov(HD ~ Pop, data = broncho.mens.trans) %>% tukey_hsd()
write.csv(HD.tk, "HD.tk.csv")
HL.tk <- aov(HL ~ Pop, data = broncho.mens.trans) %>% tukey_hsd()
write.csv(HL.tk, "HL.tk.csv")
FbL.tk <- aov(FbL ~ Pop, data = broncho.mens.trans) %>% tukey_hsd()
write.csv(FbL.tk, "FbL.tk.csv")

#perform kruskal-wallis
broncho.non.norm <- subset(broncho.mens.trans, select = c(1,2,8,10:15))

#single character
kw <- broncho.non.norm %>% kruskal_test(SVL ~ Pop)

#loop through all characters
results <- NULL
results <- as.data.frame(results)
for(col in 2:ncol(broncho.non.norm))
  { results[col-1,1] <- names(broncho.non.norm)[col]
    kw.res <- broncho.non.norm %>% kruskal_test(formula(paste0(names(broncho.non.norm)[col],"~","Pop")))
    results[col-1,2] <- kw.res$p[1]
    results[col-1,3] <- kw.res$statistic[1]}
colnames(results)[2] <- "pvalue"
colnames(results)[3] <- "F"

TL.kw <- broncho.tail.trans %>% kruskal_test(TL ~ Pop)
write.csv(TL.kw, 'TL.kw.csv')

#write results
write.csv(results, "kw_results.csv")

#perform Dunn test
SVL.dn <- broncho.non.norm %>% dunn_test(SVL ~ Pop, p.adjust.method = "bonferroni")
write.csv(SVL.dn, "SVL.dn.csv")
HdL.dn <- broncho.non.norm %>% dunn_test(HdL ~ Pop, p.adjust.method = "bonferroni")
write.csv(HdL.dn, "HdL.dn.csv")
FtL.dn <- broncho.non.norm %>% dunn_test(FtL ~ Pop, p.adjust.method = "bonferroni") 
write.csv(FtL.dn, "FtL.dn.csv")
HbL.dn <- broncho.non.norm %>% dunn_test(HbL ~ Pop, p.adjust.method = "bonferroni")
write.csv(HbL.dn, "HbL.dn.csv")
FnsH.dn <- broncho.non.norm %>% dunn_test(FnsH ~ Pop, p.adjust.method = "bonferroni")
write.csv(FnsH.dn, "FnsH.dn.csv")
SnsH.dn <- broncho.non.norm %>% dunn_test(SnsH ~ Pop, p.adjust.method = "bonferroni")
write.csv(SnsH.dn, "SnsH.dn.csv")
FnsW.dn <- broncho.non.norm %>% dunn_test(FnsW ~ Pop, p.adjust.method = "bonferroni")
write.csv(FnsW.dn, "FnsW.dn.csv")
SnsW.dn <- broncho.non.norm %>% dunn_test(SnsW ~ Pop, p.adjust.method = "bonferroni") 
write.csv(SnsW.dn, "SnsW.dn.csv")
```
**4. Perform Kurskal-Wallis on meristic characters**
```{r mer.kw}
#loop through all characters
results <- NULL
results <- as.data.frame(results)
for(col in 2:ncol(broncho.mer))
  { results[col-1,1] <- names(broncho.mer)[col]
    kw.res <- broncho.mer %>% kruskal_test(formula(paste0(names(broncho.mer)[col],"~","Pop")))
    results[col-1,2] <- kw.res$p[1]
    results[col-1,3] <- kw.res$statistic[1]}
colnames(results)[2] <- "pvalue"
colnames(results)[3] <- "F"

#write results
write.csv(results, "kw_results.mer.csv")

#perform Dunn test
UL.dn <- broncho.mer %>% dunn_test(UL ~ Pop, p.adjust.method = "bonferroni")
write.csv(UL.dn, "UL.dn.csv")
LL.dn <- broncho.mer %>% dunn_test(LL ~ Pop, p.adjust.method = "bonferroni")
write.csv(LL.dn, "LL.dn.csv")
CL.dn <- broncho.mer %>% dunn_test(CL ~ Pop, p.adjust.method = "bonferroni")
write.csv(CL.dn, "CL.dn.csv")
TdL.dn <- broncho.mer %>% dunn_test(TdL ~ Pop, p.adjust.method = "bonferroni")
write.csv(TdL.dn, "TdL.dn.csv")
FhL.dn <- broncho.mer %>% dunn_test(FhL ~ Pop, p.adjust.method = "bonferroni")
write.csv(FhL.dn, "FhL.dn.csv")
NS.dn <- broncho.mer %>% dunn_test(NS ~ Pop, p.adjust.method = "bonferroni")
write.csv(NS.dn, "NS.dn.csv")
MidSR.dn <- broncho.mer %>% dunn_test(MidSR ~ Pop, p.adjust.method = "bonferroni")
write.csv(MidSR.dn, "MidSR.dn.csv")
VerSC.dn <- broncho.mer %>% dunn_test(VerSC ~ Pop, p.adjust.method = "bonferroni")
write.csv(VerSC.dn, "VerSC.dn.csv")
DorSPU.dn <- broncho.mer %>% dunn_test(DorSPU ~ Pop, p.adjust.method = "bonferroni")
write.csv(DorSPU.dn, "DorSPU.dn.csv")
DorSPB.dn <- broncho.mer %>% dunn_test(DorSPB ~ Pop, p.adjust.method = "bonferroni")
write.csv(DorSPB.dn, "DorSPB.dn.csv")
TL.dn <- broncho.tail %>% dunn_test(TL ~ Pop, p.adjust.method = "bonferroni")
write.csv(TL.dn, "TL.dn.csv")
```
**5. Visualize post hoc test results**
```{r phoc.vis, warning=FALSE}
phoc.pval <- read.csv("posthoc_pvalue.csv", row.names = 1)#data prepared manually by pooling all p-values

library(ComplexHeatmap)

phoc.map <- ComplexHeatmap::pheatmap(phoc.pval, cluster_rows = F, cluster_cols = F,
                         display_numbers = T, #cutree_rows = 4, 
                         color = c("cornflowerblue","cornsilk"), breaks = c(0.000,0.04,1.0), 
                         border_color = "black", number_color = "black",
                         fontsize_number = 12, fontsize = 12, fontface = "bold", column_names_side = c("top"),
                         angle_col = c("0"), heatmap_legend_param = list(title = "p-value", labels_gp = gpar(fontsize = 15), 
                         title_gp = gpar(fontsize = 16, fontface = "bold")), cellwidth = 50, cellheight = 17)

phoc.map 

tiff("phoc.pheatmap.tif", res=300, width = 20.5, height = 9, unit="in")
phoc.map
dev.off()

```
**6. Character classification tree**
```{r, CTA}
library(rpart)
library(rpart.plot)
library(caret)
library(caTools)

#broncho.sl <- subset(broncho.data.mod, select = c(1,3:28))

#A. perfrom using transformed data
#read data
broncho.cta.data <- broncho.data.mod %>% mutate(PhyGrp = ifelse(Pop == "Borneo", "Borneo", 
  ifelse(Pop == "Luzon", "Luzon Clade", ifelse(Pop == "Visayas", "Luzon Clade", 
    ifelse(Pop == "Romblon", "Luzon Clade", ifelse (Pop == "Lubang", "Luzon Clade", 
      ifelse(Pop == "Palawan", "Palawan Clade", ifelse(Pop == "Mindoro", "Mindoro Clade", 
        ifelse(Pop == "Mindanao", "Mindanao+EasID Clade", 
          ifelse(Pop == "EastID", "Mindanao+EasID Clade", "No data"))))))))))

#perform CTA
broncho.cta <- rpart(as.factor(PhyGrp)~., data = broncho.cta.data[,2:26], method = "class")

#plot result
#rpart.plot(broncho.cta, main = "Classification Tree")

#check other info
#printcp(broncho.cta)

#select best split
best <- broncho.cta$cptable[which.min(broncho.cta$cptable[,"xerror"]),"CP"]
pruned_tree <- prune(broncho.cta, cp=best)
rpart.plot(pruned_tree, main = "")

#evaluate results using confusion matrix
preds <- predict(broncho.cta, broncho.cta.data[,2:26], type = "class")
confusionMatrix(as.factor(broncho.cta.data$PhyGrp), preds)

#export
#extra=106 shows probability of observation and percent of observation
tiff("class_tree.tif", res=300, width = 10, height = 5, unit="in")
rpart.plot(pruned_tree, main = "", tweak = 1.1, extra = 106, legend.x = -0.05, 
  legend.y = 1.2, split.cex=1.5, legend.cex = 1.3, boxes.include.gap = TRUE, 
  yesno = 2, split.border.col=1) 
dev.off()

#B. perform using raw data
broncho.raw.data <- subset(broncho.data, select = c(1, 3:26))

broncho.cta.raw.data <- broncho.raw.data %>% mutate(PhyGrp = ifelse(Pop == "Borneo", "Borneo", 
  ifelse(Pop == "Luzon", "Luzon Clade", ifelse(Pop == "Visayas", "Luzon Clade", 
    ifelse(Pop == "Romblon", "Luzon Clade", ifelse (Pop == "Lubang", "Luzon Clade", 
      ifelse(Pop == "Palawan", "Palawan Clade", ifelse(Pop == "Mindoro", "Mindoro Clade", 
        ifelse(Pop == "Mindanao", "Mindanao+EasID Clade", 
          ifelse(Pop == "EastID", "Mindanao+EasID Clade", "No data"))))))))))

broncho.cta <- rpart(as.factor(PhyGrp)~., data = broncho.cta.raw.data[,2:26], method = "class")
rpart.plot(broncho.cta, main = "Classification Tree (raw data)")
```
**7. Plot selected diagnostic characters**
```{r, plot.diag.char}
library(ggstatsplot)

broncho.data3 <- broncho.data.mod

broncho.data3$Pop <- factor(broncho.data3$Pop, levels = c("Luzon", "Visayas", "Romblon", "Lubang", 
  "Mindoro", "Palawan", "Mindanao", "EastID", "Borneo"))

plot_color <- c("#7f404a", "#cc1489", "#cc9666", "#408073", "#cc3d3d", 
  "#8c994d", "#7b3dcc", "#1262b3", "#000000")

ggbetweenstats(data = broncho.data3,
  x = Pop,
  y = SVL,
  type = "nonparametric", #Kruskal-Wallis
  var.equal = FALSE,
  plot.type = "box",
  pairwise.comparisons = FALSE,
  pairwise.display = "significant",
  centrality.plotting = FALSE,
  bf.message = FALSE,
  xlab = "") + ggplot2::scale_color_manual(values = plot_color) + ggplot2::theme_classic()

```

**8. Visualize data**
```{r visual, fig.width=15, fig.height=10}

#set color
plot_color <- c("#7f404a", "#cc1489", "#cc9666", "#408073", "#cc3d3d", 
  "#8c994d", "#7b3dcc", "#1262b3", "#000000")

#make a copy of the data and reorder group
broncho.data2 <- broncho.data

broncho.data2$Pop <- factor(broncho.data2$Pop, levels = c("Luzon", "Visayas", "Romblon", "Lubang", 
  "Mindoro", "Palawan", "Mindanao", "EastID", "Borneo"))

#create raincloud plots for diagnostic characters with heavy loading
p1 <- ggplot(broncho.data2,aes(x=Pop,y=SVL, fill = Pop, color = Pop))+
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = SVL), outlier.shape = NA, alpha = 0.3,
          width = 0.3, colour = "BLACK") + ylab('SVL') + xlab('') + theme_cowplot()  + 
         scale_color_manual(values = plot_color)+ scale_fill_manual(values = plot_color) + 
         theme(legend.position = "none", axis.text.x = element_text(size = 12.5), 
          axis.title = element_text(size = 14,  face ="bold")) + 
         coord_flip()

p2 <- ggplot(broncho.data2,aes(x=Pop,y=ED, fill = Pop, color = Pop))+
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = ED), outlier.shape = NA, alpha = 0.3, 
          width = 0.3, colour = "BLACK") + ylab('ED') + xlab('') + theme_cowplot() + 
         scale_color_manual(values = plot_color) + scale_fill_manual(values = plot_color) + 
         theme(legend.position = "none", axis.text.x = element_text(size = 12.5), 
          axis.title = element_text(size = 14, face ="bold")) + 
         coord_flip()

p3 <- ggplot(broncho.data2,aes(x=Pop,y=OD, fill = Pop, color = Pop))+
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = OD), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab('OD') + xlab('') + theme_cowplot() + scale_color_manual(values = plot_color) + 
         scale_fill_manual(values = plot_color) + 
         theme(legend.position = "none", axis.text.x = element_text(size = 12.5), 
          axis.title = element_text(size = 14,  face ="bold")) + 
         coord_flip()

p4 <- ggplot(broncho.data2,aes(x=Pop,y=HdL, fill = Pop, color = Pop))+
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = HdL), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab('HdL') + xlab('') + theme_cowplot() + scale_color_manual(values = plot_color) + 
         scale_fill_manual(values = plot_color) + 
         theme(legend.position = "none", axis.text.x = element_text(size = 12.5), 
          axis.title = element_text(size = 14,  face ="bold")) + 
         coord_flip()

p5 <- ggplot(broncho.data2,aes(x=Pop,y=HbL, fill = Pop, color = Pop))+
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = HbL), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab('HbL') + xlab('') + theme_cowplot() + 
         scale_color_manual(values = plot_color) + scale_fill_manual(values = plot_color) + 
         theme(legend.position = "none", axis.text.x = element_text(size = 12.5), 
          axis.title = element_text(size = 14, face ="bold")) + 
         coord_flip()

p6 <- ggplot(broncho.data2,aes(x=Pop,y=FnsH, fill = Pop, color = Pop))+
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = FnsH), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab('FnsH') + xlab('') + theme_cowplot() + scale_color_manual(values = plot_color) + 
         scale_fill_manual(values = plot_color) + 
         theme(legend.position = "none", axis.text.x = element_text(size = 12.5), 
          axis.title = element_text(size = 14, face ="bold"))+ 
         coord_flip()

p7 <- ggplot(broncho.data2,aes(x=Pop,y=FnsW, fill = Pop, color = Pop))+
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = FnsW), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab('FnsW') + xlab('') + theme_cowplot() + 
         scale_color_manual(values = plot_color) + scale_fill_manual(values = plot_color) + 
         theme(legend.position = "none", axis.text.x = element_text(size = 12.5), 
          axis.title = element_text(size = 14, face ="bold")) + 
         coord_flip()

p8 <- ggplot(broncho.data2,aes(x=Pop,y=DorSPB, fill = Pop, color = Pop))+
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = DorSPB), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab('DorSPB') + xlab('') + theme_cowplot() + 
         scale_color_manual(values = plot_color) + scale_fill_manual(values = plot_color) + 
         theme(legend.position = "none", axis.text.x = element_text(size = 12.5), 
          axis.title = element_text(size = 14, face ="bold")) + 
         coord_flip()

p9 <- ggplot(broncho.data2,aes(x=Pop,y=DorSPU, fill = Pop, color = Pop))+
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = DorSPU), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab('DorSPU') + xlab('') + theme_cowplot() + scale_color_manual(values = plot_color) + 
         scale_fill_manual(values = plot_color) + 
         theme(legend.position = "none", axis.text.x = element_text(size = 12.5), 
          axis.title = element_text(size = 14, face ="bold")) + 
         coord_flip()

p10 <- ggplot(broncho.data2,aes(x=Pop,y=MidSR, fill = Pop, color = Pop))+
         geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2, trim = TRUE)+
         geom_point(position = position_jitter(width = .10, height = 0), size = .5)+
         geom_boxplot(aes(x = as.numeric(as.factor(Pop))+0.25, y = MidSR), outlier.shape = NA, 
          alpha = 0.3, width = 0.3, colour = "BLACK") + 
         ylab('MidSR') + xlab('') + theme_cowplot() + scale_color_manual(values = plot_color) + 
         scale_fill_manual(values = plot_color) + 
         theme(legend.position = "none", axis.text.x = element_text(size = 12.5), 
          axis.title = element_text(size = 14, face ="bold")) + 
         coord_flip()

grid.arrange (p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,nrow = 2, ncol=5) 

#export plots
tiff("rcplot_pca1.tif", res=300, width = 19, height = 9, unit="in")
ggarrange (p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,nrow = 2, ncol=5, align = "hv") 
dev.off()
```
