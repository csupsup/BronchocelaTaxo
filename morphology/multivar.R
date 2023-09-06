---
title: "Bronchocela Morphology - Multivariate Analysis"
author: "Christian Supsup"
date: "24 July 2023"
output:
  html_document:
    df_print: paged
  pdf_document: default
indent: no
---
```{r setup, include=FALSE, message=TRUE}
knitr::opts_chunk$set(echo = TRUE)
library(knitr)
```
**1. Read data**
```{r read data, results='hide'}
library(mda)
library(ggplot2)
library(GroupStruct)
library(dplyr)
library(tidyverse)
library(ggbiplot)
library(gridExtra)
library(ggpubr)
library(corrplot)
library(PCAtest)
library(nzilbb.vowels)

broncho.data <- read.csv("Bronchocela_23July2023_multivar.csv", sep = ",", header = TRUE)
broncho.mens <- subset(broncho.data, select = c(1, 3:16)) #get mensural data
broncho.mens.trans <- allom(broncho.mens, type = "population1") #perform allometry correction
broncho.mer <- subset(broncho.data, select = c(1, 17:26)) #get meristic data
broncho.data.mod <- cbind(broncho.mens.trans, broncho.mer[2:11])#merge transformed mensural and meristic data
head(broncho.data.mod)
```

**2. Perform Flexible Discriminant Analysis (FDA)**
```{r perform FDA}
#perform FDA
broncho.data.trans <- broncho.data.mod
broncho.data.trans.mens <- subset(broncho.data.mod, select = c(1, 2:15))
broncho.data.trans.mer <- subset(broncho.data.mod, select = c(1, 16:25))

broncho.fda <- fda(Pop ~., data = broncho.data.trans) 
broncho.fda.mens <- fda(Pop ~., data = broncho.data.trans.mens)
broncho.fda.mer <- fda(Pop ~., data = broncho.data.trans.mer)

#generate fda confusion matrix for all data
fda.con.mat.all <- table(broncho.data.trans$Pop, predict(broncho.fda))
fda.con.mat.all
write.csv(fda.con.mat.all, "fda.con.mat.all.csv")

#check accuracy of classification
tbl <- table(broncho.data.trans$Pop, predict(broncho.fda))
sum(diag(tbl))/sum(tbl)

#generate fda confusion matrix for mensural data
fda.con.mat.mens <- table(broncho.data.trans.mens$Pop, predict(broncho.fda.mens))
fda.con.mat.mens
write.csv(fda.con.mat.mens, "fda.con.mat.men.csv")

#check accuracy of classification
tbl <- table(broncho.data.trans$Pop, predict(broncho.fda.mens))
sum(diag(tbl))/sum(tbl)

#generate fda confusion matrix for meristic data
fda.con.mat.mer <- table(broncho.data.trans.mer$Pop, predict(broncho.fda.mer))
fda.con.mat.mer
write.csv(fda.con.mat.mer, "fda.con.mat.mer.csv")

#check accuracy of classification
tbl <- table(broncho.data.trans$Pop, predict(broncho.fda.mer))
sum(diag(tbl))/sum(tbl)

#convex hull
fda.axes <- predict(broncho.fda, type = "var")
fda.axes2 <- as.data.frame(fda.axes)
phygrp.data <- broncho.data$PhyGrp
broncho.fda.hull <- cbind(phygrp.data, fda.axes2)

hull_data <- 
    broncho.fda.hull %>%
    drop_na() %>%
    group_by(broncho.data$PhyGrp) %>% 
    slice(chull(V1, V2))

#get axes values for customize plotting
#all data
fda.axes <- predict(broncho.fda, type = "var")
fda.axes2 <- as.data.frame(fda.axes)
paic.data <- broncho.data$Pop
broncho.fda.data <- cbind(paic.data, fda.axes2)

#mensural data
mens.fda.axes <- predict(broncho.fda.mens, type = "var")
mens.fda.axes2 <- as.data.frame(mens.fda.axes)
paic.data <- broncho.data$Pop
broncho.fda.data.mens <- cbind(paic.data, mens.fda.axes2)

#meristic data
mer.fda.axes <- predict(broncho.fda.mer, type = "var")
mer.fda.axes2 <- as.data.frame(mer.fda.axes)
paic.data <- broncho.data$Pop
broncho.fda.data.mer <- cbind(paic.data, mer.fda.axes2)
```
**3. View FDA results**
```{r view FDA resutls, fig.width=15}
#assign colors to each Pop
plot_color <- c("#000000", "#1262b3", "#408073", "#7f404a", "#7b3dcc", "#cc3d3d", 
  "#8c994d","#cc9666", "#cc1489")

fda0 <- ggplot(broncho.fda.data, aes(V1, V2))+ aes(fill = factor(phygrp.data)) + 
          geom_polygon(data = hull_data, alpha = 0.3) + 
          scale_fill_discrete(labels=c("Borneo", "Luzon Clade", "Mindanao + EastID Clade", 
            "Mindoro Clade", "Palawan Clade"))

fda1 <- fda0 + geom_point(aes(color=paic.data, shape = paic.data), size = 3, alpha =0.9, fill = "#cc3d3d") + 
          scale_shape_manual(values = c(8,17,11,7,19,25,15,18,1)) + 
          scale_color_manual(values = plot_color) + theme(plot.subtitle = element_text(size=15), legend.key = element_blank(), 
          legend.title = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "grey"), legend.position = c(0.85, 0.9), 
          legend.key.width = unit(1.5, "cm"), legend.key.height = unit(0.5, "cm"), 
          legend.text = element_text(size=14), plot.title = element_text(color="black", size=17, face="bold"), 
          axis.text=element_text(size=15), axis.title=element_text(size=15, face="bold")) + 
          xlab("FDA1, accuracy = 95.07%") + ylab("FDA2") + ggtitle("") + 
          guides(colour = guide_legend(override.aes = list(size=3))) + 
          scale_x_continuous(label = ~ scales::comma(.x, accuracy = 1)) + 
          labs(subtitle = "Mensural & meristic") #+ xlim(-6, 6) + ylim(-6, 6)

fda2 <- ggplot(broncho.fda.data.mens, aes(V1, V2))+ geom_point(aes(color=paic.data, shape = paic.data),size= 3, 
          alpha =0.9, fill = "#cc3d3d") + 
          scale_shape_manual(values = c(8,17,11,7,19,25,15,18,1)) + scale_color_manual(values = plot_color) + 
          theme(plot.subtitle = element_text(size=15),legend.key = element_blank(), 
          legend.title = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "grey"), legend.position = "none", 
          plot.title = element_text(color="black", size=17, face="bold"),
          axis.text=element_text(size=15), axis.title=element_text(size=15, face="bold")) + 
          xlab("FDA1, accuracy = 83.80%") + ylab("FDA2") + ggtitle("") + 
          scale_x_continuous(label = ~ scales::comma(.x, accuracy = 1)) + 
          labs(subtitle = "Mensural") #+ xlim(-6, 5) + ylim(-6, 6)


fda3 <- ggplot(broncho.fda.data.mer, aes(V1, V2))+ 
          geom_point(aes(color=paic.data, shape = paic.data), size =3, alpha =0.9, fill = "#cc3d3d") + 
          scale_shape_manual(values = c(8,17,11,7,19,25,15,18,1)) + scale_color_manual(values = plot_color) + 
          theme(plot.subtitle = element_text(size=15),legend.key = element_blank(), legend.title = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "grey"), 
          plot.title = element_text(color="black", size=17, face="bold"), legend.position = "none", 
          axis.text=element_text(size=15), axis.title=element_text(size=15, face="bold")) + 
          xlab("FDA1, accuracy = 82.39%") + ylab("FDA2") + ggtitle("") + 
          scale_y_continuous(label = ~ scales::comma(.x, accuracy = 1)) + 
          labs(subtitle = "Meristic") #+ xlim(-6, 5) + ylim(-6, 6)

#grid.arrange(fda1, fda2, fda3,nrow = 1)
fda.plot <- ggarrange(fda1, fda2, fda3, ncol =3, nrow=1, common.legend=TRUE, legend="right")
fda.plot
#export plot
tiff("broncho_fda_combined.tif", res=300, width = 16, height = 5, unit="in")
ggarrange(fda1, fda2, fda3, ncol =3, nrow=1, common.legend=TRUE, legend="right")
dev.off()
```

**5. Perform Principal Component Analysis (PCA)**

```{r perform pca}
#perform PCA
broncho.pca.all <- prcomp(broncho.data.trans[2:25], center = TRUE, scale. = TRUE)
broncho.pca.mens <- prcomp(broncho.data.trans.mens[2:15], center = TRUE, scale. = TRUE)
broncho.pca.mer <- prcomp(broncho.data.trans.mer[2:11], center = TRUE, scale. = TRUE)

#calculate variance explained for all data
eigs.all <- broncho.pca.all$sdev^2
ax.1.all <- eigs.all[1]/sum(eigs.all)*100
ax.1.all.per <- sprintf(ax.1.all, fmt = '%#.2f')
ax.2.all <- eigs.all[2]/sum(eigs.all)*100
ax.2.all.per <- sprintf(ax.2.all, fmt = '%#.2f')

ax.1.all.per
ax.2.all.per

#calculate variance explained for mensural
eigs.mens <- broncho.pca.mens$sdev^2
ax.1.mens <- eigs.mens[1]/sum(eigs.mens)*100
ax.1.mens.per <- sprintf(ax.1.mens, fmt = '%#.2f')
ax.2.mens <- eigs.mens[2]/sum(eigs.mens)*100
ax.2.mens.per <- sprintf(ax.2.mens, fmt = '%#.2f')

ax.1.mens.per
ax.2.mens.per

#calculate variance explained for meristic
eigs.mer <- broncho.pca.mer$sdev^2
ax.1.mer <- eigs.mer[1]/sum(eigs.mer)*100
ax.1.mer.per <- sprintf(ax.1.mer, fmt = '%#.2f')
ax.2.mer <- eigs.mer[2]/sum(eigs.mer)*100
ax.2.mer.per <- sprintf(ax.2.mer, fmt = '%#.2f')

ax.1.mer.per
ax.2.mer.per

#or use summary function
summary(broncho.pca.all)
summary(broncho.pca.mens)
summary(broncho.pca.mer)

#pca convex hull
#extract point loadings for all characters
pca_points <-  as_tibble(broncho.pca.all$x) %>% bind_cols(broncho.data$PhyGrp) %>% bind_cols(broncho.data$Pop)
colnames(pca_points)[25] <- "phygrp"
colnames(pca_points)[26] <- "Pop"

pca_hull <- 
  pca_points %>% 
  group_by(phygrp) %>% 
  slice(chull(PC1, PC2))

#extract variable loadings for all characters
pca_load <- as_tibble(broncho.pca.all$rotation, rownames = 'variable')
pca_load_mens <- as_tibble(broncho.pca.mens$rotation, rownames = 'variable')
pca_load_mer <- as_tibble(broncho.pca.mer$rotation, rownames = 'variable')

#extract point loadings for mensural characters
pca_points_mens <-  as_tibble(broncho.pca.mens$x) %>% bind_cols(broncho.data$PhyGrp) %>% bind_cols(broncho.data$Pop)
colnames(pca_points_mens)[15] <- "phygrp"
colnames(pca_points_mens)[16] <- "Pop"

#extract point loadings for meristic characters
pca_points_mer <-  as_tibble(broncho.pca.mer$x) %>% bind_cols(broncho.data$PhyGrp)  %>% bind_cols(broncho.data$Pop)
colnames(pca_points_mer)[11] <- "phygrp"
colnames(pca_points_mer)[12] <- "Pop"

```
**6. Assess PCA significance**
```{r pca sig}
#assess PC and variable signifance using PCAtest
broncho.pca.sig.all <- PCAtest(broncho.data.trans[2:25], 1000, 1000, 0.05, varcorr=FALSE, counter=FALSE, plot=TRUE)

broncho.pca.sig.mens <- PCAtest(broncho.data.trans.mens[2:15], 1000, 1000, 0.05, varcorr=FALSE, counter=FALSE, plot=TRUE)

broncho.pca.sig.mer <- PCAtest(broncho.data.trans.mer[2:11], 1000, 1000, 0.05, varcorr=FALSE, counter=FALSE, plot=TRUE)

#plot results using nzilbb_vowels

pca.sig.all <- pca_test(broncho.data.trans[2:25], n = 1000, scale = TRUE, variance_confint = .95, loadings_confint = 0.95)

p1 <- plot_variance_explained(pca.sig.all, pc_max = NA, percent = TRUE)
p2 <- plot_loadings(pca.sig.all, pc_no = 1, violin = FALSE, filter_boots = FALSE, quantile_threshold = FALSE)
p3 <- plot_loadings(pca.sig.all, pc_no = 2, violin = FALSE, filter_boots = FALSE, quantile_threshold = FALSE)

pca.sig.mens <- pca_test(broncho.data.trans.mens[2:15], n = 1000, scale = TRUE, variance_confint = .95, loadings_confint = 0.95)

p4 <- plot_variance_explained(pca.sig.mens, pc_max = NA, percent = TRUE)
p5 <- plot_loadings(pca.sig.mens, pc_no = 1, violin = FALSE, filter_boots = FALSE, quantile_threshold = FALSE)
p6 <- plot_loadings(pca.sig.mens, pc_no = 2, violin = FALSE, filter_boots = FALSE, quantile_threshold = FALSE)

pca.sig.mer <- pca_test(broncho.data.trans.mer[2:11], n = 1000, scale = TRUE, variance_confint = .95, loadings_confint = 0.95)

p7 <- plot_variance_explained(pca.sig.mer, pc_max = NA, percent = TRUE)
p8 <- plot_loadings(pca.sig.mer, pc_no = 1, violin = FALSE, filter_boots = FALSE, quantile_threshold = FALSE)
p9 <- plot_loadings(pca.sig.mer, pc_no = 2, violin = FALSE, filter_boots = FALSE, quantile_threshold = FALSE)

plot.sig <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9,ncol =3, nrow=3, common.legend=TRUE, legend="right")
plot.sig

tiff("pca_sig.tif", res=300, width = 11, height = 8, unit="in")
ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9,ncol =3, nrow=3, common.legend=TRUE, legend="right")
dev.off()
```

**7. View PCA results**

```{r view PCA results, fig.width=15}

p1 <- ggplot(pca_points, aes(x = PC1, y = PC2)) + 
        geom_polygon(data = pca_hull, aes(fill = phygrp), alpha = 0.3, show.legend = FALSE) + 
        geom_point(aes(color=Pop, shape = Pop), size = 3, fill = "#cc3d3d", show.legend = FALSE) + 
        scale_shape_manual(name = pca_points$Pop, values = c(8,17,11,7,19,25,15,18,1)) + 
        scale_color_manual(name = pca_points$Pop, values = plot_color) + labs(x = "PC1", y = "PC2") + 
        theme(legend.key=element_blank(), plot.subtitle = element_text(size=13),plot.margin = margin(0.2, 0.5, 0, 0.2, "cm"), 
          legend.title = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "grey"), 
          legend.position = c(0.85, 0.9), legend.key.width = unit(1, "cm"), legend.key.height = unit(0.3, "cm"), 
          legend.text = element_text(size=12), plot.title = element_text(color="black", size=17, face="bold"), 
          axis.text=element_text(size=15), axis.title=element_text(size=15, face="bold")) + 
        ggtitle("") + guides(colour = guide_legend(override.aes = list(size=3))) + 
        scale_y_continuous(label = ~ scales::comma(.x, accuracy = 1)) + 
        labs(subtitle = "") + xlim(-5.5, 5.5) + ylim(-8, 6)

p2 <- ggplot(pca_points_mens, aes(x = PC1, y = PC2)) + 
        geom_point(aes(color=broncho.data$Pop, shape = broncho.data$Pop), size = 3, fill = "#cc3d3d") + 
        scale_shape_manual(name = broncho.data$Pop, values = c(8,17,11,7,19,25,15,18,1)) + 
        scale_color_manual(name = broncho.data$Pop, values = plot_color) + labs(x = "PC1", y = "PC2") + 
        theme(plot.subtitle = element_text(size=13), plot.margin = margin(0.2, 0.5, 0, 0.2, "cm"), 
          legend.title = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "grey"),
          legend.position = "none", legend.key.width = unit(1, "cm"), legend.key.height = unit(0.1, "cm"), 
          legend.text = element_text(size=12), plot.title = element_text(color="black", size=17, face="bold"), 
          axis.text=element_text(size=15), axis.title=element_text(size=15, face="bold"))  + 
        ggtitle("") + guides(colour = guide_legend(override.aes = list(size=3))) + 
        scale_x_continuous(label = ~ scales::comma(.x, accuracy = 1)) + 
        labs(subtitle = "") + xlim(-5.5, 5.5) + ylim(-8, 6)

p3 <- ggplot(pca_points_mer, aes(x = PC1, y = PC2)) + 
        geom_point(aes(color=broncho.data$Pop, shape = broncho.data$Pop), size = 3, fill = "#cc3d3d") + 
        scale_shape_manual(name = broncho.data$Pop, values = c(8,17,11,7,19,25,15,18,1)) + 
        scale_color_manual(name = broncho.data$Pop, values = plot_color) + labs(x = "PC1", y = "PC2") + 
        theme(plot.subtitle = element_text(size=13), plot.margin = margin(0.2, 0.5, 0, 0.2, "cm"), 
          legend.title = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "grey"), 
          legend.position = "none", legend.key.width = unit(1, "cm"), 
          legend.key.height = unit(0.1, "cm"), legend.text = element_text(size=12), 
          plot.title = element_text(color="black", size=17, face="bold"), axis.text=element_text(size=15), 
          axis.title=element_text(size=15, face="bold"))  + 
        ggtitle("") + guides(colour = guide_legend(override.aes = list(size=3))) + 
        scale_y_continuous(label = ~ scales::comma(.x, accuracy = 1)) + 
        labs(subtitle = "") + xlim(-5.5, 5.5) + ylim(-8, 6)

a1 <- ggplot(pca_points, aes(x = PC1, y = PC2)) + 
        geom_polygon(data = pca_hull, aes(fill = phygrp), alpha = 0.3, show.legend = FALSE)+ 
        geom_segment(data = pca_load, aes(x = 0, y = 0, xend = PC1*8, yend = PC2*8), 
          arrow = arrow(length = unit(1/2, 'picas'))) + 
        annotate('text', x = (pca_load$PC1*8), y = (pca_load$PC2*8.4),label = pca_load$variable, size = 3.5) + 
        labs(x = "PC1, 23.02% variance explained", y = "PC2, 13.94% variance explained") + 
        theme(plot.subtitle = element_text(size=13),plot.margin = margin(0.2, 0.5, 0, 0.2, "cm"), 
        legend.title = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "grey"), legend.position = "none", 
        legend.key.width = unit(1, "cm"), legend.key.height = unit(0.1, "cm"), 
        legend.text = element_text(size=12), plot.title = element_text(color="black", size=17, face="bold"),
         axis.text=element_text(size=15), axis.title=element_text(size=15, face="bold"))  + 
         ggtitle("") + guides(colour = guide_legend(override.aes = list(size=3))) + 
         scale_y_continuous(label = ~ scales::comma(.x, accuracy = 1)) + 
         xlim(-5.5, 5.5) + ylim(-8, 6)  #+ labs(subtitle = "Mensural and meristic")

a2 <- ggplot(pca_points_mens, aes(x = PC1, y = PC2)) + 
        geom_segment(data = pca_load_mens, aes(x = 0, y = 0, xend = PC1*8, yend = PC2*8), 
        arrow = arrow(length = unit(1/2, 'picas'))) + 
        annotate('text', x = (pca_load_mens$PC1*8), y = (pca_load_mens$PC2*8.4),label = pca_load_mens$variable, size = 3.5) +  
        labs(x = "PC1, 32.07% variance explained", y = "PC2, 20.13% variance explained") + 
        theme(plot.subtitle = element_text(size=13), plot.margin = margin(0.2, 0.5, 0, 0.2, "cm"), 
        legend.title = element_blank(), panel.background = element_blank(), 
        axis.line = element_line(colour = "grey"),legend.position = "none", 
        legend.key.width = unit(1, "cm"), legend.key.height = unit(0.1, "cm"), 
        legend.text = element_text(size=12), plot.title = element_text(color="black", size=17, face="bold"), 
        axis.text=element_text(size=15), axis.title=element_text(size=15, face="bold"))  + 
        ggtitle("") + guides(colour = guide_legend(override.aes = list(size=3))) + 
        scale_x_continuous(label = ~ scales::comma(.x, accuracy = 1)) + 
        xlim(-5.5, 5.5) + ylim(-8, 6) #+ labs(subtitle = "Mensural")

a3 <- ggplot(pca_points_mer, aes(x = PC1, y = PC2)) + 
        geom_segment(data = pca_load_mer, aes(x = 0, y = 0, xend = PC1*8, yend = PC2*8), 
        arrow = arrow(length = unit(1/2, 'picas'))) + 
        annotate('text', x = (pca_load_mer$PC1*8), y = (pca_load_mer$PC2*8.4),label = pca_load_mer$variable, size = 3.5) + 
        labs(x = "PC1, 22.72% variance explained", y = "PC2, 15.03% variance explained") + 
        theme(plot.subtitle = element_text(size=13), plot.margin = margin(0.2, 0.5, 0, 0.2, "cm"), 
        legend.title = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "grey"), 
        legend.position = "none", legend.key.width = unit(1, "cm"), legend.key.height = unit(0.1, "cm"), 
        legend.text = element_text(size=12), plot.title = element_text(color="black", size=17, face="bold"), 
        axis.text=element_text(size=15), axis.title=element_text(size=15, face="bold"))  + 
        ggtitle("") + guides(colour = guide_legend(override.aes = list(size=3))) + 
        scale_y_continuous(label = ~ scales::comma(.x, accuracy = 1)) + 
        xlim(-5.5, 5.5) + ylim(-8, 6) #+ labs(subtitle = "Meristic")

#grid.arrange(p1, p2, p3, a1, a2, a3,nrow = 2)

ggarrange(p1, p2, p3, ncol =3, nrow=1, common.legend=TRUE, legend="right")
ggarrange(a1, a2, a3, ncol =3, nrow=1)

pca.plot <- ggarrange(p1, p2, p3, a1, a2, a3, ncol =3, nrow=2, common.legend=TRUE, legend="none")

#export plot
tiff("broncho_pca_combined.tif", res=300, width = 15, height = 10, unit="in")
ggarrange(p1, p2, p3, a1, a2, a3, ncol =3, nrow=2, common.legend=TRUE, legend="none")
dev.off()

tiff("broncho_pca_fda_combined.tif", res=300, width = 17, height = 15, unit="in")
ggarrange(fda1, fda2, fda3, p1, p2, p3, a1, a2, a3, ncol =3, nrow=3, common.legend=TRUE, legend="right", align = "hv")
dev.off()

#plot variable loadings for all characters - PC1
#get loadings
pc.load.sub <- pca_load %>% select(variable, PC1)
#plot to check loadings with high values
pc.load.sub <- pc.load.sub %>% arrange(PC1) %>% mutate(char = factor(variable, unique(variable))) #sor first
pc.load.sub <- pc.load.sub %>% mutate(PC1.r = sprintf('%#.2f', PC1))

ggplot(pc.load.sub, aes(x = PC1, y = char)) + geom_point()

#get high values for costum plot
pc1.high.neg <- pc.load.sub %>% filter(PC1 <= -0.33)
pc1.high.pos <- pc.load.sub %>% filter(PC1 >= 0.22)

v1 <- ggplot(pc.load.sub, aes(x = PC1, y = char, label = PC1.r), color=col) + 
        geom_point(alpha=0.8, size = 7, color = "#00BFC4") + 
        geom_point(data=pc1.high.neg, aes(x=PC1, y=char), color = "#FC4E07", size = 7) + 
        geom_point(data=pc1.high.pos, aes(x=PC1, y=char), color = "#FC4E07", size = 7) + 
        geom_segment(aes(x = 0, y = char, xend = PC1, yend = char), color = "grey") + 
        labs(x = "PC1, 23.02% variance explained", y = "", subtitle = "Mensural & Meristic") + 
        ggtitle("") + theme(panel.background = element_blank(), plot.subtitle = element_text(size=13),
          plot.margin = margin(0.2, 0.5, 0, 0.2, "cm"), legend.title = element_blank(), axis.line = element_line(colour = "black"), 
          legend.position = "none", legend.key.width = unit(1, "cm"), legend.key.height = unit(0.1, "cm"), 
          legend.text = element_text(size=12), plot.title = element_text(color="black", size=17, face="bold"), 
          axis.text=element_text(size=15), axis.title=element_text(size=15, face="bold")) + 
        geom_vline(lty = 2, xintercept = 0) + scale_x_continuous(breaks = seq(-0.45, 0.45, by = 0.2)) + 
        geom_text(color = 'black', size = 3) 

#plot variable loadings for all characters - PC2
#get loadings
pc.load.sub <- pca_load %>% select(variable, PC2)
#plot to check loadings with high values
pc.load.sub <- pc.load.sub %>% arrange(PC2) %>% mutate(char = factor(variable, unique(variable))) #sor first
pc.load.sub <- pc.load.sub %>% mutate(PC1.r = sprintf('%#.2f', PC2))

ggplot(pc.load.sub, aes(x = PC2, y = char)) + geom_point()

#get high values for costum plot
pc2.high.neg <- pc.load.sub %>% filter(PC2 <= -0.18)
pc2.high.pos <- pc.load.sub %>% filter(PC2 >= 0.39)

v2 <- ggplot(pc.load.sub, aes(x = PC2, y = char, label = PC1.r), color=col) + 
        geom_point(alpha=0.8, size = 7, color = "#00BFC4") + 
        geom_point(data=pc2.high.neg, aes(x=PC2, y=char), color = "#FC4E07", size = 7)+ 
        geom_point(data=pc2.high.pos, aes(x=PC2, y=char), color = "#FC4E07", size = 7) + 
        geom_segment(aes(x = 0, y = char, xend = PC2, yend = char), color = "grey") + 
        labs(x = "PC2, 13.94% variance explained", y = "", subtitle = "") + ggtitle("") + 
        theme(panel.background = element_blank(), plot.subtitle = element_text(size=13),
          plot.margin = margin(0.2, 0.5, 0, 0.2, "cm"), legend.title = element_blank(), 
          axis.line = element_line(colour = "black"), legend.position = "none", 
          legend.key.width = unit(1, "cm"), legend.key.height = unit(0.1, "cm"), 
          legend.text = element_text(size=12), plot.title = element_text(color="black", size=17, face="bold"), 
          axis.text=element_text(size=15), axis.title=element_text(size=15, face="bold")) + 
        geom_vline(lty = 2, xintercept = 0) + 
        scale_x_continuous(breaks = seq(-0.45, 0.45, by = 0.2)) + 
        geom_text(color = 'black', size = 3)

#plot variable loadings for mensural characters - PC1
#get loadings
pc.load.sub <- pca_load_mens %>% select(variable, PC1)
#plot to check loadings with high values
pc.load.sub <- pc.load.sub %>% arrange(PC1) %>% mutate(char = factor(variable, unique(variable))) #sor first
pc.load.sub <- pc.load.sub %>% mutate(PC1.r = sprintf('%#.2f', PC1))

ggplot(pc.load.sub, aes(x = PC1, y = char)) + geom_point()

#get high values for costum plot
pc1.high.neg <- pc.load.sub %>% filter(PC1 <= -0.38)
pc1.high.pos <- pc.load.sub %>% filter(PC1 >= 0.20)

v3 <- ggplot(pc.load.sub, aes(x = PC1, y = char, label = PC1.r), color=col) + 
        geom_point(alpha=0.8, size = 7, color = "#00BFC4") + 
        geom_point(data=pc1.high.neg, aes(x=PC1, y=char), color = "#FC4E07", size = 7)+ 
        geom_point(data=pc1.high.pos, aes(x=PC1, y=char), color = "#FC4E07", size = 7) + 
        geom_segment(aes(x = 0, y = char, xend = PC1, yend = char), color = "grey") + 
        labs(x = "PC1, 32.07% variance explained", y = "", subtitle = "Mensural") + ggtitle("") + 
        theme(panel.background = element_blank(), plot.subtitle = element_text(size=13),
          plot.margin = margin(0.2, 0.5, 0, 0.2, "cm"), legend.title = element_blank(), 
          axis.line = element_line(colour = "black"), legend.position = "none", 
          legend.key.width = unit(1, "cm"), legend.key.height = unit(0.1, "cm"), 
          legend.text = element_text(size=12), plot.title = element_text(color="black", size=17, face="bold"), 
          axis.text=element_text(size=15), axis.title=element_text(size=15, face="bold")) + 
        geom_vline(lty = 2, xintercept = 0) + scale_x_continuous(breaks = seq(-0.45, 0.45, by = 0.2)) + 
        geom_text(color = 'black', size = 3) 

#plot variable loadings for mensural characters - PC2
#get loadings
pc.load.sub <- pca_load_mens %>% select(variable, PC2)
#plot to check loadings with high values
pc.load.sub <- pc.load.sub %>% arrange(PC2) %>% mutate(char = factor(variable, unique(variable))) #sor first
pc.load.sub <- pc.load.sub %>% mutate(PC1.r = sprintf('%#.2f', PC2))

ggplot(pc.load.sub, aes(x = PC2, y = char)) + geom_point()

#get high values for costum plot
pc2.high.neg <- pc.load.sub %>% filter(PC2 <= 0.04)
pc2.high.pos <- pc.load.sub %>% filter(PC2 >= 0.42)

v4 <- ggplot(pc.load.sub, aes(x = PC2, y = char, label = PC1.r), color=col) + 
        geom_point(alpha=0.8, size = 7, color = "#00BFC4") + 
        geom_point(data=pc2.high.neg, aes(x=PC2, y=char), color = "#FC4E07", size = 7) + 
        geom_point(data=pc2.high.pos, aes(x=PC2, y=char), color = "#FC4E07", size = 7) + 
        geom_segment(aes(x = 0, y = char, xend = PC2, yend = char), color = "grey") + 
        labs(x = "PC2, 20.13% variance explained", y = "", subtitle = "") + 
        ggtitle("") + theme(panel.background = element_blank(), plot.subtitle = element_text(size=13),
          plot.margin = margin(0.2, 0.5, 0, 0.2, "cm"), legend.title = element_blank(), 
          axis.line = element_line(colour = "black"), legend.position = "none", 
          legend.key.width = unit(1, "cm"), legend.key.height = unit(0.1, "cm"), 
          legend.text = element_text(size=12), plot.title = element_text(color="black", size=17, face="bold"), 
          axis.text=element_text(size=15), axis.title=element_text(size=15, face="bold")) + 
        geom_vline(lty = 2, xintercept = 0) + 
        scale_x_continuous(breaks = seq(-0.45, 0.45, by = 0.2)) + 
        geom_text(color = 'black', size = 3)

#plot variable loadings for meristic characters - PC1
#get loadings
pc.load.sub <- pca_load_mer %>% select(variable, PC1)
#plot to check loadings with high values
pc.load.sub <- pc.load.sub %>% arrange(PC1) %>% mutate(char = factor(variable, unique(variable))) #sor first
pc.load.sub <- pc.load.sub %>% mutate(PC1.r = sprintf('%#.2f', PC1))

ggplot(pc.load.sub, aes(x = PC1, y = char)) + geom_point()

#get high values for costum plot
pc1.high.neg <- pc.load.sub %>% filter(PC1 <= -0.45)
pc1.high.pos <- pc.load.sub %>% filter(PC1 >= 0.26)

v5 <- ggplot(pc.load.sub, aes(x = PC1, y = char, label = PC1.r), color=col) + 
        geom_point(alpha=0.8, size = 7, color = "#00BFC4") + 
        geom_point(data=pc1.high.neg, aes(x=PC1, y=char), color = "#FC4E07", size = 7) + 
        geom_point(data=pc1.high.pos, aes(x=PC1, y=char), color = "#FC4E07", size = 7) + 
        geom_segment(aes(x = 0, y = char, xend = PC1, yend = char), color = "grey") + 
        labs(x = "PC1, 22.72% variance explained", y = "", subtitle = "Meristic") + ggtitle("") + 
        theme(panel.background = element_blank(), plot.subtitle = element_text(size=13),plot.margin = margin(0.2, 0.5, 0, 0.2, "cm"), 
          legend.title = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none", 
          legend.key.width = unit(1, "cm"), legend.key.height = unit(0.1, "cm"), 
          legend.text = element_text(size=12), plot.title = element_text(color="black", size=17, face="bold"), 
          axis.text=element_text(size=15), axis.title=element_text(size=15, face="bold")) + 
        geom_vline(lty = 2, xintercept = 0) + 
        scale_x_continuous(breaks = seq(-0.45, 0.45, by = 0.2)) + 
        geom_text(color = 'black', size = 3) 

#plot variable loadings for mertistic characters - PC2
#get loadings
pc.load.sub <- pca_load_mer %>% select(variable, PC2)
#plot to check loadings with high values
pc.load.sub <- pc.load.sub %>% arrange(PC2) %>% mutate(char = factor(variable, unique(variable))) #sor first
pc.load.sub <- pc.load.sub %>% mutate(PC1.r = sprintf('%#.2f', PC2))

ggplot(pc.load.sub, aes(x = PC2, y = char)) + geom_point()

#get high values for costum plot
pc2.high.neg <- pc.load.sub %>% filter(PC2 <= -0.31)
pc2.high.pos <- pc.load.sub %>% filter(PC2 >= 0.41)

v6 <- ggplot(pc.load.sub, aes(x = PC2, y = char, label = PC1.r), color=col) + 
        geom_point(alpha=0.8, size = 7, color = "#00BFC4") + 
        geom_point(data=pc2.high.neg, aes(x=PC2, y=char), color = "#FC4E07", size = 7) + 
        geom_point(data=pc2.high.pos, aes(x=PC2, y=char), color = "#FC4E07", size = 7) + 
        geom_segment(aes(x = 0, y = char, xend = PC2, yend = char), color = "grey") + 
        labs(x = "PC2, 15.03% variance explained", y = "", subtitle = "") + ggtitle("") + 
        theme(panel.background = element_blank(), plot.subtitle = element_text(size=13),
          plot.margin = margin(0.2, 0.5, 0, 0.2, "cm"), legend.title = element_blank(), 
          axis.line = element_line(colour = "black"), legend.position = "none", 
          legend.key.width = unit(1, "cm"), legend.key.height = unit(0.1, "cm"), 
          legend.text = element_text(size=12), plot.title = element_text(color="black", size=17, face="bold"), 
          axis.text=element_text(size=15), axis.title=element_text(size=15, face="bold")) + 
        geom_vline(lty = 2, xintercept = 0) + 
        scale_x_continuous(breaks = seq(-0.45, 0.45, by = 0.2)) + 
        geom_text(color = 'black', size = 3)

#export plot
tiff("broncho_pca_loadings_combined.tif", res=300, width = 10, height = 15, unit="in")
ggarrange(v1, v2, v3, v4, v5, v6, ncol =2, nrow=3, common.legend=TRUE, legend="none", align = "hv")
dev.off()

tiff("broncho_pca_loadings1.tif", res=300, width = 10, height = 5, unit="in")
ggarrange(v1, v2, ncol =2, nrow=1, common.legend=TRUE, legend="none")
dev.off()

tiff("broncho_pca_loadings2.tif", res=300, width = 10, height = 5, unit="in")
ggarrange(v3, v4, ncol =2, nrow=1, common.legend=TRUE, legend="none")
dev.off()

tiff("broncho_pca_loadings3.tif", res=300, width = 10, height = 5, unit="in")
ggarrange(v5, v6, ncol =2, nrow=1, common.legend=TRUE, legend="none")
dev.off()
```
