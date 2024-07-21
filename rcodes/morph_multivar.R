---
title: "Bronchocela Morphology - Multivariate Analysis"
author: "Christian Supsup"
date: "26 March 2024"
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
library(PCAtest) #devtools::install_github("arleyc/PCAtest")
library(nzilbb.vowels) #devtools::install_github("csupsup/nzilbb_vowels")
library(ggrepel)
library(rpart)
library(rpart.plot)
library(caret)
library(caTools)

#read data
broncho.data.raw <- read.csv("Bronchocela_21March2024_multivar.csv", sep = ",", header = TRUE)
#subset data
broncho.data <- subset(broncho.data.raw, select = c(6, 12, 14:37))
broncho.mens <- subset(broncho.data, select = c(1, 3:16)) #get mensural data
broncho.mens.trans <- allom(broncho.mens, type = "population1") #perform allometry correction
names(broncho.mens.trans)[1] = "Pop" #rename columns
broncho.mens.pld <- cbind(broncho.data$PhyGrp, broncho.mens.trans)
names(broncho.mens.pld)[1] = "PhyGrp"
#check data
head(broncho.mens.pld)

#male data
broncho.male.raw <- subset(broncho.data.raw, select = c(6, 12:37))
colnames(broncho.male.raw)[1] <- "Pop"
broncho.male.sub <- broncho.male.raw[broncho.male.raw$Sex=="m",]
broncho.male.mens <- subset(broncho.male.sub, select = c(1, 4:17))
broncho.male.trans <- allom(broncho.male.mens, type = "population1")
broncho.mens.male <- cbind(broncho.male.sub$PhyGrp, broncho.male.trans)
names(broncho.mens.male)[1] = "PhyGrp"
head(broncho.mens.male)

#female data
broncho.female.raw <- subset(broncho.data.raw, select = c(6, 12:37))
colnames(broncho.female.raw)[1] <- "Pop"
broncho.female.sub <- broncho.female.raw[broncho.male.raw$Sex=="f",]
broncho.female.mens <- subset(broncho.female.sub, select = c(1, 4:17))
broncho.female.trans <- allom(broncho.female.mens, type = "population1")
broncho.mens.female <- cbind(broncho.female.sub$PhyGrp, broncho.female.trans)
names(broncho.mens.female)[1] = "PhyGrp"
head(broncho.mens.female)
```

**2. Perform Flexible Discriminant Analysis (FDA) with mensural dataset**
```{r, perform mensural FDA }
set.seed(122588)
#perform FDA
broncho.fda.pld <- fda(Pop ~., data = broncho.mens.trans) 
broncho.fda.male <- fda(Pop ~., data = broncho.male.trans)
broncho.fda.female <- fda(Pop ~., data = broncho.female.trans)

#generate fda confusion matrix for combined data
fda.con.mat <- table(broncho.mens.trans$Pop, predict(broncho.fda.pld))
fda.con.mat
write.csv(fda.con.mat, "fda.mens.pld.csv")

#check accuracy of classification
tbl <- table(broncho.mens.trans$Pop, predict(broncho.fda.pld))
sum(diag(tbl))/sum(tbl)

#generate fda confusion matrix for male data
fda.con.mat <- table(broncho.male.trans$Pop, predict(broncho.fda.male))
fda.con.mat
write.csv(fda.con.mat, "fda.mens.male.csv")

#check accuracy of classification
tbl <- table(broncho.male.trans$Pop, predict(broncho.fda.male))
sum(diag(tbl))/sum(tbl)

#generate fda confusion matrix for female data
fda.con.mat <- table(broncho.female.trans$Pop, predict(broncho.fda.female))
fda.con.mat
write.csv(fda.con.mat, "fda.mens.female.csv")

#check accuracy of classification
tbl <- table(broncho.female.trans$Pop, predict(broncho.fda.female))
sum(diag(tbl))/sum(tbl)

#convex hull
fda.axes <- predict(broncho.fda.pld, type = "var")
fda.axes2 <- as.data.frame(fda.axes)
phygrp.data <- broncho.mens.pld$PhyGrp
broncho.fda.hull <- cbind(phygrp.data, fda.axes2)

hull_data <- 
    broncho.fda.hull %>%
    drop_na() %>%
    group_by(broncho.mens.pld$PhyGrp) %>% 
    slice(chull(V1, V2))

#get axes values for customize plotting
#all data
fda.axes <- predict(broncho.fda.pld, type = "var")
fda.axes2 <- as.data.frame(fda.axes)
paic.data <- broncho.mens.pld$Pop
broncho.fda.data.pld <- cbind(paic.data, fda.axes2)

#male data
male.fda.axes <- predict(broncho.fda.male, type = "var")
male.fda.axes2 <- as.data.frame(male.fda.axes)
paic.data <- broncho.mens.male$Pop
broncho.fda.data.male <- cbind(paic.data, male.fda.axes2)

#female data
female.fda.axes <- predict(broncho.fda.female, type = "var")
female.fda.axes2 <- as.data.frame(female.fda.axes)
paic.data <- broncho.mens.female$Pop
broncho.fda.data.female <- cbind(paic.data, female.fda.axes2)

#get mensural character fda coefficients
coef.pld <- as.data.frame(broncho.fda.pld$fit$coefficients)
coef.pld <- rownames_to_column(coef.pld, var = "var")
head(coef.pld)
write.csv(coef.pld, "coef_mens_pld.csv")

coef.male <- as.data.frame(broncho.fda.male$fit$coefficients)
coef.male <- rownames_to_column(coef.male, var = "var")
head(coef.male)
write.csv(coef.male, "coef_mens_male.csv")

coef.female <- as.data.frame(broncho.fda.female$fit$coefficients)
coef.female <- rownames_to_column(coef.female, var = "var")
head(coef.female)
write.csv(coef.female, "coef_mens_female.csv")
```

**3. Prepare FDA plots for mensural data**
```{r view FDA resutls}
#assign colors to each Pop
plot_color <- c("#000000", "#408073", "#7f404a", "#7b3dcc", "#cc3d3d", "#d970cb", 
                "#8c994d","#cc9666", "#cc1489")

plot_color.m <- c("#000000", "#7f404a", "#7b3dcc", "#cc3d3d", "#d970cb", 
                "#8c994d", "#cc1489")

point_shape <-  c(8,11,7,16,10,12,3,14,15)

point_shape.m <-  c(8,7,16,10,12,3,15)

fda0 <- ggplot(broncho.fda.data.pld, aes(V1, V2)) + 
          theme_bw() +
          geom_polygon(data = hull_data, alpha = 0.3, aes(fill = factor(phygrp.data))) + 
          scale_fill_discrete(labels=c("Borneo", "Luzon", 
            "Mindanao", "Moluccas", "Mindoro", "Palawan"))

fda1 <- fda0 + geom_point(aes(color=paic.data, shape = paic.data), size = 3, fill = "#cc3d3d") + 
          scale_shape_manual(values = point_shape) + 
          scale_color_manual(values = plot_color) + theme(plot.subtitle = element_text(size=15, face = "bold", hjust = 0.5), 
            legend.key = element_blank(), 
          legend.title = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "grey"), legend.position = "right", 
          legend.key.width = unit(1, "cm"), legend.key.height = unit(0.3, "cm"), 
          legend.text = element_text(size=15), plot.title = element_text(color="black", size=17, face="bold"), 
          axis.text=element_text(size=15), axis.title=element_text(size=15)) + 
          xlab("FDA1, accuracy = 87.45%") + ylab("FDA2") + ggtitle("") + 
          guides(colour = guide_legend(override.aes = list(size=3))) + 
          scale_x_continuous(label = ~ scales::comma(.x, accuracy = 1)) + 
          labs(subtitle = "Combined") + xlim(-6, 6) + ylim(-6, 6) + ggtitle("Mensural")

fda2 <- ggplot(broncho.fda.data.male, aes(V1, V2))+ 
          geom_point(aes(color=paic.data, shape = paic.data),size= 3, fill = "#cc3d3d") + 
          theme_bw() +
          scale_shape_manual(values = point_shape.m) + scale_color_manual(values = plot_color.m) + 
          theme(plot.subtitle = element_text(size=15, face = "bold", hjust = 0.5),legend.key = element_blank(), 
          legend.title = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "grey"), legend.position = "none", 
          plot.title = element_text(color="black", size=17, face="bold"),
          axis.text=element_text(size=15), axis.title=element_text(size=15)) + 
          xlab("FDA1, accuracy = 98.98%") + ylab("FDA2") + ggtitle("") + 
          scale_x_continuous(label = ~ scales::comma(.x, accuracy = 1)) + 
          labs(subtitle = "Male") + xlim(-8, 6) + ylim(-8, 6)

fda3 <- ggplot(broncho.fda.data.female, aes(V1, V2)) + 
          theme_bw() +
          geom_point(aes(color=paic.data, shape = paic.data), size =3, fill = "#cc3d3d") + 
          scale_shape_manual(values = point_shape) + scale_color_manual(values = plot_color) + 
          theme(plot.subtitle = element_text(size=15, face = "bold", hjust = 0.5),
          legend.key = element_blank(), legend.title = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "grey"), 
          plot.title = element_text(color="black", size=17, face="bold"), legend.position = "none", 
          axis.text=element_text(size=15), axis.title=element_text(size=15)) + 
          xlab("FDA1, accuracy = 88.29%") + ylab("FDA2") + ggtitle("") + 
          scale_y_continuous(label = ~ scales::comma(.x, accuracy = 1)) + 
          labs(subtitle = "Female") + xlim(-6, 6) + ylim(-16, 6)
```
**3. FDA with meristic data**
```{r perform meristic FDA}
broncho.mer.pld <- subset(broncho.data, select = c(1,17:26))
colnames(broncho.mer.pld)[1] <- "Pop"
broncho.mer.pld2 <- subset(broncho.data, select = c(1,2,17:26)) #with phylo group 
colnames(broncho.mer.pld2)[1] <- "Pop"
broncho.mer.male <- subset(broncho.male.sub, select = c(1, 18:27))
broncho.mer.male2 <- subset(broncho.male.sub, select = c(1,2,18:27)) #with phylo group
colnames(broncho.mer.male2)[1] <- "Pop"
broncho.mer.female <- subset(broncho.female.sub, select = c(1, 18:27))
broncho.mer.female2 <- subset(broncho.female.sub, select = c(1,2,18:27)) #with phylo group
colnames(broncho.mer.female2)[1] <- "Pop"

set.seed(122588)
#perform FDA
broncho.fda.pld <- fda(Pop ~., data = broncho.mer.pld) 
broncho.fda.male <- fda(Pop ~., data = broncho.mer.male)
broncho.fda.female <- fda(Pop ~., data = broncho.mer.female)

#generate fda confusion matrix for combined data
fda.con.mat <- table(broncho.mer.pld$Pop, predict(broncho.fda.pld))
fda.con.mat
write.csv(fda.con.mat, "fda.mer.pld.csv")

#check accuracy of classification
tbl <- table(broncho.mer.pld$Pop, predict(broncho.fda.pld))
sum(diag(tbl))/sum(tbl)

#generate fda confusion matrix for male data
fda.con.mat <- table(broncho.mer.male$Pop, predict(broncho.fda.male))
fda.con.mat
write.csv(fda.con.mat, "fda.mer.male.csv")

#check accuracy of classification
tbl <- table(broncho.mer.male$Pop, predict(broncho.fda.male))
sum(diag(tbl))/sum(tbl)

#generate fda confusion matrix for female data
fda.con.mat <- table(broncho.mer.female$Pop, predict(broncho.fda.female))
fda.con.mat
write.csv(fda.con.mat, "fda.mer.female.csv")

#check accuracy of classification
tbl <- table(broncho.mer.female$Pop, predict(broncho.fda.female))
sum(diag(tbl))/sum(tbl)

#convex hull
fda.axes <- predict(broncho.fda.pld, type = "var")
fda.axes2 <- as.data.frame(fda.axes)
phygrp.data <- broncho.mer.pld2$PhyGrp
broncho.fda.hull <- cbind(phygrp.data, fda.axes2)

hull_data <- 
    broncho.fda.hull %>%
    drop_na() %>%
    group_by(broncho.mer.pld2$PhyGrp) %>% 
    slice(chull(V1, V2))

#get axes values for customize plotting
#all data
fda.axes <- predict(broncho.fda.pld, type = "var")
fda.axes2 <- as.data.frame(fda.axes)
paic.data <- broncho.mer.pld$Pop
broncho.fda.data.pld <- cbind(paic.data, fda.axes2)

#male data
male.fda.axes <- predict(broncho.fda.male, type = "var")
male.fda.axes2 <- as.data.frame(male.fda.axes)
paic.data <- broncho.mer.male$Pop
broncho.fda.data.male <- cbind(paic.data, male.fda.axes2)

#female data
female.fda.axes <- predict(broncho.fda.female, type = "var")
female.fda.axes2 <- as.data.frame(female.fda.axes)
paic.data <- broncho.mer.female$Pop
broncho.fda.data.female <- cbind(paic.data, female.fda.axes2)

#get mensural character fda coefficients
coef.pld <- as.data.frame(broncho.fda.pld$fit$coefficients)
coef.pld <- rownames_to_column(coef.pld, var = "var")
head(coef.pld)
write.csv(coef.pld, "coef_mer_pld.csv")

coef.male <- as.data.frame(broncho.fda.male$fit$coefficients)
coef.male <- rownames_to_column(coef.male, var = "var")
head(coef.male)
write.csv(coef.male, "coef_mer_male.csv")

coef.female <- as.data.frame(broncho.fda.female$fit$coefficients)
coef.female <- rownames_to_column(coef.female, var = "var")
head(coef.female)
write.csv(coef.female, "coef_mer_female.csv")
```

**4. Prepare FDA plots for meristic data**
```{r view meristic FDA resutls}
#assign colors to each Pop
plot_color <- c("#000000", "#408073", "#7f404a", "#7b3dcc", "#cc3d3d", "#d970cb", 
                "#8c994d","#cc9666", "#cc1489")

plot_color.m <- c("#000000", "#7f404a", "#7b3dcc", "#cc3d3d", "#d970cb", 
                "#8c994d", "#cc1489")

point_shape <-  c(8,11,7,16,10,12,3,14,15)

point_shape.m <-  c(8,7,16,10,12,3,15)

fda0 <- ggplot(broncho.fda.data.pld, aes(V1, V2)) + aes(fill = factor(phygrp.data)) + 
          theme_bw() +
          geom_polygon(data = hull_data, alpha = 0.3) + 
          scale_fill_discrete(labels=c("Borneo", "Luzon", 
            "Mindanao", "Moluccas", "Mindoro", "Palawan"))

fda4 <- fda0 + geom_point(aes(color=paic.data, shape = paic.data), size = 3, fill = "#cc3d3d") + 
          scale_shape_manual(values = point_shape) + 
          scale_color_manual(values = plot_color) + 
          theme(plot.subtitle = element_text(size=15, face = "bold", hjust = 0.5),legend.key = element_blank(), 
          legend.title = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "grey"), legend.position = "right", 
          legend.key.width = unit(1, "cm"), legend.key.height = unit(0.3, "cm"), 
          legend.text = element_text(size=15), plot.title = element_text(color="black", size=17, face="bold"), 
          axis.text=element_text(size=15), axis.title=element_text(size=15)) + 
          xlab("FDA1, accuracy = 81.88%") + ylab("FDA2") + ggtitle("") + 
          guides(colour = guide_legend(override.aes = list(size=3))) + 
          scale_x_continuous(label = ~ scales::comma(.x, accuracy = 1)) + 
          labs(subtitle = "Combined") + xlim(-5, 5) + ylim(-5, 5) + ggtitle("Meristic")

fda5 <- ggplot(broncho.fda.data.male, aes(V1, V2))+ geom_point(aes(color=paic.data, shape = paic.data),size= 3, 
          alpha =0.9, fill = "#cc3d3d") + 
          theme_bw() +
          scale_shape_manual(values = point_shape.m) + scale_color_manual(values = plot_color.m) + 
          theme(plot.subtitle = element_text(size=15, face = "bold", hjust = 0.5),legend.key = element_blank(), 
          legend.title = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "grey"), legend.position = "none", 
          plot.title = element_text(color="black", size=17, face="bold"),
          axis.text=element_text(size=15), axis.title=element_text(size=15)) + 
          xlab("FDA1, accuracy = 91.91%") + ylab("FDA2") + ggtitle("") + 
          scale_x_continuous(label = ~ scales::comma(.x, accuracy = 1)) + 
          labs(subtitle = "Male") + xlim(-5, 6) + ylim(-5, 5)

fda6 <- ggplot(broncho.fda.data.female, aes(V1, V2)) + 
          theme_bw() +
          geom_point(aes(color=paic.data, shape = paic.data), size =3, fill = "#cc3d3d") + 
          scale_shape_manual(values = point_shape) + scale_color_manual(values = plot_color) + 
          theme(plot.subtitle = element_text(size=15, face = "bold", hjust = 0.5),
          legend.key = element_blank(), legend.title = element_blank(), 
          panel.background = element_blank(), axis.line = element_line(colour = "grey"), 
          plot.title = element_text(color="black", size=17, face="bold"), legend.position = "none", 
          axis.text=element_text(size=15), axis.title=element_text(size=15)) + 
          xlab("FDA1, accuracy = 80.31%") + ylab("FDA2") + ggtitle("") + 
          scale_y_continuous(label = ~ scales::comma(.x, accuracy = 1)) + 
          labs(subtitle = "Female") + xlim(-6, 6) + ylim(-5, 6)
```

**5. View FDA results**
```{r fda plots, fig.width=15, fig.height=10}
ggarrange(fda1, fda2, fda3, fda4, fda5, fda6, ncol =3, nrow=2, common.legend=TRUE, legend="right", align = "hv")


#export
#combine fda results (tif format)
tiff("broncho_fda_comb.tif", res=300, width = 17, height = 10, unit="in")
ggarrange(fda1, fda2, fda3, fda4, fda5, fda6, ncol =3, nrow=2, common.legend=TRUE, legend="right", align = "hv")
dev.off()

#combine fda results (pdf format)
pdf("broncho_fda_comb.pdf", width = 17, height = 10)
ggarrange(fda1, fda2, fda3, fda4, fda5, fda6, ncol =3, nrow=2, common.legend=TRUE, legend="right", align = "hv")
dev.off()
```

**6. Perform Principal Component Analysis (PCA) with mensural dataset**
```{r perform mensural pca}
#perform PCA
broncho.pca.pld <- prcomp(broncho.mens.pld[3:16], center = TRUE, scale. = TRUE)
broncho.pca.male <- prcomp(broncho.mens.male[3:16], center = TRUE, scale. = TRUE)
broncho.pca.female <- prcomp(broncho.mens.female[3:16], center = TRUE, scale. = TRUE)

#calculate variance explained for combined data
eigs.all <- broncho.pca.pld$sdev^2
ax.1.all <- eigs.all[1]/sum(eigs.all)*100
ax.1.all.per <- sprintf(ax.1.all, fmt = '%#.2f')
ax.2.all <- eigs.all[2]/sum(eigs.all)*100
ax.2.all.per <- sprintf(ax.2.all, fmt = '%#.2f')

ax.1.all.per
ax.2.all.per

#calculate variance explained for male
eigs.male <- broncho.pca.male$sdev^2
ax.1.male <- eigs.male[1]/sum(eigs.male)*100
ax.1.male.per <- sprintf(ax.1.male, fmt = '%#.2f')
ax.2.male <- eigs.male[2]/sum(eigs.male)*100
ax.2.male.per <- sprintf(ax.2.male, fmt = '%#.2f')

ax.1.male.per
ax.2.male.per

#calculate variance explained for female
eigs.female <- broncho.pca.female$sdev^2
ax.1.female <- eigs.female[1]/sum(eigs.female)*100
ax.1.female.per <- sprintf(ax.1.female, fmt = '%#.2f')
ax.2.female <- eigs.female[2]/sum(eigs.female)*100
ax.2.female.per <- sprintf(ax.2.female, fmt = '%#.2f')

ax.1.female.per
ax.2.female.per

#or use summary function
summary(broncho.pca.pld)
summary(broncho.pca.male)
summary(broncho.pca.female)

#pca convex hull
#extract point loadings for all characters
pca_points <-  as_tibble(broncho.pca.pld$x) %>% bind_cols(broncho.mens.pld$PhyGrp) %>% bind_cols(broncho.mens.pld$Pop)
colnames(pca_points)[15] <- "phygrp"
colnames(pca_points)[16] <- "Pop"

pca_hull <- 
  pca_points %>% 
  group_by(phygrp) %>% 
  slice(chull(PC1, PC2))

#extract variable loadings for all characters
pca_load_pld <- as_tibble(broncho.pca.pld$rotation, rownames = 'variable')
pca_load_male <- as_tibble(broncho.pca.male$rotation, rownames = 'variable')
pca_load_female <- as_tibble(broncho.pca.female$rotation, rownames = 'variable')

head(pca_load_pld)
head(pca_load_male)
head(pca_load_female)

write.csv(pca_load_pld, "pca_load_pld.mens.csv")
write.csv(pca_load_male, "pca_load_male.mens.csv")
write.csv(pca_load_female, "pca_load_female.mens.csv")

#extract point loadings for mensural characters
pca_points_male <-  as_tibble(broncho.pca.male$x) %>% bind_cols(broncho.mens.male$PhyGrp) %>% bind_cols(broncho.mens.male$Pop)
colnames(pca_points_male)[15] <- "phygrp"
colnames(pca_points_male)[16] <- "Pop"

#extract point loadings for meristic characters
pca_points_female <-  as_tibble(broncho.pca.female$x) %>% bind_cols(broncho.mens.female$PhyGrp)  %>% bind_cols(broncho.mens.female$Pop)
colnames(pca_points_female)[15] <- "phygrp"
colnames(pca_points_female)[16] <- "Pop"
```

**7. Assess PCA significance with mensural data**
```{r mensural pca sig}
set.seed(122588)
#assess PC and variable signifance using PCAtest
broncho.pca.sig.pld <- PCAtest(broncho.mens.pld[3:16], 1000, 1000, 0.05, varcorr=FALSE, counter=FALSE, plot=TRUE)

broncho.pca.sig.male <- PCAtest(broncho.mens.male[3:16], 1000, 1000, 0.05, varcorr=FALSE, counter=FALSE, plot=TRUE)

broncho.pca.sig.female <- PCAtest(broncho.mens.female[3:16], 1000, 1000, 0.05, varcorr=FALSE, counter=FALSE, plot=TRUE)

#plot results using nzilbb_vowels
pca.sig.pld <- pca_test(broncho.mens.pld[3:16], n = 1000, scale = TRUE, variance_confint = .95, loadings_confint = 0.95)

p1 <- plot_variance_explained(pca.sig.pld, pc_max = NA, percent = TRUE)
p2 <- plot_loadings(pca.sig.pld, pc_no = 1, violin = FALSE, filter_boots = FALSE, quantile_threshold = FALSE)
p3 <- plot_loadings(pca.sig.pld, pc_no = 2, violin = FALSE, filter_boots = FALSE, quantile_threshold = FALSE)

pca.sig.male <- pca_test(broncho.mens.male[3:16], n = 1000, scale = TRUE, variance_confint = .95, loadings_confint = 0.95)

p4 <- plot_variance_explained(pca.sig.male, pc_max = NA, percent = TRUE)
p5 <- plot_loadings(pca.sig.male, pc_no = 1, violin = FALSE, filter_boots = FALSE, quantile_threshold = FALSE)
p6 <- plot_loadings(pca.sig.male, pc_no = 2, violin = FALSE, filter_boots = FALSE, quantile_threshold = FALSE)

pca.sig.female <- pca_test(broncho.mens.female[3:16], n = 1000, scale = TRUE, variance_confint = .95, loadings_confint = 0.95)

p7 <- plot_variance_explained(pca.sig.female, pc_max = NA, percent = TRUE)
p8 <- plot_loadings(pca.sig.female, pc_no = 1, violin = FALSE, filter_boots = FALSE, quantile_threshold = FALSE)
p9 <- plot_loadings(pca.sig.female, pc_no = 2, violin = FALSE, filter_boots = FALSE, quantile_threshold = FALSE)

ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9,ncol =3, nrow=3, common.legend=TRUE, legend="right")

tiff("pca_sig_mensural.tif", res=300, width = 17, height = 12, unit="in")
ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9,ncol =3, nrow=3, common.legend=TRUE, legend="right")
dev.off()

pdf("pca_sig_mensural.pdf",width = 17, height = 12 )
ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9,ncol =3, nrow=3, common.legend=TRUE, legend="right")
dev.off()
```

**8. Prepare PCA plots for mensural data**
```{r view mensural PCA results}
mens0 <- ggplot(pca_points, aes(x = PC1, y = PC2)) + theme_bw() +
        geom_polygon(data = pca_hull, aes(fill = factor(phygrp)), alpha = 0.3) + 
        scale_fill_discrete(labels=c("Borneo", "Luzon", 
            "Mindanao", "Moluccas", "Mindoro", "Palawan")) 
mens1 <- mens0 + geom_point(aes(color=Pop, shape = Pop), size = 3, fill = "#cc3d3d") + 
        scale_shape_manual(name = pca_points$Pop, values = point_shape) + 
        scale_color_manual(name = pca_points$Pop, values = plot_color) + labs(x = "PC1", y = "PC2") + 
        theme(legend.key=element_blank(), plot.subtitle = element_text(size=15, face = "bold", hjust = 0.5),
              plot.margin = margin(0.2, 0.5, 0, 0.2, "cm"), 
              legend.title = element_blank(), panel.background = element_blank(), 
              axis.line = element_line(colour = "grey"), legend.position = "right", 
              legend.key.width = unit(1, "cm"), legend.key.height = unit(0.3, "cm"), 
              legend.text = element_text(size=15), plot.title = element_text(color="black", size=17, face="bold"),
              axis.text=element_text(size=15), axis.title=element_text(size=15)) + 
        ggtitle("") + guides(colour = guide_legend(override.aes = list(size=3))) + 
        scale_y_continuous(label = ~ scales::comma(.x, accuracy = 1)) + 
        labs(x = "PC1, 32.07% variance explained", y = "PC2, 20.24% variance explained") + xlim(-6, 6) + ylim(-7, 6) +
        geom_segment(data = pca_load_pld, aes(x = 0, y = 0, xend = PC1*10, yend = PC2*10), 
          arrow = arrow(length = unit(1/2, 'picas'))) + 
        geom_text_repel(data = pca_load_pld, aes(x = PC1*10, y = PC2*11),
                label = pca_load_pld$variable, size = 5) + 
        labs(subtitle = "Combined") + ggtitle("Mensural")

mens2 <- ggplot(pca_points_male, aes(x = PC1, y = PC2)) + theme_bw() +
        geom_point(aes(color=Pop, shape = Pop), size = 3, fill = "#cc3d3d") + 
        scale_shape_manual(name = pca_points_male$Pop, values = point_shape.m) + 
        scale_color_manual(name = pca_points_male$Pop, values = plot_color.m) + labs(x = "PC1", y = "PC2") + 
        theme(plot.subtitle = element_text(size=15, face = "bold", hjust = 0.5), plot.margin = margin(0.2, 0.5, 0, 0.2, "cm"), 
          legend.title = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "grey"),
          legend.position = "none", legend.key.width = unit(1, "cm"), legend.key.height = unit(0.1, "cm"), 
          legend.text = element_text(size=15), plot.title = element_text(color="black", size=17, face="bold"), 
          axis.text=element_text(size=15), axis.title=element_text(size=15))  + 
        ggtitle("") + guides(colour = guide_legend(override.aes = list(size=3))) + 
        scale_x_continuous(label = ~ scales::comma(.x, accuracy = 1)) + 
        labs(x = "PC1, 32.23% variance explained", y = "PC2, 15.81% variance explained") + xlim(-6, 6) + ylim(-7, 6) +
        geom_segment(data = pca_load_male, aes(x = 0, y = 0, xend = PC1*10, yend = PC2*10), 
          arrow = arrow(length = unit(1/2, 'picas'))) + 
        geom_text_repel(data = pca_load_male, aes(x = PC1*10, y = PC2*11),
                label = pca_load_male$variable, size = 5) +
        labs(subtitle = "Male")

mens3 <- ggplot(pca_points_female, aes(x = PC1, y = PC2)) + theme_bw() +
        geom_point(aes(color=Pop, shape = Pop), size = 3, fill = "#cc3d3d") + 
        scale_shape_manual(name = pca_points_female$Pop, values = point_shape) + 
        scale_color_manual(name = pca_points_female$Pop, values = plot_color) + labs(x = "PC1", y = "PC2") + 
        theme(plot.subtitle = element_text(size=15, face = "bold", hjust = 0.5), plot.margin = margin(0.2, 0.5, 0, 0.2, "cm"), 
          legend.title = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "grey"), 
          legend.position = "none", legend.key.width = unit(1, "cm"), 
          legend.key.height = unit(0.1, "cm"), legend.text = element_text(size=12), 
          plot.title = element_text(color="black", size=17, face="bold"), axis.text=element_text(size=15), 
          axis.title=element_text(size=15))  + 
        ggtitle("") + guides(colour = guide_legend(override.aes = list(size=3))) + 
        scale_y_continuous(label = ~ scales::comma(.x, accuracy = 1)) + 
        labs(x = "PC1, 33.79% variance explained", y = "PC2, 19.82% variance explained")+ xlim(-6, 6) + ylim(-7, 6) +
        geom_segment(data = pca_load_female, aes(x = 0, y = 0, xend = PC1*10, yend = PC2*10), 
          arrow = arrow(length = unit(1/2, 'picas'))) + 
        geom_text_repel(data = pca_load_female, aes(x = PC1*10, y = PC2*11),
                label = pca_load_female$variable, size = 5) +
        labs(subtitle = "Female")
```

**9. PCA with meristic data**
```{r perform meristic pca}
#perform PCA
broncho.pca.pld <- prcomp(broncho.mer.pld[2:11], center = TRUE, scale. = TRUE)
broncho.pca.male <- prcomp(broncho.mer.male[2:11], center = TRUE, scale. = TRUE)
broncho.pca.female <- prcomp(broncho.mer.female[2:11], center = TRUE, scale. = TRUE)

#calculate variance explained for combined data
eigs.all <- broncho.pca.pld$sdev^2
ax.1.all <- eigs.all[1]/sum(eigs.all)*100
ax.1.all.per <- sprintf(ax.1.all, fmt = '%#.2f')
ax.2.all <- eigs.all[2]/sum(eigs.all)*100
ax.2.all.per <- sprintf(ax.2.all, fmt = '%#.2f')

ax.1.all.per
ax.2.all.per

#calculate variance explained for male
eigs.male <- broncho.pca.male$sdev^2
ax.1.male <- eigs.male[1]/sum(eigs.male)*100
ax.1.male.per <- sprintf(ax.1.male, fmt = '%#.2f')
ax.2.male <- eigs.male[2]/sum(eigs.male)*100
ax.2.male.per <- sprintf(ax.2.male, fmt = '%#.2f')

ax.1.male.per
ax.2.male.per

#calculate variance explained for female
eigs.female <- broncho.pca.female$sdev^2
ax.1.female <- eigs.female[1]/sum(eigs.female)*100
ax.1.female.per <- sprintf(ax.1.female, fmt = '%#.2f')
ax.2.female <- eigs.female[2]/sum(eigs.female)*100
ax.2.female.per <- sprintf(ax.2.female, fmt = '%#.2f')

ax.1.female.per
ax.2.female.per

#or use summary function
summary(broncho.pca.pld)
summary(broncho.pca.male)
summary(broncho.pca.female)

#pca convex hull
#extract point loadings for all characters
pca_points <-  as_tibble(broncho.pca.pld$x) %>% bind_cols(broncho.mer.pld2$PhyGrp) %>% bind_cols(broncho.mer.pld2$Pop)
colnames(pca_points)[11] <- "phygrp"
colnames(pca_points)[12] <- "Pop"

pca_hull <- 
  pca_points %>% 
  group_by(phygrp) %>% 
  slice(chull(PC1, PC2))

#extract variable loadings for all characters
pca_load_pld <- as_tibble(broncho.pca.pld$rotation, rownames = 'variable')
pca_load_male <- as_tibble(broncho.pca.male$rotation, rownames = 'variable')
pca_load_female <- as_tibble(broncho.pca.female$rotation, rownames = 'variable')

write.csv(pca_load_pld, "pca_load_pld.mer.csv")
write.csv(pca_load_male, "pca_load_male.mer.csv")
write.csv(pca_load_female, "pca_load_female.mer.csv")

#extract point loadings for meristic characters
pca_points_male <-  as_tibble(broncho.pca.male$x) %>% bind_cols(broncho.mer.male2$PhyGrp) %>% bind_cols(broncho.mer.male2$Pop)
colnames(pca_points_male)[11] <- "phygrp"
colnames(pca_points_male)[12] <- "Pop"

#extract point loadings for meristic characters
pca_points_female <-  as_tibble(broncho.pca.female$x) %>% bind_cols(broncho.mer.female2$PhyGrp)  %>% bind_cols(broncho.mer.female2$Pop)
colnames(pca_points_female)[11] <- "phygrp"
colnames(pca_points_female)[12] <- "Pop"
```

**10. Assess PCA significance with meristic data**
```{r meristic pca sig}
set.seed(122588)
#assess PC and variable signifance using PCAtest
broncho.pca.sig.pld <- PCAtest(broncho.mer.pld[2:11], 1000, 1000, 0.05, varcorr=FALSE, counter=FALSE, plot=TRUE)

broncho.pca.sig.male <- PCAtest(broncho.mer.male[2:11], 1000, 1000, 0.05, varcorr=FALSE, counter=FALSE, plot=TRUE)

broncho.pca.sig.female <- PCAtest(broncho.mer.female[2:11], 1000, 1000, 0.05, varcorr=FALSE, counter=FALSE, plot=TRUE)

#plot results using nzilbb_vowels
pca.sig.pld <- pca_test(broncho.mer.pld[2:11], n = 1000, scale = TRUE, variance_confint = .95, loadings_confint = 0.95)

p1 <- plot_variance_explained(pca.sig.pld, pc_max = NA, percent = TRUE)
p2 <- plot_loadings(pca.sig.pld, pc_no = 1, violin = FALSE, filter_boots = FALSE, quantile_threshold = FALSE)
p3 <- plot_loadings(pca.sig.pld, pc_no = 2, violin = FALSE, filter_boots = FALSE, quantile_threshold = FALSE)

pca.sig.male <- pca_test(broncho.mer.male[2:11], n = 1000, scale = TRUE, variance_confint = .95, loadings_confint = 0.95)

p4 <- plot_variance_explained(pca.sig.male, pc_max = NA, percent = TRUE)
p5 <- plot_loadings(pca.sig.male, pc_no = 1, violin = FALSE, filter_boots = FALSE, quantile_threshold = FALSE)
p6 <- plot_loadings(pca.sig.male, pc_no = 2, violin = FALSE, filter_boots = FALSE, quantile_threshold = FALSE)

pca.sig.female <- pca_test(broncho.mer.female[2:11], n = 1000, scale = TRUE, variance_confint = .95, loadings_confint = 0.95)

p7 <- plot_variance_explained(pca.sig.female, pc_max = NA, percent = TRUE)
p8 <- plot_loadings(pca.sig.female, pc_no = 1, violin = FALSE, filter_boots = FALSE, quantile_threshold = FALSE)
p9 <- plot_loadings(pca.sig.female, pc_no = 2, violin = FALSE, filter_boots = FALSE, quantile_threshold = FALSE)

ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9,ncol =3, nrow=3, common.legend=TRUE, legend="right")

tiff("pca_sig_meristic.tif", res=300, width = 17, height = 12, unit="in")
ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9,ncol =3, nrow=3, common.legend=TRUE, legend="right")
dev.off()

pdf("pca_sig_meristic.pdf",width = 17, height = 12 )
ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9,ncol =3, nrow=3, common.legend=TRUE, legend="right")
dev.off()
```

**11. Prepare PCA plots for meristic data**
```{r view meristic PCA results}
mer4 <- ggplot(pca_points, aes(x = PC1, y = PC2)) + theme_bw() +
        geom_polygon(data = pca_hull, aes(fill = phygrp), alpha = 0.3, show.legend = FALSE) + 
        geom_point(aes(color=Pop, shape = Pop), size = 3, fill = "#cc3d3d", show.legend = FALSE) + 
        scale_shape_manual(name = pca_points$Pop, values = point_shape) + 
        scale_color_manual(name = pca_points$Pop, values = plot_color) + labs(x = "PC1", y = "PC2") + 
        theme(legend.key=element_blank(), plot.subtitle = element_text(size=15, face = "bold", hjust = 0.5),
              plot.margin = margin(0.2, 0.5, 0, 0.2, "cm"), 
              legend.title = element_blank(), panel.background = element_blank(), 
              axis.line = element_line(colour = "grey"), 
          legend.position = "right", legend.key.width = unit(1, "cm"), legend.key.height = unit(0.3, "cm"), 
          legend.text = element_text(size=15), plot.title = element_text(color="black", size=17, face="bold"), 
          axis.text=element_text(size=15), axis.title=element_text(size=15)) + 
        ggtitle("Meristic") + guides(colour = guide_legend(override.aes = list(size=3))) + 
        scale_y_continuous(label = ~ scales::comma(.x, accuracy = 1)) + 
        labs(x = "PC1, 22.71% variance explained", y = "PC2, 14.96% variance explained") + xlim(-6, 5) + ylim(-7, 6) +
        geom_segment(data = pca_load_pld, aes(x = 0, y = 0, xend = PC1*10, yend = PC2*10), 
          arrow = arrow(length = unit(1/2, 'picas'))) + 
        geom_text_repel(data = pca_load_pld, aes(x = PC1*10, y = PC2*11),
                label = pca_load_pld$variable, size = 5) + labs(subtitle = "Combined")

mer5 <- ggplot(pca_points_male, aes(x = PC1, y = PC2)) + theme_bw() +
        geom_point(aes(color=Pop, shape = Pop), size = 3, fill = "#cc3d3d") + 
        scale_shape_manual(name = pca_points_male$Pop, values = point_shape.m) + 
        scale_color_manual(name = pca_points_male$Pop, values = plot_color.m) + labs(x = "PC1", y = "PC2") + 
        theme(plot.subtitle = element_text(size=15, face = "bold", hjust = 0.5), plot.margin = margin(0.2, 0.5, 0, 0.2, "cm"), 
          legend.title = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "grey"),
          legend.position = "none", legend.key.width = unit(1, "cm"), legend.key.height = unit(0.1, "cm"), 
          legend.text = element_text(size=15), plot.title = element_text(color="black", size=17, face="bold"), 
          axis.text=element_text(size=15), axis.title=element_text(size=15))  + 
        ggtitle("") + guides(colour = guide_legend(override.aes = list(size=3))) + 
        scale_x_continuous(label = ~ scales::comma(.x, accuracy = 1)) + 
        labs(x = "PC1, 23.64% variance explained", y = "PC2, 15.82% variance explained") + xlim(-6, 5) + ylim(-7, 6) +
        geom_segment(data = pca_load_male, aes(x = 0, y = 0, xend = PC1*10, yend = PC2*10), 
          arrow = arrow(length = unit(1/2, 'picas'))) + 
        geom_text_repel(data = pca_load_male, aes(x = PC1*10, y = PC2*10),
                label = pca_load_male$variable, size = 5) + labs(subtitle = "Male")

mer6 <- ggplot(pca_points_female, aes(x = PC1, y = PC2)) + theme_bw() +
        geom_point(aes(color=Pop, shape = Pop), size = 3, fill = "#cc3d3d") + 
        scale_shape_manual(name = pca_points_female$Pop, values = point_shape) + 
        scale_color_manual(name = pca_points_female$Pop, values = plot_color) + labs(x = "PC1", y = "PC2") + 
        theme(plot.subtitle = element_text(size=15, face = "bold", hjust = 0.5), plot.margin = margin(0.2, 0.5, 0, 0.2, "cm"), 
          legend.title = element_blank(), panel.background = element_blank(), 
          axis.line = element_line(colour = "grey"), 
          legend.position = "none", legend.key.width = unit(1, "cm"), 
          legend.key.height = unit(0.1, "cm"), legend.text = element_text(size=12), 
          plot.title = element_text(color="black", size=17, face="bold"), axis.text=element_text(size=15), 
          axis.title=element_text(size=15))  + 
        ggtitle("") + guides(colour = guide_legend(override.aes = list(size=3))) + 
        scale_y_continuous(label = ~ scales::comma(.x, accuracy = 1)) + 
        labs(x = "PC1, 22.77% variance explained", y = "PC2, 16.06% variance explained")+ xlim(-6, 5) + ylim(-7, 6) +
        geom_segment(data = pca_load_female, aes(x = 0, y = 0, xend = PC1*10, yend = PC2*10), 
          arrow = arrow(length = unit(1/2, 'picas'))) + 
        geom_text_repel(data = pca_load_female, aes(x = PC1*10, y = PC2*11),
                label = pca_load_female$variable, size = 5) + labs(subtitle = "Female")
```

**12. View PCA biplots**
```{r pca biplots,  fig.width=15, fig.height=10}
ggarrange(mens1, mens2, mens3, mer4, mer5, mer6, ncol =3, nrow=2, common.legend=TRUE, legend="right", align = "hv")

#export
#combine pca results
tiff("broncho_pca_comb.tif", res=300, width = 17, height = 10, unit="in")
ggarrange(mens1, mens2, mens3, mer4, mer5, mer6, ncol =3, nrow=2, common.legend=TRUE, legend="right", align = "hv")
dev.off()

#combine pca results (pdf format)
pdf("broncho_pca_comb.pdf", width = 17, height = 10)
ggarrange(mens1, mens2, mens3, mer4, mer5, mer6, ncol =3, nrow=2, common.legend=TRUE, legend="right", align = "hv")
dev.off()
```
