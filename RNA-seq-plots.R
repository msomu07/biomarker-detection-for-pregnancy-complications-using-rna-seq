#BiocManager::install("edgeR", force = TRUE)
library(edgeR)
library(dplyr)
## read in datasets (txt files)
phenotypes_df = read.table("~/Desktop/HORIZONS-Pregnancy-Data/phenotypes.txt", header = TRUE)
read_counts_df = read.table("~/Desktop/HORIZONS-Pregnancy-Data/read_count_data.txt", header = TRUE)
read_counts = data.matrix(read_counts_df)
head(read_counts_df)
phenotype_groups = phenotypes_df$phenotype

d <- DGEList(counts=read_counts,group=factor(phenotype_groups)) ## number of groups 
d
##phenotype_groups
##ncol(read_counts)

##filtering the data
dim(d)
d.full <- d
head(d$counts)
head(cpm(d))
apply(d$counts, 2, sum)
keep <- rowSums(cpm(d)> quantile(cpm(d), c(.10))) >= 2 ## figure out a more accurate threshold
d <- d[keep,]
dim(d)

##mean(rowSums(cpm(d)))
##quantile(cpm(d), c(.10))

d$samples$lib.size <- colSums(d$counts)
d$samples

##normalizing the data
d <- calcNormFactors(d)
d

## data exploration
plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)
## drop outliers

#estimating the dispersion
d1 <- estimateCommonDisp(d, verbose=T)
names(d1)
d1 <- estimateTagwiseDisp(d1)
names(d1)
plotBCV(d1)

## creating a generalized linear model fit
design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="power") 
## can change power
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)

## ------------------- EVERYTHING VS NORMAL ----------------------

## 1 = hypertension 4 = normal
et14 <- exactTest(d1, pair=c(1, 4))
topTags(et14)
write.table(et14, "normal_vs_hypertension.txt")
de14 <- decideTestsDGE(et14, adjust.method="BH", p.value=0.05)
summary(de14)
de1tags14 <- rownames(d1)[as.logical(de14)] 
plotSmear(et14, de.tags=de1tags14)
abline(h = c(-2, 2), col = "blue")

## 2 = diabetes 4 = normal 
et24 <- exactTest(d1, pair=c(2, 4))
topTags(et24)
write.table(et24, "normal_vs_diabetes.txt")
de24 <- decideTestsDGE(et24, adjust.method="BH", p.value=0.05)
summary(de24)
de1tags24 <- rownames(d1)[as.logical(de24)] 
plotSmear(et24, de.tags=de1tags24)
abline(h = c(-2, 2), col = "blue")

## 5 = preeclampsia 4 = normal
et54 <- exactTest(d1, pair=c(5, 4))
topTags(et54)
write.table(et54, "normal_vs_preeclampsia.txt")
de54 <- decideTestsDGE(et54, adjust.method="BH", p.value=0.05)
summary(de54)
de1tags54 <- rownames(d1)[as.logical(de54)] 
plotSmear(et54, de.tags=de1tags54)
abline(h = c(-2, 2), col = "blue")

## 3 = non pregnant 4 = normal
et34 <- exactTest(d1, pair=c(3, 4))
topTags(et34)
write.table(et34, "normal_vs_nonpregnant.txt")
de34 <- decideTestsDGE(et34, adjust.method="BH", p.value=0.05)
summary(de34)
de1tags34 <- rownames(d1)[as.logical(de34)] 
plotSmear(et34, de.tags=de1tags34)
abline(h = c(-2, 2), col = "blue")


## ------------------- EVERYTHING VS NON-PREGNANT ------------------- 


## 1 = hypertension 3 = non pregnant
et13 <- exactTest(d1, pair=c(1, 3))
topTags(et13)
write.table(et13, "nonpregnant_vs_hypertension.txt")
de13 <- decideTestsDGE(et13, adjust.method="BH", p.value=0.05)
summary(de13)
de1tags13 <- rownames(d1)[as.logical(de13)] 
plotSmear(et13, de.tags=de1tags13)
abline(h = c(-2, 2), col = "blue")

## 2 = diabetes 3 = non pregnant 
et23 <- exactTest(d1, pair=c(2, 3))
topTags(et23)
write.table(et23, "nonpregnant_vs_diabetes.txt")
de23 <- decideTestsDGE(et23, adjust.method="BH", p.value=0.05)
summary(de23)
de1tags23 <- rownames(d1)[as.logical(de23)] 
plotSmear(et23, de.tags=de1tags23)
abline(h = c(-2, 2), col = "blue")

## 5 = preeclampsia 3 = non pregnant 
et53 <- exactTest(d1, pair=c(5, 3))
topTags(et53)
write.table(et53, "normal_vs_preeclampsia.txt")
de53 <- decideTestsDGE(et53, adjust.method="BH", p.value=0.05)
summary(de53)
de1tags53 <- rownames(d1)[as.logical(de53)] 
plotSmear(et53, de.tags=de1tags53)
abline(h = c(-2, 2), col = "blue")

## ------------------- DISEASES VS DISEASES ----------------------

## 2 = diabetes 1 = hypertension 
et21 <- exactTest(d1, pair=c(2, 1))
topTags(et21)
write.table(et21, "hypertension_vs_diabetes.txt")
de21 <- decideTestsDGE(et21, adjust.method="BH", p.value=0.05)
summary(de21)
de1tags21 <- rownames(d1)[as.logical(de21)] 
plotSmear(et21, de.tags=de1tags21)
abline(h = c(-2, 2), col = "blue")

## 5 = preeclampsia 1 = hypertension 
et51 <- exactTest(d1, pair=c(5, 1))
topTags(et51)
write.table(et51, "hypertension_vs_preeclampsia.txt")
de51 <- decideTestsDGE(et51, adjust.method="BH", p.value=0.05)
summary(de51)
de1tags51 <- rownames(d1)[as.logical(de51)] 
plotSmear(et51, de.tags=de1tags51)
abline(h = c(-2, 2), col = "blue")

## 5 = preeclampsia 2 = diabetes  
et52 <- exactTest(d1, pair=c(5, 2))
topTags(et52)
write.table(et52, "diabetes_vs_preeclampsia.txt")
de52 <- decideTestsDGE(et52, adjust.method="BH", p.value=0.05)
summary(de52)
de1tags52 <- rownames(d1)[as.logical(de52)] 
plotSmear(et52, de.tags=de1tags52)
abline(h = c(-2, 2), col = "blue")

