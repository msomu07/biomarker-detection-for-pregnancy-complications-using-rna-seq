BiocManager::install("edgeR", force = TRUE)
library(edgeR)
library(dplyr)
## copy previous r code and find the differentially expressed 
## genes in each trimester rather than in the entire
## preeclampsia table

third_trimester_data_df = read.table("~/Desktop/HORIZONS-Pregnancy-Data/third_trimester_data.txt", header = TRUE)
third_trimester_data = data.matrix(third_trimester_data_df)
#head(third_trimester_data)
diseases_df = read.table("~/Desktop/HORIZONS-Pregnancy-Data/third_trimester_diseases.txt", header = TRUE)
diseases = diseases_df$disease

d <- DGEList(counts=third_trimester_data,group=factor(diseases)) 
d
##phenotype_groups
##ncol(read_counts)

##filtering the data
dim(d)
d.full <- d
head(d$counts)
head(cpm(d))
apply(d$counts, 2, sum)
keep <- rowSums(cpm(d)> quantile(cpm(d), c(.10))) >= 2 
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

## 3 = normal 4 = preeclampsia
et43 <- exactTest(d1, pair=c(4, 3))
topTags(et43)
write.table(et43, "preeclamspsia_diffe_thirdT.txt")
de43 <- decideTestsDGE(et43, adjust.method="BH", p.value=0.05)
summary(de43)
de1tags43 <- rownames(d1)[as.logical(de43)] 
plotSmear(et43, de.tags=de1tags43)
abline(h = c(-2, 2), col = "blue")

## 3 = normal 1 = chronic hypertension
et13 <- exactTest(d1, pair=c(1, 3))
topTags(et13)
write.table(et13, "chronic_hypertension_diffe_thirdT.txt")
de13 <- decideTestsDGE(et13, adjust.method="BH", p.value=0.05)
summary(de13)
de1tags13 <- rownames(d1)[as.logical(de13)] 
plotSmear(et13, de.tags=de1tags13)
abline(h = c(-2, 2), col = "blue")

## 3 = normal 2 = gestational diabetes
et23 <- exactTest(d1, pair=c(2, 3))
topTags(et23)
write.table(et23, "gestational_diabetes_diffe_thirdT.txt")
de23 <- decideTestsDGE(et23, adjust.method="BH", p.value=0.05)
summary(de23)
de1tags23 <- rownames(d1)[as.logical(de23)] 
plotSmear(et23, de.tags=de1tags23)
abline(h = c(-2, 2), col = "blue")
