library(minfi)

idat_Path_2 <- "C:/Users/beyra/Desktop/Minfi Package/illumina 850K EPIC2/Illumina 850k EPIC2/Illumina 850k EPIC2"

list.files(idat_Path_2, pattern = "idat$")

rgSet_2 <- read.metharray.exp(base = idat_Path_2)



mars_csv <- "C:/Users/beyra/Desktop/Minfi Package/illumina 850K EPIC1/Illumina 850K EPIC1/MARS_Dataset_Preview.csv"
mars <- read.csv(mars_csv, header = TRUE, stringsAsFactors = FALSE)

mars$Sentrix_ID <- gsub("ID", "", mars$Sentrix_ID)
rg_samples <- sampleNames(rgSet_2)

mars_filtered_2 <- mars[mars$Sentrix_ID %in% rg_samples, ]
pd_2 <- pData(rgSet_2)
pd_2 <- cbind(pd_2, mars_filtered_2)
pData(rgSet_2) <- pd_2
phenoData_2 <- pData(rgSet_2)
View(phenoData_2)    


phenoData_2 <- pData(rgSet_2)
phenoData_2[,1:6]

manifest_2 <- getManifest(rgSet_2)
manifest_2

head(getProbeInfo(manifest_2))
getProbeInfo(manifest_2)


MSet_2 <- preprocessRaw(rgSet_2) 
MSet_2

head(getMeth(MSet_2)[,1:3])
getMeth(MSet_2)[,1:3]
getMeth(MSet_2)

head(getUnmeth(MSet_2)[,1:3])


RSet_2 <- ratioConvert(MSet_2, what = "both", keepCN = TRUE)
RSet_2

beta_2 <- getBeta(RSet_2)



# GenomicRatioSet
GRset_2 <- mapToGenome(RSet_2)
GRset_2

beta_2 <- getBeta(GRset_2)
M_2 <- getM(GRset_2)
CN_2 <- getCN(GRset_2)

sampleNames_2 <- sampleNames(GRset_2)
probeNames_2 <- featureNames(GRset_2)
pheno_2 <- pData(GRset_2)

gr_2 <- granges(GRset_2)
head(gr_2, n= 3)

annotation_2 <- getAnnotation(GRset_2)
names(annotation_2)


head(getMeth(MSet_2)[,1:3])

head(getUnmeth(MSet_2)[,1:3])


qc_2 <- getQC(MSet_2)
head(qc_2)

plotQC(qc_2)
title(main = "EPIC2 QC Plot")

densityPlot(MSet_2, sampGroups = phenoData_2$Sample_Group)


densityBeanPlot(MSet_2, sampGroups = phenoData_2$Sample_Group)


#Control probes plot
controlStripPlot(rgSet_2, controls="BISULFITE CONVERSION II")

qcReport(rgSet_2, pdf= "qcReport.pdf")

#Sex prediction

predictedSex_2 <- getSex(GRset_2, cutoff = -2)
head(predictedSex$predictedSex)
#plotSex(GRset_2)


#Preprocessing and normalization
MSet.swan_2 <- preprocessSWAN(rgSet_2)


#preprocessQuantile
GRset.quantile_2 <- preprocessQuantile(rgSet_2, fixOutliers = TRUE,
                                     removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                     quantileNormalize = TRUE, stratified = TRUE, 
                                     mergeManifest = FALSE, sex = NULL)





#SNPs
snps_2 <- getSnpInfo(GRset_2)
head(snps_2,10)

GRset_2 <- addSnpInfo(GRset_2)

GRset_2 <- dropLociWithSnps(GRset_2, snps=c("SBE","CpG"), maf=0)

#Cell type composition
library(FlowSorted.Blood.EPIC)
cellCounts_2 <- estimateCellCounts2(rgSet_2,
                                  compositeCellType = "Blood",
                                  processMethod = "auto")


head(cellCounts_2)

View(cellCounts_2)

props_2 <- cellCounts_2$prop

barplot(t(props_2), 
        beside = FALSE, 
        col = rainbow(ncol(props_2)), 
        legend.text = colnames(props_2),
        args.legend = list(x = "topright"),
        xlab = "Samples", 
        ylab = "Proportion", 
        main = "Cell Type Composition per Sample_EPIC2")


