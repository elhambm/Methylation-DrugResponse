library(minfi)

idat_path_1 <- "C:/Users/beyra/Desktop/Minfi Package/illumina 850K EPIC1/Illumina 850k EPIC1/Illumina 850k EPIC1"


list.files(idat_path_1, pattern = "idat$")

rgSet_1 <- read.metharray.exp(base = idat_path_1)

#R.version.string

mars_csv <- "C:/Users/beyra/Desktop/Minfi Package/illumina 850K EPIC1/Illumina 850K EPIC1/MARS_Dataset_Preview.csv"
mars <- read.csv(mars_csv, header = TRUE, stringsAsFactors = FALSE)

mars$Sentrix_ID <- gsub("ID", "", mars$Sentrix_ID)
rg_samples <- sampleNames(rgSet_1)

mars_filtered_1 <- mars[mars$Sentrix_ID %in% rg_samples, ]
pd_1 <- pData(rgSet_1)
pd_1 <- cbind(pd_1, mars_filtered_1)
pData(rgSet_1) <- pd_1
phenoData_1 <- pData(rgSet_1)
View(phenoData_1)    



manifest_1 <- getManifest(rgSet_1)
manifest_1


head(getProbeInfo(manifest_1))


MSet_1 <- preprocessRaw(rgSet_1) 
MSet_1

head(getMeth(MSet_1)[,1:3])


head(getUnmeth(MSet_1)[,1:3])

RSet_1 <- ratioConvert(MSet_1, what = "both", keepCN = TRUE)
RSet_1


beta_1 <- getBeta(RSet_1)

# GenomicRatioSet
GRset_1 <- mapToGenome(RSet_1)
GRset_1

beta_1 <- getBeta(GRset_1)
M_1 <- getM(GRset_1)
CN_1 <- getCN(GRset_1)


sampleNames_1 <- sampleNames(GRset_1)
probeNames_1 <- featureNames(GRset_1)
pheno_1 <- pData(GRset_1)


gr_1 <- granges(GRset_1)
head(gr_1, n= 3)


annotation_1 <- getAnnotation(GRset_1)
names(annotation_1)


# Quality control
head(getMeth(MSet_1)[,1:3])

head(getUnmeth(MSet_1)[,1:3])


qc_1 <- getQC(MSet_1)
head(qc_1)

plotQC(qc_1)
title(main = "EPIC1 QC Plot")


densityPlot(MSet_1, sampGroups = phenoData_1$Sample_Group)
title(main = "Density plot_EPIC1")

densityBeanPlot(MSet_1, sampGroups = phenoData_1$Sample_Group)
title(main = "Density Bean plot_EPIC1")

#Control probes plot
controlStripPlot(rgSet_1, controls="BISULFITE CONVERSION II")

qcReport(rgSet_1, pdf= "qcReport.pdf")

#Sex prediction

predictedSex_1 <- getSex(GRset_1, cutoff = -2)
head(predictedSex_1$predictedSex)
#plotSex(GRset_1)


#Preprocessing and normalization:Designed specifically to correct bias between Type I and Type II probes on Illumina arrays.
MSet.swan_1 <- preprocessSWAN(rgSet_1)


#preprocessQuantile: Makes the entire distribution of signal intensities identical across all samples.
GRset.quantile_1 <- preprocessQuantile(rgSet_1, fixOutliers = TRUE,
                                     removeBadSamples = TRUE, badSampleCutoff = 10.5,
                                     quantileNormalize = TRUE, stratified = TRUE, 
                                     mergeManifest = FALSE, sex = NULL)
#SNPs
snps_1 <- getSnpInfo(GRset_1)
head(snps_1,10)

GRset_1 <- addSnpInfo(GRset_1)

GRset_1 <- dropLociWithSnps(GRset_1, snps=c("SBE","CpG"), maf=0)
###############################################################################################################
#Cell type composition
library(FlowSorted.Blood.EPIC)
cellCounts_1 <- estimateCellCounts2(rgSet_1,
                                    compositeCellType = "Blood",
                                    processMethod = "auto")


head(cellCounts_1)

View(cellCounts_1)

props <- cellCounts_1$prop

barplot(t(props), 
        beside = FALSE, 
        col = rainbow(ncol(props)), 
        legend.text = colnames(props),
        args.legend = list(x = "topright"),
        xlab = "Samples", 
        ylab = "Proportion", 
        main = "Cell Type Composition per Sample_EPIC1")

###############################################################################################
#extract only the intersecting probes (common features) between EPIC1 & EPIC2
shared_probes <- intersect(rownames(beta_1), rownames(beta_2))
shared_probes # only the names of the shared probes

beta_1_sub <- beta_1[shared_probes,] #creating matrix
beta_2_sub <- beta_2[shared_probes,]

cbind(beta_1_sub, beta_2_sub)

#convert column into row
cbind_t <- t(as.matrix(cbind(beta_1_sub, beta_2_sub)))
cbind_t

combined_all <- as.data.frame(cbind_t)





# IDs to keep = rownames of cbind_t
ids <- rownames(cbind_t)
ids <- unique(trimws(ids))                    # drop dups/whitespace
ids <- ids[!is.na(ids) & ids != ""]          # drop empties

# make sure the Sentrix column name matches your data ("Sentrix_ID" vs "SentrixID")
sentrix_col <- if ("Sentrix_ID" %in% names(mars)) "Sentrix_ID" else "SentrixID"

# (optional) coerce to character to avoid factor/type issues
mars[[sentrix_col]] <- as.character(mars[[sentrix_col]])

# filter mars by those IDs
mars_sub <- mars[mars[[sentrix_col]] %in% ids, , drop = FALSE]

library(dplyr)

subset <- mars_sub %>%
  select(Sentrix_ID, HMDme_00, HMDme_01, HMDme_02, HMDme_03,
         HMDme_04, HMDme_05, HMDme_06, HMDme_07, HMDme_08)

#################################################################
#################################################################
#treatment
library(tidyverse)

#NA in mars_sub replace to 0
mars_sub_1 <- subset%>%
  mutate(across(starts_with("HMDme_"), ~ as.numeric(as.character(.)))) %>%
  mutate(across(starts_with("HMDme_"), ~ replace_na(., 0)))


weeks <- 0:8

subset$response_slope <- apply(
  subset[, paste0("HMDme_0", 0:8)],
  1,
  function(x) {
    coef(lm(x ~ weeks))[2] 
  }
)



library(dplyr)
cols  <- paste0("HMDme_0", 0:8)
weeks <- 0:8

# 1)
mars_sub_1 <- subset %>%
  mutate(across(all_of(cols), ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(across(all_of(cols), ~ replace_na(., 0)))

# 2)
mars_sub_1$response_slope <- apply(
  mars_sub_1[, cols, drop = FALSE],
  1,
  function(x) coef(lm(x ~ weeks))[2]
)
















##########################################################################################
analysis_all <- list()  #A list for storing analysis results
methylation_all <- list()


drug_names <- c("TCA", "SSRI", "SNRI", "NASSA", "NARI", "SSRE", "OTH", "NL", "PP", "LT", "BZD", "SLP")


sentrix_col <- if ("Sentrix_ID" %in% names(mars_sub_1)) "Sentrix_ID" else "SentrixID"

response_df <- mars_sub_1 %>%
  dplyr::select(all_of(sentrix_col), response_slope) %>%
  distinct()

rownames(response_df) <- response_df[[sentrix_col]]




for (r in drug_names) {
  
  med <- mars_sub[, grep(r, colnames(mars_sub))] %>%
    mutate(
      NID   = mars_sub$NID,
      resp  = mars_sub$SKRSME0W,
      week  = mars_sub$Sadweeks
    )
  
  med[is.na(med)] <- 0
  med$med_in <- rowSums(med[1:9])
  
  
  filtered_med <- med[
    rowSums(med[, 7:9], na.rm = TRUE) >= 1 &
      rowSums(med[, 1:9], na.rm = TRUE) >= 4,
  ]
  
  med_nam <- rownames(filtered_med)
  med_nam  
  
  
  response_med <- response_df[response_df[[sentrix_col]] %in% med_nam, , drop = FALSE]

  methylation_all[[r]] <- combined_all[
    rownames(combined_all) %in% rownames(response_med),
  ]
  
  analysis_all[[r]] <- response_med[
    rownames(response_med) %in% rownames(methylation_all[[r]]),
  ]
  
 }



























