#=========================================================================================
# This script is use to analyze the 5xFAD and Mettl3 KO project in MagdaLab during rotation
# Issue report on zhanghouyu@cibr.ac.cn
# Copyright (c) 2022 __KoziolLab@CIBR__. All rights reserved.
#=========================================================================================

#=========================================================================================
# Project 1. WT & Mettl3 KO MS-imaging data
#=========================================================================================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("MS_Modifications.R")
setwd("2022-03-04 mettl3-KO-loxP1/")

### Step 0. Generate excel sheets to store m/z information
GenrateExcelforMSdata(SamplesNames = c("loxP1","loxP2","loxP3","KO1","KO2","KO3"),
                      BrainRegionNames = c("cerebellum","pons&medulla","midbrain",
                                           "hippocampus","thalamus","hypothalamus","fornix",
                                           "caudate putamen","basal forebrain",
                                           "ventral striatum","anterior olfactory",
                                           "olfactory bulbs","cortex","corpus callosum"),
                      SectionNames = c("01","02","03"))

### Step 1. Collapse modification intensity
CalculateModificationIntensity(MSFilePath = "0_Negative rawdata/",
                               ModificationReferenceFile = "./Modification_reference_Neg.csv")
CalculateModificationIntensity(MSFilePath = "0_Positive rawdata/",
                               ModificationReferenceFile = "./Modification_reference_Pos.csv")

### Step 2. Merge samples
MergeModificationIntensity(CSVPath = "1_Negative refined", run_RunMetaboAnalystR = F)
MergeModificationIntensity(CSVPath = "1_Positive refined", run_RunMetaboAnalystR = F)

### Step 3. Normalization
# Normalization across each slide
NormalizationIntensitySlides(MergeModificationIntensityFile = "1_Negative refined-merged.csv",
                             SectionIndex = c("-01","-02","-03","-04","-05","-06"))

NormalizationIntensitySlides(MergeModificationIntensityFile = "1_Positive refined-merged.csv",
                             SectionIndex = c("-01","-02","-03","-04","-05","-06"))

### Step 4. Pick certain brain regions and modification types for analyses

#Total brain regions and modification features for this project
BrainRegionNames <- paste0("-",c("cerebellum","pons&medulla","midbrain","hippocampus","thalamus",
                                 "hypothalamus","fornix","caudate","basal","ventral",
                                 "anterior","olfactory","cortex","corpus"))

ModificationsList <- list()
ModificationsList[["NegMode"]] <- read_csv("Modification_reference_Neg.csv", show_col_types = F) %>% pull(rNs) %>% unique()
ModificationsList[["PosMode"]] <- read_csv("Modification_reference_Pos.csv", show_col_types = F) %>% pull(rNs) %>% unique()
ModificationsList[["NegMode_methylated"]] <- c("methylated U","methylated C+hm5dC","ca5dC","m6dA","m5dC","methylated A","methylated G")
ModificationsList[["PosMode_methylated"]] <- c("methylated A","m6Am","methylated C+hm5dC","ac4C","G+8-oxo-dG","G+9-oxo-dG",
                                               "methylated G","m22G","methylated U","m5CMP","m6AMP","m5CTP","m6dA","m5dC",
                                               "ca5dC","m5dCTP","m6dATP", "f5dCTP","8-oxo-dGTP")

IntFile <- c("1_Negative refined-merged.csv",
             "1_Negative refined-merged-NormalizedSections.csv",
             "1_Negative refined-merged-NormalizedSlides.csv",
             "1_Negative refined-merged-NormalizedPercentage.csv",
             "1_Positive refined-merged.csv",
             "1_Positive refined-merged-NormalizedSections.csv",
             "1_Positive refined-merged-NormalizedSlides.csv",
             "1_Positive refined-merged-NormalizedPercentage.csv")

for (BrainRegionName in BrainRegionNames){
  cat("Processing", BrainRegionName,"...\n")
  PickBrainRegion(MergeModificationIntensityFile = IntFile[4],
                  PickBrainRegionNames = BrainRegionName,
                  run_RunMetaboAnalystR = T,
                  FeatureList = ModificationsList[["NegMode_methylated"]])
}

for (i in IntFile){RunMetaboAnalystR(pktablePath = i, rowNormMet = "NULL")}

### Step 5. Box Plot analyses for each modification in each brain region
for (i in IntFile[1:4]){
  BoxModification(MergeModificationIntensityFile = i,
                  SampleIndex = c("Neg-mettl3-loxP1-","Neg-mettl3-loxP2-","Neg-mettl3-loxP3-",
                                  "Neg-mettl3-KO1-","Neg-mettl3-KO2-","Neg-mettl3-KO3-"))
}
for (i in IntFile[5:8]){
  BoxModification(MergeModificationIntensityFile = i,
                  SampleIndex = c("Pos-mettl3-loxP1-","Pos-mettl3-loxP2-","Pos-mettl3-loxP3-",
                                  "Pos-mettl3-KO1-","Pos-mettl3-KO2-","Pos-mettl3-KO3-"))
}

### Step 6. Merge technical replicates for Box Plot analyses for each modification in each brain region
MergeTechnicalReplicates(MergeModificationIntensityFile = "1_Negative refined-merged-NormalizedPercentage.csv")
BoxModification(MergeModificationIntensityFile = "1_Negative refined-merged-NormalizedPercentage-mergeTR.csv",
                SampleIndex = c("Neg-mettl3-loxP1-","Neg-mettl3-loxP2-","Neg-mettl3-loxP3-",
                                "Neg-mettl3-KO1-","Neg-mettl3-KO2-","Neg-mettl3-KO3-"))

#Filter certain samples
MergeModificationIntensityFile = "1_Negative refined-merged-NormalizedPercentage.csv"
headers <- read.csv(MergeModificationIntensityFile, header = F, nrows = 1, as.is = T)[-1]
mergedMS <- read.csv(MergeModificationIntensityFile, header = F, row.names = 1, check.names = F, skip = 2)
colnames(mergedMS) <- headers
mergedMS1 <- mergedMS %>% select(contains("mid")) %>% select(!contains("-04"))
t.test(log2(as.numeric(mergedMS1[12,c(1:18)])), log2(as.numeric(mergedMS1[12,19:36])))
t.test(log2(as.numeric(mergedMS1[12,c(1:15)])), log2(as.numeric(mergedMS1[12,16:30])))

#=========================================================================================
# Project 2. Analyze WT & 5-FAD MS-imaging data for 6 biological replicates
#=========================================================================================
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("2022-03-02 WT vs 5-FAD mouse brain RNA modification/")
### Step 1. Collapse modification intensity
CalculateModificationIntensity(MSFilePath = "0_Negative rawdata/",
                               ModificationReferenceFile = "./Modification_reference_Neg.csv")
CalculateModificationIntensity(MSFilePath = "0_Positive rawdata/",
                               ModificationReferenceFile = "./Modification_reference_Pos.csv")

### Step 2. Merge samples
MergeModificationIntensity(CSVPath = "1_Negative refined", run_RunMetaboAnalystR = F)
MergeModificationIntensity(CSVPath = "1_Positive refined", run_RunMetaboAnalystR = F)

### Step 3. Normalization
# Normalization across each section (Not recommended!!!)
NormalizationIntensitySections(MergeModificationIntensityFile = "1_Negative refined-merged.csv",
                               SampleIndex = c("FAD1","FAD2","FAD3","FAD4","FAD5","FAD6","WT1","WT2","WT3","WT4","WT5","WT6"),
                               SectionIndex = c("-01","-02","-03"))

NormalizationIntensitySections(MergeModificationIntensityFile = "1_Positive refined-merged.csv",
                               SampleIndex = c("FAD1","FAD2","FAD3","FAD4","FAD5","FAD6","WT1","WT2","WT3","WT4","WT5","WT6"),
                               SectionIndex = c("-01","-02","-03"))

# Normalization across each slide
# Because the Sample 1,2,3 were measured on the same slide, 4,5,6 were another slide,
# We need to seperate they to P1 and P2 to do Slide-based normalization and merged them
NormalizationIntensitySlides(IntFileName = "2022-03-02 WT vs 5-FAD mouse brain RNA modification/1_Negative-merged P1.csv",
                             SectionIndex = c("-01","-02","-03"))
NormalizationIntensitySlides(IntFileName = "2022-03-02 WT vs 5-FAD mouse brain RNA modification/1_Negative-merged P2.csv",
                             SectionIndex = c("-01","-02","-03"))


# Normalization using the percentage
NormalizationIntensityPercentage(MergeModificationIntensityFile = "1_Negative refined-merged.csv",
                                 ModificationReferenceFile = "./Modification_reference_Neg.csv")

NormalizationIntensityPercentage(MergeModificationIntensityFile = "1_Positive refined-merged.csv",
                                 ModificationReferenceFile = "./Modification_reference_Pos.csv")

### Step 4. Pick certain brain regions and modification types for analyses

#Total brain regions and modification features for this project
BrainRegionNames <- paste0("-",c("cerebellum","pons&medulla","midbrain","hippocampus","thalamus",
                                 "hypothalamus","fornix","caudate putamen","basal forebrain","ventral striatum",
                                 "anterior olfactory","olfactory bulbs","cortex","corpus callosum"))

ModificationsList <- list()
ModificationsList[["NegMode"]] <- read_csv("Modification_reference_Neg.csv", show_col_types = F) %>% pull(rNs) %>% unique()
ModificationsList[["PosMode"]] <- read_csv("Modification_reference_Pos.csv", show_col_types = F) %>% pull(rNs) %>% unique()
ModificationsList[["NegMode_methylated"]] <- c("methylated A","methylated G","methylated U","methylated C+hm5dC",
                                               "m6Am","ac4C","m22G","hm5CTP","m6dA","m5dC","ca5dC","m5dCTP","m6dATP","f5dCTP",
                                               "m5CMP","m6AMP","GTP+8-oxo-dGTP","G+8-oxo-dG","m5CTP+hm5dCTP")
ModificationsList[["PosMode_methylated"]] <- c("m6Am","ac4C","m22G","m6dA","m5dC","ca5dC","m5dCTP","m6dATP",
                                               "f5dCTP","m5CMP","m6AMP")

IntFile <- c("1_Negative refined-merged.csv",
             "1_Negative refined-merged-NormalizedSections.csv",
             "1_Negative refined-merged-NormalizedSlides.csv",
             "1_Negative refined-merged-NormalizedPercentage.csv",
             "1_Positive refined-merged.csv",
             "1_Positive refined-merged-NormalizedSections.csv",
             "1_Positive refined-merged-NormalizedSlides.csv",
             "1_Positive refined-merged-NormalizedPercentage.csv")

for (BrainRegionName in BrainRegionNames){
  cat("Processing", BrainRegionName,"...\n")
  PickBrainRegion(MergeModificationIntensityFile = IntFile[4],
                  PickBrainRegionNames = BrainRegionName,
                  run_RunMetaboAnalystR = T,
                  FeatureList = ModificationsList[["NegMode"]])
}

for (i in IntFile){RunMetaboAnalystR(pktablePath = i, rowNormMet = "NULL")}

### Step 5. Box Plot analyses for each modification in each brain region
for (i in IntFile[c(1:4)]){
  BoxModification(MergeModificationIntensityFile = i,
                  SampleIndex = c("WT1-","WT2-","WT3-","WT4-","WT5-","WT6-",
                                  "5-FAD1-","5-FAD2-","5-FAD3-","5-FAD4-","5-FAD5-","5-FAD6-"))
}

for (i in IntFile[5:8]){
  BoxModification(MergeModificationIntensityFile = i,
                  SampleIndex = c("WT1-","WT2-","WT3-","WT4-","WT5-","WT6-",
                                  "5-FAD1-","5-FAD2-","5-FAD3-","5-FAD4-","5-FAD5-","5-FAD6-"))
}
