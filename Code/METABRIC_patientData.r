##################################################
## Create a file with all the required data for patients/tumors in the METABRIC cohort
##
## Ensure that all the required data is downloaded, the required frameworks are installed, and that the working directory is set to /path/to/ClaudinLow/
##################################################

# Load required libraries
library(gtools)
library(genefu)
library(estimate)

print("Analyzing the METABRIC cohort")

##################################################
## Read and format clinical data

# Read clinical data from Curtis et al. Nature 2012, supplementary files 2 and 3 (includes PAM50 without claudin-low)
print("Reading clinical data files")
mb_first <- read.table(file = "./Data/table_S2_revised.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
mb_second <- read.table(file = "./Data/table_S3_revised.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Bind the two data frames and sort by patient ID
curtis_mb <- rbind(mb_first, mb_second)
curtis_mb <- curtis_mb[mixedorder(curtis_mb$METABRIC_ID, decreasing = TRUE), ]
rownames(curtis_mb) <- NULL

# Read clinical data from Pereira et al. Nature Communications 2016
pereira_mb <- read.table(file = "./Data/brca_metabric/data_clinical_patient.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Read gene expression data (Note: not all patients have gene expression data)
print("Reading gene expression data")
exprs <- read.table(file = "./Data/brca_metabric/data_expression_median.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
expr_samples <- colnames(exprs)[3:length(exprs)]
expr_samples <- gsub("\\.", "-", expr_samples)

# Find which samples have data for all
allPtIDs <- intersect(intersect(expr_samples, pereira_mb$PATIENT_ID), curtis_mb$METABRIC_ID)
allPtIDs <- sort(allPtIDs)

# Clinical data frame with only patients with complete information
row.names(curtis_mb) <- curtis_mb$METABRIC_ID
curtis_mb_complete <- curtis_mb[allPtIDs, ]

row.names(pereira_mb) <- pereira_mb$PATIENT_ID
pereira_mb_complete <- pereira_mb[allPtIDs, ]


##################################################
## Correctly format gene expression data

# Remove records with no Entrez gene ID
tmp_exprs <- exprs[!(is.na(exprs$Entrez_Gene_Id)), ]
rownames(tmp_exprs) <- NULL

# Find duplicated Entrez gene IDs
duplicated_entrez_id <- tmp_exprs$Entrez_Gene_Id[duplicated(tmp_exprs$Entrez_Gene_Id)]

print(paste("There are", length(unique(duplicated_entrez_id)), "duplicated entrez IDs in the gene expression data."))

# This step is a manual curation of duplicates based on METABRIC-data downloaded on 28.03.2019 (files last modified 19.02.2019) - see the script/list for details. If different files are used, this manual curation may need to be repeated, or at least verified.
source(file = "./ReferenceFiles/RemoveDuplicateEntrezIDsMETABRIC_expression.r")

# Remove duplicated rows
tmp_exprs_noDuplicates <- tmp_exprs[-rowsToRemove, ]

duplicatedNow <- tmp_exprs_noDuplicates$Entrez_Gene_Id[duplicated(tmp_exprs_noDuplicates$Entrez_Gene_Id)]

# Check if duplicated rows successfully removed
if(length(duplicatedNow) > 0){
  stop("There are/is still ", length(duplicatedNow), " duplicated row(s). The curated list of rows to remove is incompatible with your gene expression data. See the documentation in './ReferenceFiles/RemoveDuplicateEntrezIDsMETABRIC_expression.r' and './Data/README.md'. Stopping analysis now.")
}

if(length(duplicatedNow) == 0){
  print("Duplicated rows successfully removed")
}

# Remove gene name and ID columns, and set the row names to the entrez IDs
exprs_clean <- tmp_exprs_noDuplicates[ ,-(1:2)]
row.names(exprs_clean) <- tmp_exprs_noDuplicates$Entrez_Gene_Id

# Write gene expression data with unique EntrezIDs
print("Writing gene expression data to file")
write.table(x = exprs_clean, file = "./Output/METABRIC_geneExpression_uniqueEntrez.txt", sep = "\t", quote = FALSE)

##################################################
## Claudin-low classification with the nine-cell line predictor

print("Identifying claudin-low tumors using the nine-cell line predictor (Prat et al. 2010)")

# Find entrez IDs for claudin low genes and identify those available in METABRIC data
entrezID_CLgenes <- claudinLowData$fnames
overlappingCL_entrezID <- intersect(entrezID_CLgenes, row.names(exprs_clean))

# Select relevant rows
exprs_CLGenes <- exprs_clean[row.names(exprs_clean) %in% overlappingCL_entrezID, ]

# Train centroids based on available genes
trainingData <- claudinLowData
trainingData$xd <- medianCtr(trainingData$xd)
trainingData$xd <- trainingData$xd[rownames(trainingData$xd) %in% rownames(exprs_CLGenes), ]

# Scale METABRIC data and training data
exprs_CLGenes_scaled <- t(scale(t(exprs_CLGenes), scale = TRUE, center = TRUE))
exprs_CLGenes_scaled <- exprs_CLGenes_scaled[rownames(trainingData$xd), ]

trainingData$xd <- t(scale(t(trainingData$xd), scale = TRUE, center = TRUE))

# Determine claudin-low status
cl_class <- claudinLow(x = trainingData$xd, classes = as.matrix(trainingData$classes$Group, ncol = 1), y = exprs_CLGenes_scaled, distm = "euclidean")

print("Finished identifying claudin-low tumors")

##################################################
## Integrate claudin-low classifications into clinical data table

print("Merging claudin-low classification with clinical data")

# Get the PAM50/intrinsic subtypes from Curtis et al. Nature 2012, supplementary tables 2 and 3
PAM50_subtypes <- data.frame(PAM50 = curtis_mb_complete$Pam50Subtype, stringsAsFactors = FALSE)

# Get the claudin-low status and distances to centroids from the nine-cell line predictor
nineCellLine_ClaudinLow <- cl_class$predictions
nineCellLine_ClaudinLow <- data.frame(nineCellLine_ClaudinLow,
                              CL_distance = cl_class$distances[, "euclidian distance to Claudin-low"],
                              Other_distance = cl_class$distances[, "euclidian distance to Others"])
nineCellLine_ClaudinLow$Samples <- gsub("\\.", "-", nineCellLine_ClaudinLow$Samples)
nineCellLine_ClaudinLow <- nineCellLine_ClaudinLow[mixedorder(nineCellLine_ClaudinLow$Samples, decreasing = TRUE), ]
nineCellLine_ClaudinLow <- data.frame(ClaudinLow = as.character(nineCellLine_ClaudinLow$Call),
                              CL_distance = nineCellLine_ClaudinLow$CL_distance,
                              Other_distance = nineCellLine_ClaudinLow$Other_distance,
                            row.names = nineCellLine_ClaudinLow$Samples,
                            stringsAsFactors = FALSE)
nineCellLine_ClaudinLow[nineCellLine_ClaudinLow$ClaudinLow == "Claudin", "ClaudinLow"] <- "ClaudinLow"
nineCellLine_ClaudinLow$ClaudinLow[nineCellLine_ClaudinLow$ClaudinLow == "Others"] <- PAM50_subtypes$PAM50[nineCellLine_ClaudinLow$ClaudinLow == "Others"]

# Bind the subtypes
allSubtypes <- cbind(PAM50_subtypes, nineCellLine_ClaudinLow)

# Find PAM50 subtype for claudin-low samples
PAM50ClaudinLow <- allSubtypes$ClaudinLow
PAM50ClaudinLow[PAM50ClaudinLow == "ClaudinLow"] <- paste(allSubtypes$PAM50[PAM50ClaudinLow == "ClaudinLow"], "ClaudinLow", sep = "")

# Gather all subtype information into main clinical data table
allSubtypes <- data.frame(allSubtypes, PAM50ClaudinLow = PAM50ClaudinLow)
patientData <- cbind(pereira_mb_complete, allSubtypes)

# Remove samples with NC subtype
patientData <- patientData[!(patientData$PAM50 == "NC"), ]

# Find number of claudin-low tumors in the cohort
clNum <- table(patientData$ClaudinLow)[["ClaudinLow"]]
print(paste("There are", clNum, "claudin-low tumors in the cohort."))


##################################################
## Add MKI67 to patientData

print("Finding MKI67 gene expression levels")

# EntrezID for MKI67
mki67Entrez <- 4288

# Subset and scale MKI67 gene expression
exprs_mki67 <- exprs_clean[as.character(mki67Entrez), ]
rownames(exprs_mki67) <- "MKI67"
exprs_mki67_rotated <- t(exprs_mki67)
rownames(exprs_mki67_rotated) <- gsub("\\.", "-", rownames(exprs_mki67_rotated))
exprs_mki67_selected <- data.frame(MKI67 = exprs_mki67_rotated[rownames(exprs_mki67_rotated) %in% patientData$PATIENT_ID, "MKI67"])
exprs_mki67_ordered <- exprs_mki67_selected[order(rownames(exprs_mki67_selected)), ]
MKI67 <- scale(exprs_mki67_ordered, scale = TRUE, center = TRUE)

# Add MKI67 gene expression to patientData
patientData <- cbind(patientData, MKI67)


##################################################
## Get the immune score, stromal score, and ESTIMATE score from each tumor using ESTIMATE (Yoshihara et al., Nature Communications 2013)

print("Inferring immune and stromal cell infiltration using ESTIMATE")

# Filter for common genes in ESTIMATE's gene list
filterCommonGenes(input.f = "./Output/METABRIC_geneExpression_uniqueEntrez.txt", output.f = "./Output/METABRIC_geneExpression_uniqueEntrez_filterCommonGenes.txt", id="EntrezID")

# Get ESTIMATE, immune, and stromal scores scores
estimateScore(input.ds = "./Output/METABRIC_geneExpression_uniqueEntrez_filterCommonGenes.txt", output.ds = "./Output/METABRIC_ESTIMATEscores.txt", platform="illumina")

# Read ESTIMATE output
readESTIMATE <- read.table(file = "./Output/METABRIC_ESTIMATEscores.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, skip = 2)

# Merge ESTIMATE output with all clinical Data
estimate <- data.frame(readESTIMATE[, 3:(length(readESTIMATE) - 2)], row.names = readESTIMATE$NAME)
estimate <- t(estimate)
estimate <- data.frame(estimate, PATIENT_ID = gsub("\\.", "-", rownames(estimate)), row.names = gsub("\\.", "-", rownames(estimate)))
patientData <- merge(x = patientData, y = estimate, by =  "PATIENT_ID")
rownames(patientData) <- patientData$PATIENT_ID

# Remove unnecessary output files
system("rm ./Output/METABRIC_geneExpression_uniqueEntrez_filterCommonGenes.txt")
system("rm ./Output/METABRIC_ESTIMATEscores.txt")

##################################################
## Get mutation information
print("Finding mutations")
mutations <- read.table(file = "./Data/brca_metabric/data_mutations_extended.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")

# Remove mutations found in patients without all required data
mutations <- mutations[mutations$Tumor_Sample_Barcode %in% patientData$PATIENT_ID, ]

# Add columns for TP53 and PIK3CA mutation
patientData <- data.frame(patientData, TP53_Status = "WT", PIK3CA_Status = "WT", stringsAsFactors = FALSE)

# Add TP53 mutation status to patientData
tp53MutPatients <- mutations[grep(x = mutations$Hugo_Symbol, pattern = "TP53"), "Tumor_Sample_Barcode"]
patientData[tp53MutPatients, "TP53_Status"] <- "Mut"

# Add PIK3CA mutation status to patientData
pik3caMutPatients <- mutations[grep(x = mutations$Hugo_Symbol, pattern = "PIK3CA"), "Tumor_Sample_Barcode"]
patientData[pik3caMutPatients, "PIK3CA_Status"] <- "Mut"

# Total number of mutations in tumor
for(patient in patientData$PATIENT_ID){
  patientData$MutationCounts[patientData$PATIENT_ID == patient] <- length(grep(x = mutations$Tumor_Sample_Barcode, pattern = patient))
}


##################################################
## Get gene centric copy number aberrations
# Get selected relevant genes
print("Finding selected copy number aberrations")
cna <- read.table(file = "./Data/brca_metabric/data_CNA.txt", sep = "\t", header = TRUE)

# Get selected relevant genes
myc <- cna[grep(x = cna$Hugo_Symbol, pattern = "^MYC$"), ]
erbb2 <- cna[grep(x = cna$Hugo_Symbol, pattern = "^ERBB2$"), ]
mdm4 <- cna[grep(x = cna$Hugo_Symbol, pattern = "^MDM4$"), ]

# Create data frame for selected genes
cnas <- data.frame(MYC = t(myc), ERBB2 = t(erbb2), mdm4 = t(mdm4), stringsAsFactors = FALSE)
colnames(cnas) <- c("MYC", "ERBB2", "MDM4")
cnas <- cnas[3:nrow(cnas), ]

# Create simplified row to AMP/DEL/NEUTRAL
cnas$MYC_SIMPLE[cnas$MYC == 0] <- "NEUT"
cnas$MYC_SIMPLE[cnas$MYC > 0] <- "AMP"
cnas$MYC_SIMPLE[cnas$MYC < 0] <- "DEL"

cnas$ERBB2_SIMPLE[cnas$ERBB2 == 0] <- "NEUT"
cnas$ERBB2_SIMPLE[cnas$ERBB2 > 0] <- "AMP"
cnas$ERBB2_SIMPLE[cnas$ERBB2 < 0] <- "DEL"

cnas$MDM4_SIMPLE[cnas$MDM4 == 0] <- "NEUT"
cnas$MDM4_SIMPLE[cnas$MDM4 > 0] <- "AMP"
cnas$MDM4_SIMPLE[cnas$MDM4 < 0] <- "DEL"

# Merge with patientData
cnas$PATIENT_ID <- gsub("\\.", "-", row.names(cnas))
patientData <- merge(x = patientData, y = cnas, by =  "PATIENT_ID")
rownames(patientData) <- patientData$PATIENT_ID


##################################################
## Genomic instability index

print("Calculating genomic instability index")

# Read copy number segments
readSegments <- read.table(file = "./Data/ascatSegments_withoutCNVs.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Read sampleID mapping to METABRIC-ID
patientMap <- read.table(file = "./Data/tumorIdMap.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Merge IDs
segments <- merge(x = readSegments, y = patientMap, by =  "sample")
segments <- segments[with(segments, order(metabricId, chr, start)), ]
rownames(segments) <- NULL

# Remove segments for which other data levels aren't available
segments <- segments[segments$metabricId %in% patientData$PATIENT_ID, ]

# Remove patients which do not have CNA data from patientData
patientData <- patientData[patientData$PATIENT_ID %in% segments$metabricId, ]

# Read chromosome start and ends positions
chromPos <- read.table(file = "./ReferenceFiles/ChromosomeStartEnds.txt", header = TRUE)
totalNucleotides <- sum(chromPos$end - chromPos$start)

# Calculate GII for each tumor and add to patientData, add purity and ploidy
for(pt in patientData$PATIENT_ID){
  # Find segments for the given tumor
  currentSegments <- segments[grep(x = segments$metabricId, pattern = pt), ]

  # Get segments with CNA
  cnaSegments <- currentSegments[currentSegments$nMajor != 1 | currentSegments$nMinor != 1, ]

  # Calculate GII and add to patientData
  patientData$GII[patientData$PATIENT_ID == pt] <- round(sum(cnaSegments$end - cnaSegments$start)/totalNucleotides, digits = 3)


  ## Ploidy corrected GII
  # Round ploidy to nearest whole number
  roundedPloidy <- round(x = currentSegments$ploidy[1], digits = 0)

  # Get segments with a total copy number aberrant from the ploidy
  cnaSegmentsPloidyCorrected <- currentSegments[(currentSegments$nMajor + currentSegments$nMinor) != roundedPloidy, ]

  # Calculate ploidy corrected GII and add to patientData
  patientData$GII_PloidyCorrected[patientData$PATIENT_ID == pt] <- round(sum(cnaSegmentsPloidyCorrected$end - cnaSegmentsPloidyCorrected$start)/totalNucleotides, digits = 3)


  # Add purity and ploidy to patientData
  patientData$purity[patientData$PATIENT_ID == pt] <- currentSegments[1, "purity"]
  patientData$ploidy[patientData$PATIENT_ID == pt] <- currentSegments[1, "ploidy"]
  patientData$ploidyRounded[patientData$PATIENT_ID == pt] <- roundedPloidy
}

print("Finished calculating genomic instability index")

##################################################
## Write patientData
print(paste("There are", nrow(patientData), "tumors with all required data, of which", table(patientData$ClaudinLow)[["ClaudinLow"]], "are classified as claudin-low by the nine-cell line predictor."))

print("Saving patient data to file")

write.table(x = patientData, file = "./Output/METABRIC_patientData.txt", sep = "\t", quote = FALSE)

##################################################
## Write tables with subtype distributions

# PAM50
pam50table <- data.frame(table(patientData$PAM50))
pam50table <- data.frame(Count = pam50table$Freq, row.names = pam50table$Var1)

pam50proportion <- c()
for(i in 1:nrow(pam50table)){
  prop <- pam50table[i, "Count"]/(sum(pam50table$Count))
  pam50proportion <- c(pam50proportion, prop)
}

pam50table <- data.frame(pam50table, Proportion = pam50proportion)

colnames(pam50table) <- c("Count_PAM50", "Proportion_PAM50")

write.table(x = pam50table, file = "./Output/METABRIC_PAM50_distribution.txt", sep = "\t", quote = FALSE)

# Claudin-low + PAM50
claudinLow_table <- data.frame(table(patientData$ClaudinLow))
claudinLow_table <- data.frame(Count = claudinLow_table$Freq, row.names = claudinLow_table$Var1)

claudinLow_proportion <- c()
for(i in 1:nrow(claudinLow_table)){
  prop <- claudinLow_table[i, "Count"]/(sum(claudinLow_table$Count))
  claudinLow_proportion <- c(claudinLow_proportion, prop)
}

claudinLow_table <- data.frame(claudinLow_table, Proportion = claudinLow_proportion)

colnames(claudinLow_table) <- c("Count", "Proportion")

write.table(x = claudinLow_table, file = "./Output/METABRIC_PAM50+CL_distribution.txt", sep = "\t", quote = FALSE)


# Claudin-low with underlying PAM50 subtype
distrib <- table(patientData$PAM50ClaudinLow)
clPAM <- data.frame(PAM50minusCL = c(distrib[["Basal"]],
                                     distrib[["Her2"]],
                                     distrib[["LumA"]],
                                     distrib[["LumB"]],
                                     distrib[["Normal"]]),
           CLwithUnderlyingPAM50 = c(distrib[["BasalClaudinLow"]],
                                     distrib[["Her2ClaudinLow"]],
                                     distrib[["LumAClaudinLow"]],
                                     distrib[["LumBClaudinLow"]],
                                     distrib[["NormalClaudinLow"]]),
                        row.names = rownames(pam50table))

CL_prop <- c()
for(i in 1:nrow(clPAM)){
  prop <- clPAM$CLwithUnderlyingPAM50[i]/(clPAM$PAM50minusCL[i] + clPAM$CLwithUnderlyingPAM50[i])
  CL_prop <- c(CL_prop, prop)
}

clPAM <- data.frame(clPAM, CL_prop = CL_prop, row.names = rownames(clPAM))

write.table(x = clPAM, file = "./Output/METABRIC_PAM50ClaudinLow_distribution.txt", sep = "\t", quote = FALSE)

print("Finished creating METABRIC data table. Nothing needs to be done regarding 'NAs introduced by coercion' messages.")
