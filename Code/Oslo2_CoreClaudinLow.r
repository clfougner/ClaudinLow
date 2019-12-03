###############################################
## Analyze claudin-low and CoreCl tumors in the Oslo2 cohort.
##
## Ensure that the working directory is set to /path/to/ClaudinLow/
###############################################

# Load the required packages
library(genefu)
library(ComplexHeatmap)
library(circlize)
library(Biobase)
library(GEOquery)
library(estimate)
library(sigclust)

# Set colors
subtypeColors <- c("#DC362A",
                   "#ED1E78",
                   "#333A8E",
                   "#2F6DB9",
                   "#267255")

claudinLowColor <- "#F0F724"

uninterestingColor <- "black"

###############################################
## Download and format data

# Query expression and clinical data
print("Downloading Oslo2 data from GEO")
gset <- getGEO("GSE80999", GSEMatrix =TRUE, getGPL=TRUE)

# Build data frame with clinical data
geoID <- gset$GSE80999_series_matrix.txt.gz@experimentData@other$sample_id
geoID <- unlist(strsplit(x = geoID, split = " "))
sampleID <- gset$GSE80999_series_matrix.txt.gz@phenoData@data$title
sampleID <- gsub("BreastTumor_", "", sampleID)
erStatus <- gset$GSE80999_series_matrix.txt.gz@phenoData@data$`estrogen receptor (er) status:ch1`
her2status <- gset$GSE80999_series_matrix.txt.gz@phenoData@data$`human epidermal growth factor receptor 2 (her2) status:ch1`
pam50 <- c(gset$GSE80999_series_matrix.txt.gz@phenoData@data$`pam50 subgroup:ch1`[complete.cases(gset$GSE80999_series_matrix.txt.gz@phenoData@data$`pam50 subgroup:ch1`)],
           gset$GSE80999_series_matrix.txt.gz@phenoData@data$`pam50 subtype:ch1`[complete.cases(gset$GSE80999_series_matrix.txt.gz@phenoData@data$`pam50 subtype:ch1`)])

erStatus[erStatus == "NA"] <- NA
her2status[her2status == "NA"] <- NA
her2status[her2status == "0"] <- NA
pam50[pam50 == "NA"] <- NA


patientData <- data.frame(geoID = geoID,
                          sampleID = sampleID,
                          erStatus = erStatus,
                          her2status = her2status,
                          pam50 = pam50,
                      row.names = sampleID,
                      stringsAsFactors = FALSE)

###############################################
## Format gene expression data; find the mean gene expression value for genes represented by multiple probes

print("Formatting Oslo2 gene expression data")

# Get gene expression data
exprs <- exprs(gset$GSE80999_series_matrix.txt.gz)

# Get annotation for gene expression data
probes <- gset$GSE80999_series_matrix.txt.gz@featureData@data

# Find probes with an associated Entrez ID
probesWithEntrez <- probes[!is.na(probes$GENE), c("SPOT_ID", "GENE")]

# Remove gene expression values for probes not associated with an Entrez ID
exprsAll <- exprs[rownames(exprs) %in% probesWithEntrez$SPOT_ID, ]

# Bind probe data to gene expression data
exprsAll <- cbind(probesWithEntrez, exprsAll)

# Reset row names
rownames(exprsAll) <- NULL

# Identify duplicated genes
dups <- unique(exprsAll$GENE[duplicated(exprsAll$GENE)])

# Find number of rows in gene expresion data frame
sampleNum <- length(exprsAll)

# Create an empty matrix to place averaged values in
avgExprs <- matrix(data = NA, ncol = sampleNum, nrow = length(dups))
colnames(avgExprs) <- colnames(exprsAll)
avgExprs[, 2] <- dups

# Find the mean gene expression value for duplicated genes
for(j in 1:length(dups)){
   i <- dups[j]
   rows <- which(exprsAll$GENE == i)
   meanExprs <- colMeans(exprsAll[rows, 3:sampleNum])
   avgExprs[j, 3:ncol(avgExprs)] <- meanExprs
 }

# Bind the average values to the original data frame
exprsNoDups <- rbind(avgExprs, exprsAll)

# Remove rows with duplicate genes (averages are placed at the top of the data frame and are therefore the first occurence of a duplicated gene and are therefore kept)
exprsNoDups <- exprsNoDups[!duplicated(exprsNoDups$GENE), ]

# Remove probe and EntrezID column, set row names as EntrezID
exprsNoDups <- data.frame(exprsNoDups[, 3:length(exprsNoDups)], row.names = exprsNoDups$GENE)

# Write gene expression table (required for ESTIMATE)
print("Writing gene expression data to file")
write.table(x = exprsNoDups, file = "./Output/Oslo2_geneExpression_uniqueEntrez.txt", sep = "\t", quote = FALSE)

###############################################
## Claudin-low classification

print("Identifying claudin-low tumors using the nine-cell line predictor")

# Subset the relevant genes
entrezID_CLgenes <- claudinLowData$fnames
exprsCL <- exprsNoDups[rownames(exprsNoDups) %in% entrezID_CLgenes, ]

# Train centroids based on available genes
trainingData <- claudinLowData
trainingData$xd <- medianCtr(trainingData$xd)
trainingData$xd <- trainingData$xd[rownames(trainingData$xd) %in% rownames(exprsNoDups), ]

# Scale testing data and training data
exprs_CLGenes_scaled <- t(scale(t(data.matrix(exprsNoDups)), scale = TRUE, center = TRUE))
exprs_CLGenes_scaled <- data.matrix(exprs_CLGenes_scaled)[rownames(trainingData$xd), ]

trainingData$xd <- t(scale(t(trainingData$xd), scale = TRUE, center = TRUE))

# Classify claudin-low status
cl_class <- claudinLow(x = trainingData$xd, classes = as.matrix(trainingData$classes$Group, ncol = 1), y = exprs_CLGenes_scaled, distm = "euclidean")

# Add claudin-low status to patientData
patientData$ClaudinLow <- as.character(cl_class$predictions$Call)
patientData$ClaudinLow[patientData$ClaudinLow == "Claudin"] <- "ClaudinLow"
patientData$ClaudinLow[patientData$ClaudinLow == "Others"] <- "Other"


###############################################
## ESTIMATE

print("Inferring immune and stromal infiltration using ESTIMATE")

# Filter for common genes in ESTIMATE's gene list
filterCommonGenes(input.f = "./Output/Oslo2_geneExpression_uniqueEntrez.txt", output.f = "./Output/Oslo2_geneExpression_uniqueEntrez_filterCommonGenes.txt", id="EntrezID")

# Get ESTIMATE, immune, and stromal scores scores
estimateScore(input.ds = "./Output/Oslo2_geneExpression_uniqueEntrez_filterCommonGenes.txt", output.ds = "./Output/Oslo2_ESTIMATEscores.txt", platform="agilent")

# Read ESTIMATE output
estimate <- read.table(file = "./Output/Oslo2_ESTIMATEscores.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, skip = 2)

# Add scores to patientData
estimate <- data.frame(estimate[, 3:(length(estimate))], row.names = estimate$NAME)
estimate <- t(estimate)
patientData <- cbind(patientData, estimate)

# Remove unnecessary output files
system("rm ./Output/Oslo2_geneExpression_uniqueEntrez_filterCommonGenes.txt")
system("rm ./Output/Oslo2_ESTIMATEscores.txt")

###############################################
## IntClust Subtyping

print("IntClust based on gene expression data")

# Format expression data for IntClust
exprs_ic <- t(exprs(gset$GSE80999_series_matrix.txt.gz))

# Make annotation data frame
annot_ic <- data.frame(probe = probes$ID, EntrezGene.ID = probes$GENE, Gene.Symbol = probes$GENE_SYMBOL, row.names = probes$ID)

# IntClust subtyping
intClustSubtypes <- molecular.subtyping(sbt.model = "intClust", data = exprs_ic, annot = annot_ic)

# Add the IntClust subtypes to patientData
patientData <- cbind(patientData, intClust = intClustSubtypes$subtype)
patientData$IC4Simple <- patientData$intClust
patientData$IC4Simple <- as.character(patientData$IC4Simple)
patientData$IC4Simple[patientData$IC4Simple != "iC4"] <- "Other"
patientData$erStatus[is.na(patientData$erStatus)] <- "NA"
patientData$intClust <- as.character(patientData$intClust)
patientData$intClust[patientData$intClust == "iC4" & patientData$erStatus == "pos"] <- "iC4ER+"
patientData$intClust[patientData$intClust == "iC4" & patientData$erStatus == "neg"] <- "iC4ER-"
patientData$erStatus[patientData$erStatus == "NA"] <- NA

###############################################
## Core Claudin-low

print("Identifying core claudin-low tumors")

# Read condensed claudin-low gene list
claudinLowGenes <- read.table(file = "./ReferenceFiles/ClaudinLowGenes.txt", sep = "\t", header = TRUE)

# Subset gene expression data for the relevant genes
exprsCL_core <- exprsNoDups[rownames(exprsNoDups) %in% claudinLowGenes$Entrez, ]

# Center and scale
exprs_CLGenes_scaled <- t(scale(t(data.matrix(exprsCL_core)), scale = TRUE, center = TRUE))

# Write core claudin-low gene expression data to file
write.table(x = exprs_CLGenes_scaled, file = "./Output/Oslo2_CoreClaudinLowExpression.txt", sep = "\t", row.names = TRUE, quote = FALSE)

# Identify core core-claudin-low tumors
dists <- dist(x = t(exprs_CLGenes_scaled), method = "euclidean")
clusterObj <- hclust(d = dists, method = "complete")
whichClust <- cutree(clusterObj, k = 2)
clClust <- whichClust[whichClust == 2]
coreCLsamples <- names(clClust)

# Add core claudin-low status to patientData
patientData$CoreClaudinLow <- "Other"
patientData$CoreClaudinLow[patientData$geoID %in% coreCLsamples] <- "CoreClaudinLow"

print("Checking significance of core claudin-low cluster in the Oslo2 cohort")

# Check significance of cluster using Sigclust
pval <- sigclust(t(exprs_CLGenes_scaled), nsim = 10000)
pdf(file = "./Output/Figures/Oslo2_CoreClaudinLow_SigClust_nsim10000.pdf")
  plot(pval)
dev.off()

###############################################
## Heatmap annotations

print("Generating heatmap")

# PAM50
PAM50 <- patientData$pam50
PAM50Colors <- c("Basal" = subtypeColors[1],
                 "Her2" = subtypeColors[2],
                 "LumA" = subtypeColors[3],
                 "LumB" = subtypeColors[4],
                 "Normal" = subtypeColors[5])

# Claudin-low
ClaudinLow <- patientData$ClaudinLow
ClaudinLowColors <- c("Other" = uninterestingColor,
                      "ClaudinLow" = "yellow")

# Core claudin-low
CoreClaudinLow <- patientData$CoreClaudinLow
CoreClaudinLowColors <- c("Other" = uninterestingColor,
                      "CoreClaudinLow" = "yellow")

# Immune score
ImmuneScore <- patientData$ImmuneScore
ImmuneScoreColors <- colorRamp2(breaks = c(min(patientData$ImmuneScore),
                                           1000,
                                           max(patientData$ImmuneScore)),
                                colors = c("#430C53", "#20928C", "#F4E02B"))

# Stromal score
StromalScore <- patientData$StromalScore
StromalScoreColors <- colorRamp2(breaks = c(min(patientData$StromalScore),
                                            1000,
                                            max(patientData$StromalScore)),
                                 colors = c("#430C53", "#20928C", "#F4E02B"))

# ER
ER <- patientData$erStatus
ERColors <-c("pos" = "#CD3434",
             "neg" = uninterestingColor)

# HER2
HER2 <- patientData$her2status
HER2Colors <-c("pos" = "#CD3434",
               "neg" = uninterestingColor)

# IntClust4
IntClust4 <- patientData$IC4Simple
IntClust4Colors <- c("iC4" = "yellow",
               "Other" = uninterestingColor)

###############################################
## Heatmap
annotationDF <- data.frame(PAM50 = PAM50,
                           ClaudinLow = ClaudinLow,
                           CoreClaudinLow = CoreClaudinLow,
                           ImmuneScore = ImmuneScore,
                           StromalScore = StromalScore,
                           ER = ER,
                           HER2 = HER2,
                           IntClust4 = IntClust4)

heatmapAnnot <- HeatmapAnnotation(df = annotationDF,
                                 show_annotation_name = TRUE,
                                 annotation_name_side = "right",
                                 show_legend = FALSE,
                                 gap = unit(c(0.05, 0.05, 0.05, 0, 0.05, 0, 0.05), "cm"),
                                 width = 10,
                                 na_col = "black",
                                 col = list(PAM50 = PAM50Colors,
                                            ClaudinLow = ClaudinLowColors,
                                            CoreClaudinLow = CoreClaudinLowColors,
                                            ImmuneScore = ImmuneScoreColors,
                                            StromalScore = StromalScoreColors,
                                            ER = ERColors,
                                            HER2 = HER2Colors,
                                            IntClust4 = IntClust4Colors))

hm <- Heatmap(matrix = exprs_CLGenes_scaled,
              clustering_distance_columns = "euclidean",
              clustering_method_columns = "complete",
              col = colorRamp2(c(-2.5, 0, 2.5), c("#343A94", "white", "#A20828")),
              show_column_names = FALSE,
              top_annotation = heatmapAnnot,
              show_heatmap_legend = TRUE)
pdf("./Output/Figures/Oslo2_CoreClaudinLow_Heatmap.pdf", height = 5, width = 8)
  draw(hm)
dev.off()

###############################################
## Heatmap - oversize for annotation

heatmapAnnot <- HeatmapAnnotation(df = annotationDF,
                                  show_annotation_name = TRUE,
                                  annotation_name_side = "right",
                                  show_legend = TRUE,
                                  gap = unit(c(0.05, 0.05, 0.05, 0, 0.05, 0, 0.05), "cm"),
                                  width = 10,
                                  na_col = "black",
                                  col = list(PAM50 = PAM50Colors,
                                             ClaudinLow = ClaudinLowColors,
                                             CoreClaudinLow = CoreClaudinLowColors,
                                             ImmuneScore = ImmuneScoreColors,
                                             StromalScore = StromalScoreColors,
                                             ER = ERColors,
                                             HER2 = HER2Colors,
                                             IntClust4 = IntClust4Colors))

hm <- Heatmap(matrix = exprs_CLGenes_scaled,
              clustering_distance_columns = "euclidean",
              clustering_method_columns = "complete",
              col = colorRamp2(c(-2.5, 0, 2.5), c("#343A94", "white", "#A20828")),
              show_column_names = FALSE,
              top_annotation = heatmapAnnot,
              show_heatmap_legend = TRUE)
pdf("./Output/Figures/Oslo2_CoreClaudinLow_Heatmap_oversize.pdf", height = 15, width = 15)
  draw(hm)
dev.off()

print("Finished generating CoreCL heatmap")

###############################################
## Save patientData to file

# PAM50 to claudin-low
patientData$PAM50toCL <- patientData$pam50
patientData$PAM50toCL[patientData$ClaudinLow == "ClaudinLow"] <- paste(patientData$PAM50toCL[patientData$ClaudinLow == "ClaudinLow"], "ClaudinLow", sep ="")

# Claudin-low vs other claudin-low
patientData$CoreCLOtherCL <- patientData$CoreClaudinLow
patientData$CoreCLOtherCL[patientData$ClaudinLow == "ClaudinLow" & patientData$CoreCLOtherCL != "CoreClaudinLow"] <- "OtherClaudinLow"

# PAM50 to core claudin-low
patientData$PAM50toCoreCL <- patientData$pam50
patientData$PAM50toCoreCL[patientData$ClaudinLow == "ClaudinLow"] <- paste(patientData$pam50[patientData$ClaudinLow == "ClaudinLow"], "OtherClaudinLow", sep = "")
patientData$PAM50toCoreCL[patientData$CoreClaudinLow == "CoreClaudinLow"] <- paste(patientData$pam50[patientData$CoreClaudinLow == "CoreClaudinLow"], "CoreClaudinLow", sep = "")

# Wtrite patientData to file
print("Saving Oslo2 patient data to file")
write.table(x = patientData, file = "./Output/Oslo2_patientData.txt", sep = "\t", quote = FALSE)
