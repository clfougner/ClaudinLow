##################################################
## Analyses of the TCGA-BRCA cohort
##
## Make sure that path is set to /path/to/ClaudinLow/
##################################################

print("Analyzing the TCGA-BRCA cohort")

# Load the required libraries
library(genefu)
library(ComplexHeatmap)
library(circlize)
library(estimate)
library(sigclust)

# Colors
subtypeColorsSingle <- c("#DC362A",
                         "#ED1E78",
                         "#333A8E",
                         "#2F6DB9",
                         "#267255")

claudinLowColor <- "yellow"


##################################################
## Load and format the required data

print("Reading patient/tumor data")

# Read patientData
patientData <- read.table(file = "./Data/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)


# Read clinical sample file
histoData <- read.table(file = "./Data/brca_tcga_pan_can_atlas_2018/data_clinical_sample.txt",
                        header = TRUE,
                        sep = "\t",
                        stringsAsFactors = FALSE)

print("Reading gene expression data")

# Read gene expression data
exprs <- read.table(file = "./Data/brca_tcga_pan_can_atlas_2018/data_RNA_Seq_v2_expression_median.txt", sep = "\t", header = TRUE)

# Find sample IDs with gene expression data, re-format
exprsSamples <- substr(colnames(exprs[, 3:length(exprs)]), 1, nchar(colnames(exprs[, 3:length(exprs)])) - 3)
exprsSamples <- gsub("\\.", "-", exprsSamples)

# Remove patients without gene expression data from patientData
patientData <- patientData[patientData$PATIENT_ID %in% exprsSamples, ]
rownames(patientData) <- NULL

# Set EntrezID as row name
exprs_entrez <- data.frame(exprs[, 3:length(exprs)], row.names = exprs$Entrez_Gene_Id)

# Write gene expression data to file
print("Writing gene expression data to file")
write.table(x = exprs_entrez, file = "./Output/TCGA_geneExpression.txt", sep = "\t", quote = FALSE)

##################################################
## Claudin-low classification

print("Identifying claudin-low tumors using the nine-cell line predictor")

# Get EntrezIDs for claudin-low subtyping
entrezID_CLgenes <- claudinLowData$fnames

# Find overlapping genes
overlappingCL_entrezID <- intersect(entrezID_CLgenes, row.names(exprs_entrez))

# Subset the relevant gene expression data
exprs_CLGenes <- exprs_entrez[row.names(exprs_entrez) %in% overlappingCL_entrezID, ]

# Train centroids based on available genes
trainingData <- claudinLowData
trainingData$xd <- medianCtr(trainingData$xd)
trainingData$xd <- trainingData$xd[rownames(trainingData$xd) %in% rownames(exprs_CLGenes), ]

# Scale data and training data
exprs_CLGenes_scaled <- t(scale(t(exprs_CLGenes), scale = TRUE, center = TRUE))
exprs_CLGenes_scaled <- exprs_CLGenes_scaled[rownames(trainingData$xd), ]

trainingData$xd <- t(scale(t(trainingData$xd), scale = TRUE, center = TRUE))

# Determine claudin-low status
cl_class <- claudinLow(x = trainingData$xd, classes = as.matrix(trainingData$classes$Group, ncol = 1), y = exprs_CLGenes_scaled, distm = "euclidean")

# Add claudin-low status to patientData
patientData$ClaudinLow <- as.character(cl_class$predictions$Call)
patientData$ClaudinLow[patientData$ClaudinLow == "Claudin"] <- "ClaudinLow"
patientData$ClaudinLow[patientData$ClaudinLow == "Others"] <- "Other"

patientData$CL_Dist <- cl_class$distances[, "euclidian distance to Claudin-low"]
patientData$Other_Dist <- cl_class$distances[, "euclidian distance to Others"]

print("Finished identifying claudin-low tumors")

####################################
## ESTIMATE

print("Inferring immune and stromal infiltration using ESTIMATE")

# Filter for common genes in the ESTIMATE gene list
filterCommonGenes(input.f = "./Output/TCGA_geneExpression.txt", output.f = "./Output/TCGA_geneExpression_filterCommonGenes.txt", id="EntrezID")

# Get ESTIMATE, immune, and stromal scores scores
estimateScore(input.ds = "./Output/TCGA_geneExpression_filterCommonGenes.txt", output.ds = "./Output/TCGA_ESTIMATEscores.txt", platform="agilent")

# Read ESTIMATE output
estimate <- read.table(file = "./Output/TCGA_ESTIMATEscores.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, skip = 2)

# Add ESTIMATE output to patientData
estimate <- data.frame(estimate[, 3:(length(estimate))], row.names = estimate$NAME)
estimate <- t(estimate)
patientData <- cbind(patientData, estimate)

# Remove unnecessary output file
system("rm ./Output/TCGA_geneExpression_filterCommonGenes.txt")
system("rm ./Output/TCGA_ESTIMATEscores.txt")

####################################
## IntClust subtyping

print("Performing gene expression-based IntClust subtyping")

# Format expression data correctly
exprs_ic <- t(exprs_entrez)

# Create annotation data frame
annot_ic <- data.frame(EntrezGene.ID = exprs$Entrez_Gene_Id, Gene.Symbol = exprs$Hugo_Symbol)

# Perform IntClust subtyping
intClustSubtypes <- molecular.subtyping(sbt.model = "intClust", data = exprs_ic, annot = annot_ic, do.mapping = TRUE)

# Add IntClust subtypes to patientData
patientData <- cbind(patientData, intClust = intClustSubtypes$subtype)
patientData$IC4Simple <- patientData$intClust
patientData$IC4Simple <- as.character(patientData$IC4Simple)
patientData$IC4Simple[patientData$IC4Simple != "iC4"] <- "Other"

####################################
## Histology

# Add histological classification to patientData
for(i in 1:nrow(histoData)){
  patientData$Histology[patientData$PATIENT_ID == histoData$PATIENT_ID[i]] <- histoData$TUMOR_TYPE[otherPatientData$PATIENT_ID == histoData$PATIENT_ID[i]]
}


####################################
## Core ClaudinLow

print("Identifying CoreCL tumors using the condensed gene list")

# Read the condensed claudin-low gene list
claudinLowGenes <- read.table(file = "./ReferenceFiles/ClaudinLowGenes.txt", sep = "\t", header = TRUE)

# The dataset uses an outdated EntrezID for OCLN; fix.
claudinLowGenes$Entrez[claudinLowGenes$Hugo_Gene == "OCLN"] <- 4950

# Find entrez IDs for claudin low genes and identify those available in the dataset
entrezID_CLgenes <- claudinLowGenes$Entrez
overlappingCL_entrezID <- intersect(entrezID_CLgenes, row.names(exprs_entrez))

# Subset relevant rows
exprs_CLGenes <- exprs_entrez[row.names(exprs_entrez) %in% overlappingCL_entrezID, ]

# Scale METABRIC data and training data
exprs_CLGenes_scaled <- t(scale(t(exprs_CLGenes), scale = TRUE, center = TRUE))

# Order row names (n.b. numbers are treated as characters here, ordering is not by number value)
exprs_CLGenes_scaled <- exprs_CLGenes_scaled[order(rownames(exprs_CLGenes_scaled)), ]

# Set rownames to HUGO gene names
rownames(exprs_CLGenes_scaled) <- as.character(claudinLowGenes$Hugo_Gene[order(as.character(claudinLowGenes$Entrez))])

# Write core claudin-low gene expression data to file
write.table(x = exprs_CLGenes_scaled, file = "./Output/TCGA_CoreClaudinLowExpression.txt", sep = "\t", row.names = TRUE, quote = FALSE)


###############################################
## Heatmap annotations

# Colors
subtypeColors <- c("#DC362A",
                   "#ED1E78",
                   "#333A8E",
                   "#2F6DB9",
                   "#267255")

claudinLowColor <- "#F0F724"
uninterestingColor <- "black"

patientData$SUBTYPE[patientData$SUBTYPE == "BRCA_Basal"] <- "Basal-like"
patientData$SUBTYPE[patientData$SUBTYPE == "BRCA_Her2"] <- "HER2-enriched"
patientData$SUBTYPE[patientData$SUBTYPE == "BRCA_LumA"] <- "LumA"
patientData$SUBTYPE[patientData$SUBTYPE == "BRCA_LumB"] <- "LumB"
patientData$SUBTYPE[patientData$SUBTYPE == "BRCA_Normal"] <- "Normal-like"

# PAM50
PAM50 <- patientData$SUBTYPE
PAM50Colors <- c("Basal-like" = subtypeColors[1],
                 "HER2-enriched" = subtypeColors[2],
                 "LumA" = subtypeColors[3],
                 "LumB" = subtypeColors[4],
                 "Normal-like" = subtypeColors[5])

# Claudin Low
ClaudinLow <- patientData$ClaudinLow
ClaudinLowColors <- c("Other" = uninterestingColor,
                      "ClaudinLow" = claudinLowColor)

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

# IntClust4
IC4 <- patientData$IC4Simple
IC4Colors <- c("iC4" = claudinLowColor,
               "Other" = uninterestingColor)

# Create annotation data frame
annotationDF <- data.frame(PAM50 = PAM50,
                           ClaudinLow = ClaudinLow,
                           ImmuneScore = ImmuneScore,
                           StromalScore = StromalScore,
                           IntClust4 = IC4)

# Create annotation object
heatmapAnnot <- HeatmapAnnotation(df = annotationDF,
                                  show_annotation_name = TRUE,
                                  annotation_name_side = "right",
                                  show_legend = FALSE,
                                  gap = unit(c(0.05, 0.05, 0, 0.05), "cm"),
                                  width = 10,
                                  na_col = "black",
                                  col = list(PAM50 = PAM50Colors,
                                             ClaudinLow = ClaudinLowColors,
                                             ImmuneScore = ImmuneScoreColors,
                                             StromalScore = StromalScoreColors,
                                             IntClust4 = IC4Colors))

# Generate heatmap
hm <- Heatmap(matrix = exprs_CLGenes_scaled,
              clustering_distance_columns = "euclidean",
              clustering_method_columns = "complete",
              col = colorRamp2(c(-2.5, 0, 2.5), c("#343A94", "white", "#A20828")),
              show_column_names = FALSE,
              top_annotation = heatmapAnnot,
              show_heatmap_legend = TRUE)

pdf("./Output/Figures/TCGA_CoreClaudinLow_Heatmap.pdf", height = 5, width = 8)
  draw(hm)
dev.off()


###############################################
## Heatmap - oversize for annotation

# Generate heatmap annotation oject for oversized heatmap
heatmapAnnot <- HeatmapAnnotation(df = annotationDF,
                                  show_annotation_name = TRUE,
                                  annotation_name_side = "right",
                                  show_legend = TRUE,
                                  gap = unit(c(0.05, 0.05, 0, 0.05), "cm"),
                                  width = 10,
                                  na_col = "black",
                                  col = list(PAM50 = PAM50Colors,
                                             ClaudinLow = ClaudinLowColors,
                                             ImmuneScore = ImmuneScoreColors,
                                             StromalScore = StromalScoreColors,
                                             IntClust4 = IC4Colors))

# Generate oversized heatmap
hm <- Heatmap(matrix = exprs_CLGenes_scaled,
              clustering_distance_columns = "euclidean",
              clustering_method_columns = "complete",
              col = colorRamp2(c(-2.5, 0, 2.5), c("#343A94", "white", "#A20828")),
              show_column_names = FALSE,
              top_annotation = heatmapAnnot,
              show_heatmap_legend = TRUE)

pdf("./Output/Figures/TCGA_CoreClaudinLow_Heatmap_oversize.pdf", height = 15, width = 15)
  draw(hm)
dev.off()

print("Finished generating heatmap using the condensed claudin-low gene list")

###############################################
## Check significance of cluster using Sigclust

print("Checking significance of cluster in TCGA-BRCA")

pval <- sigclust(t(exprs_CLGenes_scaled), nsim = 10000)
pdf(file = "./Output/Figures/TCGA_CoreClaudinLow_SigClust_nsim10000.pdf")
  plot(pval)
dev.off()

###############################################
## Save patientData to file
print("Saving TCGA-BRCA patient data to file")
write.table(x = patientData, file = "./Output/TCGA-BRCA_patientData.txt", sep = "\t", quote = FALSE)

print("Finished analyses of the TCGA-BRCA cohort")