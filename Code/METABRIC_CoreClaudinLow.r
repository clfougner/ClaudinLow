##################################################
## Claudin-low subtyping using a condensed gene list
##
## Requires the output from ./Code/METABRIC_patientData.r. Make sure that the path is set to /path/to/ClaudinLow/
##################################################

print("Identifiying core claudin-low tumors in the METABRIC cohort")

# Load the required packages
library(ComplexHeatmap)
library(circlize)
library(sigclust)
library(ggsci)

##################################################
## Read  data

# Read patientData
patientData <- read.table(file = "./Output/METABRIC_patientData.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Selected claudin-low related genes
clGenes <- read.table(file = "./ReferenceFiles/ClaudinLowGenes.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Colors
subtypeColors <- c("#DC362A",
                   "#ED1E78",
                   "#333A8E",
                   "#2F6DB9",
                   "#267255")

claudinLowColor <- "yellow"

uninterestingColor <- "black"


##################################################
## Correctly format the gene expression data

print("Reading gene expression data")

# Load gene expression data output from ./Code/METABRIC_patientData.r
exprs <- read.table(file = "./Output/METABRIC_geneExpression_uniqueEntrez.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

print("Formatting gene expression data")

# Subset the selected claudin-low genes
exprs_CL <- exprs[as.character(clGenes$Entrez), ]

# Rename to HUGO names
rownames(exprs_CL) <- clGenes$Hugo_Gene

# Rotate (for scaling)
exprs_rotated <- t(exprs_CL)

# Rename with "-" instead of "."
rownames(exprs_rotated) <- gsub("\\.", "-", rownames(exprs_rotated))

# Remove patients without the other required data
exprs_selected <- exprs_rotated[rownames(exprs_rotated) %in% patientData$PATIENT_ID, ]

# Order the gene expression data so that they are in the same order as ./Output/METABRIC_patientData.txt
exprs_ordered <- exprs_selected[order(rownames(exprs_selected)), ]

# Center and scale the data
exprs_centered_scaled <- scale(exprs_ordered, scale = TRUE, center = TRUE)

# Return to correct format for heatmap
exprs_centered_scaled <- t(exprs_centered_scaled)

# Write core claudin-low gene expression data to file
write.table(x = exprs_centered_scaled, file = "./Output/METABRIC_CoreClaudinLowExpression.txt", sep = "\t", row.names = TRUE, quote = FALSE)


##################################################
## Find core claudin-low samples

dists <- dist(x = t(exprs_centered_scaled), method = "euclidean")
clusterObj <- hclust(d = dists, method = "complete")
whichClust <- cutree(clusterObj, k = 2)
clClust <- whichClust[whichClust == 2]
coreCLsamples <- names(clClust)

patientData$CoreCL <- "Other"
patientData$CoreCL[patientData$PATIENT_ID %in% coreCLsamples] <- "CoreClaudinLow"

# CoreCL vs OtherCL
patientData$CoreCLOtherCL <- patientData$CoreCL
patientData$CoreCLOtherCL[patientData$ClaudinLow == "ClaudinLow" & patientData$CoreCLOtherCL != "CoreClaudinLow"] <- "OtherClaudinLow"

# CoreCL in the context of PAM50 subtypes
patientData$PAM50toCoreCL <- patientData$PAM50ClaudinLow
patientData$PAM50toCoreCL[grep(x = patientData$PAM50toCoreCL, pattern = "ClaudinLow")] <- paste(patientData$PAM50[grep(x = patientData$PAM50toCoreCL, pattern = "ClaudinLow")], "OtherClaudinLow", sep = "")
patientData$PAM50toCoreCL[patientData$CoreCL == "CoreClaudinLow"] <- paste(patientData$PAM50[patientData$CoreCL == "CoreClaudinLow"], "CoreClaudinLow", sep = "")

# Write core claudin-low to file
write.table(x = patientData, file = "./Output/METABRIC_patientData_CoreClaudinLow.txt", sep = "\t", quote = FALSE)

### Check significance of cluster using Sigclust
print("Testing significance of the CoreCL cluster using SigClust")

pval <- sigclust(t(exprs_centered_scaled), nsim = 10000)
pdf(file = "./Output/Figures/METABRIC_CoreClaudinLow_SigClust_nsim10000.pdf")
  plot(pval)
dev.off()


##################################################
## Heatmap annotations

# PAM50/Intrinsic subtypes
PAM50 <- patientData$PAM50
PAM50Colors <- c("Basal" = subtypeColors[1],
                 "Her2" = subtypeColors[2],
                 "LumA" = subtypeColors[3],
                 "LumB" = subtypeColors[4],
                 "Normal" = subtypeColors[5])
# Claudin-low
ClaudinLow <- patientData$ClaudinLow
ClaudinLow[ClaudinLow != "ClaudinLow"] <- "Other"

ClaudinLowColors <- c("Other" = uninterestingColor,
                      "ClaudinLow" = claudinLowColor)

# Core claudin-low
CoreClaudinLow <- patientData$CoreCL
CoreClaudinLowColors <- c("Other" = uninterestingColor,
                          "CoreClaudinLow" = claudinLowColor)

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

# Histology - requires data from Mukherjee et al. npj Breast Cancer 2018 (see METABRIC_patientData.r)
#Histology <- patientData$TUMOR_TYPE_SIMPLE
#
#colPal <- pal_nejm()(7)
#HistologyColors <- c("NST" = colPal[2],
#                     "MixedNST/SpecialType" = colPal[1],
#                     "Lobular" = colPal[3],
#                     "Tubular" = colPal[4],
#                     "Mucinous" = colPal[5],
#                     "Medullary-like" = colPal[6],
#                     "SpecialType" = colPal[7])

# TP53 mutation
TP53 <- patientData$TP53_Status
TP53Colors <-c("Mut" = "#792E85",
                "WT" = uninterestingColor)

# PIK3CA mutation
PIK3CA <- patientData$PIK3CA_Status
PIK3CAColors <-c("Mut" = "#792E85",
                 "WT" = uninterestingColor)

# MYC gain
MYC <- patientData$MYC_SIMPLE
MYCColors <-c("AMP" = "#792E85",
              "DEL" = uninterestingColor,
              "NEUT" = uninterestingColor)

# MDM4 gain
MDM4 <- patientData$MDM4_SIMPLE
MDM4Colors <-c("AMP" = "#792E85",
              "DEL" = uninterestingColor,
              "NEUT" = uninterestingColor)

# GII
GII <- patientData$GII_PloidyCorrected
GIIColors <- colorRamp2(breaks = c(0, 0.2, 1),
                        colors = c("#F0F724",
                                   "#CE4975",
                                   "#16198A"))

# IntClust
IntClust4 <- patientData$INTCLUST
IntClust4[IntClust4 == "4ER+"] <- "IntClust4"
IntClust4[IntClust4 == "4ER-"] <- "IntClust4"
IntClust4[IntClust4 != "IntClust4"] <- "Other"

IntClust4Colors <- c("Other" = uninterestingColor,
                     "IntClust4" = claudinLowColor)

# ER status
ER <- patientData$ER_IHC
ER[ER == "Positve"] <- "Positive"
ERColors <-c("Positive" = "#CD3434",
             "Negative" = uninterestingColor)

# HER2 gain
HER2 <- patientData$ERBB2_SIMPLE
HER2Colors <-c("AMP" = "#CD3434",
               "DEL" = uninterestingColor,
               "NEUT" = uninterestingColor)

## Gather all annotations in a data frame
annotationDF <- data.frame(PAM50 = PAM50,
                           ClaudinLow = ClaudinLow,
                           CoreClaudinLow = CoreClaudinLow,
                           ImmuneScore = ImmuneScore,
                           StromalScore = StromalScore,
                           Histology = Histology,
                           TP53 = TP53,
                           PIK3CA = PIK3CA,
                           MYC = MYC,
                           MDM4 = MDM4,
                           GII = GII,
                           IntClust4 = IntClust4,
                           ER = ER, 
                           HER2 = HER2)

## Make annotation object
heatmapAnnot <- HeatmapAnnotation(df = annotationDF,
                                  show_annotation_name = TRUE,
                                  annotation_name_side = "right",
                                  show_legend = FALSE,
                                  na_col = "black",
                                  #gap = unit(c(0.05, 0.05, 0.05, 0, 0.05, 0.05, 0, 0, 0, 0, 0.05, 0.05, 0, 0), "cm"), # If histology data is available
                                  gap = unit(c(0.05, 0.05, 0.05, 0, 0.05, 0, 0, 0, 0, 0.05, 0.05, 0, 0), "cm"),
                                  col = list(PAM50 = PAM50Colors,
                                             ClaudinLow = ClaudinLowColors,
                                             CoreClaudinLow = CoreClaudinLowColors,
                                             ImmuneScore = ImmuneScoreColors,
                                             StromalScore = StromalScoreColors,
                                             #Histology = HistologyColors, # Un-comment if histology data is available
                                             TP53 = TP53Colors,
                                             PIK3CA = PIK3CAColors,
                                             MYC = MYCColors,
                                             MDM4 = MDM4Colors,
                                             GII = GIIColors,
                                             IntClust4 = IntClust4Colors,
                                             ER = ERColors,
                                             HER2 = HER2Colors))

# Draw heatmap with core claudin-low annotation
print("Generating core claudin-low heatmap")
hm <- Heatmap(matrix = exprs_centered_scaled,
              clustering_distance_columns = "euclidean",
              clustering_method_columns = "complete",
              col = colorRamp2(c(-2.5, 0, 2.5), c("#343A94", "white", "#A20828")),
              show_column_names = FALSE,
              top_annotation = heatmapAnnot,
              show_heatmap_legend = TRUE)

pdf("./Output/Figures/METABRIC_CoreClaudinLowHeatmap.pdf", width = 9, height = 8)
  draw(hm)
dev.off()

## Make annotation object for oversized heatmap
heatmapAnnot = HeatmapAnnotation(df = annotationDF,
                       show_annotation_name = TRUE,
                       annotation_name_side = "right",
                       na_col = "black",
                       #gap = unit(c(0.05, 0.05, 0.05, 0, 0.05, 0.05, 0, 0, 0, 0, 0.05, 0.05, 0, 0), "cm"), # If histology data is available
                       gap = unit(c(0.05, 0.05, 0.05, 0, 0.05, 0, 0, 0, 0, 0.05, 0.05, 0, 0), "cm"),
                       col = list(PAM50 = PAM50Colors,
                                  ClaudinLow = ClaudinLowColors,
                                  CoreClaudinLow = CoreClaudinLowColors,
                                  ImmuneScore = ImmuneScoreColors,
                                  StromalScore = StromalScoreColors,
                                  #Histology = HistologyColors, # Un-comment if histology is available
                                  TP53 = TP53Colors,
                                  MYC = MYCColors,
                                  MDM4 = MDM4Colors,
                                  GII = GIIColors,
                                  IntClust4 = IntClust4Colors,
                                  ER = ERColors,
                                  HER2 = HER2Colors,
                                  PIK3CA = PIK3CAColors))

# Draw heatmap oversized so that annotations don't get cut off
hm <- Heatmap(matrix = exprs_centered_scaled,
                 clustering_distance_rows = "euclidean",
                 clustering_method_rows = "complete",
                 col = colorRamp2(c(-2.5, 0, 2.5), c("#343A94", "white", "#A20828")),
                 show_column_names = FALSE,
                 top_annotation = heatmapAnnot)

pdf("./Output/Figures/METABRIC_CoreClaudinLowHeatmap_oversize.pdf", width = 15, height = 15)
  draw(hm)
dev.off()

print("Finished generating core claudin-low heatmap")

