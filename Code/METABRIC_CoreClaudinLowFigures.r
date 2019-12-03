##################################################
## Core claudin-low figures
##
## Requires the output from ./Code/METABRIC_patientData.r (./Output/METABRIC_geneExpression_uniqueEntrez.txt) and METABRIC_CoreClaudinLow.r (METABRIC_patientData_coreClaudinLow.txt). Ensure that the working directory is set to /path/to/ClaudinLow/
##################################################

print("Generating figures for core claudin-low tumors")

# Load the required packages
library(survival)
library(survminer)
library(ggplot2)
library(ggsignif)
library(gridExtra)

##################################################
## Read required data

# Read patient data 
patientData <- read.table(file = "./Output/METABRIC_patientData_CoreClaudinLow.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# PAM50 colors
subtypeColorsSingle <- c("#DC362A",
                         "#ED1E78",
                         "#333A8E",
                         "#2F6DB9",
                         "#267255")

# Claudin-low colors
coreCLColor <- "gold"
otherCLColor <- "darkorange2"
allOthersColor <- "#F2F2F2"
claudinColors <- c(otherCLColor, coreCLColor)

# Get condensed claudin-low gene list
clGenes <- read.table(file = "./ReferenceFiles/ClaudinLowGenes.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)


##################################################
## Load and format gene expression data

print("Reading gene expression data")
exprs <- read.table(file = "./Output/METABRIC_geneExpression_uniqueEntrez.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

print("Formatting gene expression data")
exprs_CL <- exprs[as.character(clGenes$Entrez), ]
rownames(exprs_CL) <- clGenes$Hugo_Gene
exprs_rotated <- t(exprs_CL)
rownames(exprs_rotated) <- gsub("\\.", "-", rownames(exprs_rotated))
exprs_selected <- exprs_rotated[rownames(exprs_rotated) %in% patientData$PATIENT_ID, ]
exprs_ordered <- exprs_selected[order(rownames(exprs_selected)), ]
exprs_centered_scaled <- scale(exprs_ordered, scale = TRUE, center = TRUE)

# Add gene expression for the selected genes to patientData
patientDataExprs <- cbind(patientData, exprs_centered_scaled)
patientDataExprs$CoreCLOtherCL <- factor(patientDataExprs$CoreCLOtherCL, levels = c("CoreClaudinLow", "OtherClaudinLow", "Other"))

##################################################
## Create gene expression boxplots for CoreCL and OtherCL

pdf("./Output/Figures/METABRIC_CoreClaudinLow_ExpressionBoxplots_CoreCLOtherCL.pdf", pointsize = 7, height = 15, width = 12)
p <- list()
  for (i in (length(patientDataExprs) - nrow(clGenes) + 1):length(patientDataExprs)){
    currentGene <- as.character(colnames(patientDataExprs)[i])
    idx <- i - (length(patientDataExprs) - nrow(clGenes))
    
    minExpr <- min(eval(parse(text = paste("patientDataExprs$", currentGene, sep = ""))))
    maxExpr <- max(eval(parse(text = paste("patientDataExprs$", currentGene, sep = ""))))
    
    p[[idx]] <- ggplot(patientDataExprs, aes_string(x = "CoreCLOtherCL",
                                        y = currentGene)) + 
      geom_boxplot(aes(fill = CoreCLOtherCL)) +
      theme_bw() +
      theme(legend.position = "none",
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black")) +
      scale_fill_manual(values = c("CoreClaudinLow" = coreCLColor,
                                   "OtherClaudinLow" = otherCLColor,
                                   "Other" = allOthersColor)) +
      labs(title = currentGene, x = "") +
      ylim(minExpr - 0.5, maxExpr + 3.5) +
      geom_signif(comparisons = list(c("Other", "OtherClaudinLow"),
                                     c("Other", "CoreClaudinLow"),
                                     c("OtherClaudinLow", "CoreClaudinLow")), 
                  map_signif_level=TRUE,
                  y_position = c(maxExpr + 1,
                                 maxExpr + 2,
                                 maxExpr + 3.1))
  }

grid.arrange(grobs = p, ncol = 5)

dev.off()

##################################################
## Gene expression plots for basal-like only

patientDataExprs$PAM50toCoreCL<- factor(patientDataExprs$PAM50toCoreCL , levels = levels(as.factor(patientDataExprs$PAM50toCoreCL))[c(2, 3, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)])

pdf("./Output/Figures/METABRIC_CoreClaudinLow_ExpressionBoxplots_BasalLike.pdf", pointsize = 7, height = 15, width = 12)
p <- list()
  for (i in (length(patientDataExprs) - nrow(clGenes) + 1):length(patientDataExprs)){
    currentGene <- as.character(colnames(patientDataExprs)[i])
    idx <- i - (length(patientDataExprs) - nrow(clGenes))
    
    minExpr <- min(eval(parse(text = paste("patientDataExprs$", currentGene, sep = ""))))
    maxExpr <- max(eval(parse(text = paste("patientDataExprs$", currentGene, sep = ""))))
    
    p[[idx]] <- ggplot(patientDataExprs[patientDataExprs$PAM50toCoreCL == "BasalCoreClaudinLow" | patientDataExprs$PAM50toCoreCL == "BasalOtherClaudinLow" | patientDataExprs$PAM50toCoreCL == "Basal", ], aes_string(x = "PAM50toCoreCL",
                                            y = currentGene)) + 
      geom_boxplot(aes(fill = PAM50toCoreCL)) +
      theme_bw() +
      theme(legend.position = "none",
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            axis.text.x = element_text(angle = 30, hjust = 1, color = "black"),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black")) +
      theme(plot.margin=unit(c(0.1, 0.1, 0, 1),"cm")) +
      scale_fill_manual(values = c("BasalCoreClaudinLow" = coreCLColor,
                                   "BasalOtherClaudinLow" = otherCLColor,
                                   "Basal" = subtypeColorsSingle[1])) +
      labs(title = currentGene, x = "") +
      ylim(minExpr - 0.5, maxExpr + 3.5) +
      geom_signif(comparisons = list(c("Basal", "BasalOtherClaudinLow"),
                                     c("Basal", "BasalCoreClaudinLow"),
                                     c("BasalOtherClaudinLow", "BasalCoreClaudinLow")), 
                  map_signif_level=TRUE,
                  y_position = c(maxExpr + 1,
                                 maxExpr + 2.2,
                                 maxExpr + 3.2))
    
  }

grid.arrange(grobs = p, ncol = 5)

dev.off()

##################################################
## Generate Fig. 4
pdf("./Output/Figures/METABRIC_CoreClaudinLow_Fig4.pdf", width = 6, height = 4, pointsize = 10)
par(mar=c(4, 4, 3, 3), mfrow = c(2, 3))

##################################################
## Distribution of subtypes

# Get number of tumors by subtype and claudin-low status
coreCLDistribution <- table(patientData$PAM50[patientData$CoreCLOtherCL == "CoreClaudinLow"])
coreCLDistributionCounts <- data.frame(CoreClaudinLow = c(coreCLDistribution[["Basal"]],
                                                          coreCLDistribution[["Her2"]],
                                                          coreCLDistribution[["LumA"]],
                                                          0, #tumorDistribution[["LumB"]], # No CoreCL LumBs
                                                          coreCLDistribution[["Normal"]]), 
                                       row.names = c("Basal",
                                                     "Her2",
                                                     "LumA",
                                                     "LumB",
                                                     "Normal"))

otherCLDistribution <- table(patientData$PAM50[patientData$CoreCLOtherCL == "OtherClaudinLow"])
otherCLDistributionCounts <- data.frame(OtherClaudinLow = c(otherCLDistribution[["Basal"]],
                                                            0, #otherCLDistribution[["Her2"]], # No OtherCL HER2Es
                                                            otherCLDistribution[["LumA"]],
                                                            otherCLDistribution[["LumB"]],
                                                            otherCLDistribution[["Normal"]]), 
                                        row.names = c("Basal",
                                                      "Her2",
                                                      "LumA",
                                                      "LumB",
                                                      "Normal"))

clDistribution <- as.matrix(cbind(coreCLDistributionCounts, otherCLDistributionCounts))

# Make plot
barplot(clDistribution,
        col = subtypeColorsSingle,
        horiz = TRUE,
        xlab = "Number of tumors",
        xlim = c(0, 80),
        las = 1,
        axisnames = TRUE,
        axes = TRUE,
        legend = TRUE,
        args.legend = list(x = "top",
                           bty = "n",
                           ncol = 5,
                           inset = c(0, -1),
                           cex = 0.75,
                           xpd = NA))

##################################################
## Genomic instability index

# Select relevant data
currentData <- data.frame(patientID = patientData$PATIENT_ID, GII = patientData$GII_PloidyCorrected, ClaudinLow = patientData$PAM50toCoreCL, PAM50 = patientData$PAM50)

# Remove LumB, LumA, HER2-enriched and normal-like
currentData <- currentData[!currentData$PAM50 == "LumB", ]
currentData <- currentData[!currentData$PAM50 == "LumA", ]
currentData <- currentData[!currentData$PAM50 == "Her2", ]
currentData <- currentData[!currentData$PAM50 == "Normal", ]

# Drop unused levels
currentData <- droplevels(currentData)

# Set as factor so that the subtypes are ordered correctly in the figure
currentData$ClaudinLow <- factor(currentData$ClaudinLow , levels = levels(as.factor(currentData$ClaudinLow))[c(2, 3, 1)])

# Make plot
bp <- boxplot(GII ~ ClaudinLow,
              data = currentData,
              ylab = "Genomic instability index",
              xlab = "",
              ylim = c(0, 1),
              main = "Genomic instability index",
              outline = FALSE,
              las = 2,
              col = c(subtypeColorsSingle[1],
                      subtypeColorsSingle[1],
                      subtypeColorsSingle[1]),
              frame = FALSE)

rect(xleft = c(0.6),
     ybottom = bp$stats[2, c(1)],
     xright = c(1.4),
     ytop = bp$stats[4, c(1)],
     angle = 60,
     density = 20,
     col = "yellow",
     border = "black")

rect(xleft = c(1.6),
     ybottom = bp$stats[2, c(2)],
     xright = c(2.4),
     ytop = bp$stats[4, c(2)],
     angle = 60,
     density = 5,
     col = "orange",
     border = "black")

rect(xleft = c(0.6, 1.6),
     ybottom = bp$stats[3, c(1, 2)],
     xright = c(1.4, 2.4),
     ytop = bp$stats[3, c(1, 2)],
     col = "black",
     lwd = 2)

##################################################
## Mutation count plot for Normal, Basal

# Subset relevant data
currentData <- data.frame(patientID = patientData$PATIENT_ID, MutationCounts = patientData$MutationCounts, ClaudinLow = patientData$PAM50toCoreCL, PAM50 = patientData$PAM50)

# Remove LumB, LumA, HER2-enriched and normal-like
currentData <- currentData[!currentData$PAM50 == "LumB", ]
currentData <- currentData[!currentData$PAM50 == "LumA", ]
currentData <- currentData[!currentData$PAM50 == "Her2", ]
currentData <- currentData[!currentData$PAM50 == "Normal", ]

# Drop unused levels
currentData <- droplevels(currentData)

# Set levels so that subtypes are ordered correctly
currentData$ClaudinLow <- factor(currentData$ClaudinLow , levels = levels(as.factor(currentData$ClaudinLow))[c(2, 3, 1)])

# Make plots
bp <- boxplot(log2(MutationCounts + 1) ~ ClaudinLow,
              data = currentData,
              ylab = "Log2(Mutations + 1)",
              xlab = "",
              ylim = c(1, 5),
              main = "Mutations",
              outline = FALSE,
              las = 2,
              col = c(subtypeColorsSingle[1],
                      subtypeColorsSingle[1],
                      subtypeColorsSingle[1]),
              frame = FALSE)

rect(xleft = c(0.6),
     ybottom = bp$stats[2, c(1)],
     xright = c(1.4),
     ytop = bp$stats[4, c(1)],
     angle = 60,
     density = 20,
     col = "yellow",
     border = "black")

rect(xleft = c(1.6),
     ybottom = bp$stats[2, c(2)],
     xright = c(2.4),
     ytop = bp$stats[4, c(2)],
     angle = 60,
     density = 5,
     col = "orange",
     border = "black")

rect(xleft = c(0.6, 1.6),
     ybottom = bp$stats[3, c(1, 2)],
     xright = c(1.4, 2.4),
     ytop = bp$stats[3, c(1, 2)],
     col = "black",
     lwd = 2)

##################################################
## IntClust4 barplot

# Get proportions in IntClust4
interesting <- c("BasalCoreClaudinLow", "BasalOtherClaudinLow", "Basal")
allProps <- c()
for(subtype in interesting){
  totalPatients <- sum(table(patientData$INTCLUST[patientData$PAM50toCoreCL == subtype]))
    pos <- table(patientData$INTCLUST[patientData$PAM50toCoreCL == subtype])[["4ER+"]] + table(patientData$INTCLUST[patientData$PAM50toCoreCL == subtype])[["4ER-"]]
  proportion <- pos/totalPatients
  allProps <- c(allProps, proportion)
}

names(allProps) <- interesting

# Make plot
barplot(allProps, 
        col = c(subtypeColorsSingle[1],
                subtypeColorsSingle[1],
                subtypeColorsSingle[1]),
        ylab = "Proportion",
        ylim = c(0, 1),
        main = "IntClust4",
        las = 2)


barplot(allProps,
        add = TRUE,
        las = 2,
        angle = 60,
        density = c(20, 0, 0),
        col = "yellow")

barplot(allProps,
        add = TRUE,
        las = 2,
        angle = 60,
        density = c(0, 5, 0),
        col = "orange")

##################################################
## MKI67

# Select relevant data
currentData <- data.frame(patientID = patientData$PATIENT_ID, MKI67 = patientData$MKI67, ClaudinLow = patientData$PAM50toCoreCL, PAM50 = patientData$PAM50)

# Remove LumB, LumA, HER2-enriched and normal-like
currentData <- currentData[!currentData$PAM50 == "LumB", ]
currentData <- currentData[!currentData$PAM50 == "LumA", ]
currentData <- currentData[!currentData$PAM50 == "Her2", ]
currentData <- currentData[!currentData$PAM50 == "Normal", ]

# Drop unused levels
currentData <- droplevels(currentData)

# Set as factor so that the subtypes are ordered correctly in the figure
currentData$ClaudinLow <- factor(currentData$ClaudinLow , levels = levels(as.factor(currentData$ClaudinLow))[c(2, 3, 1)])

# Make plot
bp <- boxplot(MKI67 ~ ClaudinLow,
              data = currentData,
              ylab = "MKI67 gene expression",
              xlab = "",
              main = "MKI67",
              outline = FALSE,
              las = 2,
              col = c(subtypeColorsSingle[1],
                      subtypeColorsSingle[1],
                      subtypeColorsSingle[1]),
              frame = FALSE)

rect(xleft = c(0.6),
     ybottom = bp$stats[2, c(1)],
     xright = c(1.4),
     ytop = bp$stats[4, c(1)],
     angle = 60,
     density = 20,
     col = "yellow",
     border = "black")

rect(xleft = c(1.6),
     ybottom = bp$stats[2, c(2)],
     xright = c(2.4),
     ytop = bp$stats[4, c(2)],
     angle = 60,
     density = 5,
     col = "orange",
     border = "black")

rect(xleft = c(0.6, 1.6),
     ybottom = bp$stats[3, c(1, 2)],
     xright = c(1.4, 2.4),
     ytop = bp$stats[3, c(1, 2)],
     col = "black",
     lwd = 2)

dev.off()

##################################################
## Kaplan-Meier plot

# Create a data frame for the relevant data
survivalDF_root <- data.frame(sample = patientData$PATIENT_ID,
                              subtype = patientData$PAM50toCoreCL,
                              survival = patientData$OS_MONTHS,
                              status = patientData$OS_STATUS,
                              DiseaseSpecificSurvival = patientData$VITAL_STATUS)


# Subset basal-like and claudin-low
survivalDF <- survivalDF_root[grep(x = survivalDF_root$subtype, pattern = "Basal"), ]

# Make Kaplan-Meier plot
survivalDF$SurvObj <- Surv(survivalDF$survival, event = survivalDF$DiseaseSpecificSurvival== "Died of Disease")
survivalByType <- survfit(SurvObj ~ subtype, data = survivalDF)
survivalPlotBasal <- ggsurvplot(survivalByType,
                                data = survivalDF,
                                palette = c(subtypeColorsSingle[1],
                                            coreCLColor,
                                            otherCLColor),
                                pval = TRUE,
                                legend = "bottom",
                                censor = FALSE,
                                xlim = c(0, 200),
                                break.x.by = 50,
                                xlab = "Months",
                                ylab = "",
                                title = "Basal-like")

pdf("./Output/Figures/METABRIC_CoreClaudinLow_basalSurvival.pdf", width = 4, height = 4, onefile = FALSE) 
  print(survivalPlotBasal)
dev.off()


##################################################
## TP53, MYC, MDM4 and PIK3CA mutations and amplifications

# Select subtypes
interesting <- c("BasalCoreClaudinLow", "BasalOtherClaudinLow", "Basal")

pdf("./Output/Figures/METABRIC_CoreClaudinLow_MutationsAmplifications.pdf", width = 6, height = 2.6, pointsize = 10)
par(mar=c(10,4,4,2), mfrow = c(1, 4))

##### TP53
# Find proportions
allProps <- c()
for(subtype in interesting){
  totalPatients <- sum(table(patientData$TP53_Status[patientData$PAM50toCoreCL == subtype]))
  pos <- table(patientData$TP53_Status[patientData$PAM50toCoreCL == subtype])[["Mut"]]
  proportion <- pos/totalPatients
  allProps <- c(allProps, proportion)
}
names(allProps) <- interesting

# Make plot
barplot(allProps, 
        col = c(subtypeColorsSingle[1],
                subtypeColorsSingle[1],
                subtypeColorsSingle[1]),
        ylab = "Proportion with mutation",
        ylim = c(0, 1),
        main = "TP53",
        las = 2)

barplot(allProps,
        add = TRUE,
        las = 2,
        angle = 60,
        density = c(20, 0, 0),
        col = "yellow")

barplot(allProps,
        add = TRUE,
        las = 2,
        angle = 60,
        density = c(0, 5, 0),
        col = "orange")

###### MYC
# Find proportions
allProps <- c()
for(subtype in interesting){
  totalPatients <- sum(table(patientData$MYC_SIMPLE[patientData$PAM50toCoreCL == subtype]))
  pos <- table(patientData$MYC_SIMPLE[patientData$PAM50toCoreCL == subtype])[["AMP"]]
  proportion <- pos/totalPatients
  allProps <- c(allProps, proportion)
}
names(allProps) <- interesting

# Make plot
barplot(allProps, 
        col = c(subtypeColorsSingle[1],
                subtypeColorsSingle[1],
                subtypeColorsSingle[1]),
        ylab = "Proportion with gain",
        ylim = c(0, 0.8),
        main = "MYC",
        las = 2)


barplot(allProps,
        add = TRUE,
        las = 2,
        angle = 60,
        density = c(20, 0, 0),
        col = "yellow")

barplot(allProps,
        add = TRUE,
        las = 2,
        angle = 60,
        density = c(0, 5, 0),
        col = "orange")

##### MDM4
# Find proportions
allProps <- c()
for(subtype in interesting){
  totalPatients <- sum(table(patientData$MDM4_SIMPLE[patientData$PAM50toCoreCL == subtype]))
  pos <- table(patientData$MDM4_SIMPLE[patientData$PAM50toCoreCL == subtype])[["AMP"]]
  proportion <- pos/totalPatients
  allProps <- c(allProps, proportion)
}
names(allProps) <- interesting

# Make plot
barplot(allProps, 
        col = c(subtypeColorsSingle[1],
                subtypeColorsSingle[1],
                subtypeColorsSingle[1]),
        ylab = "Proportion with gain",
        ylim = c(0, 0.6),
        main = "MDM4",
        las = 2)


barplot(allProps,
        add = TRUE,
        las = 2,
        angle = 60,
        density = c(20, 0, 0),
        col = "yellow")

barplot(allProps,
        add = TRUE,
        las = 2,
        angle = 60,
        density = c(0, 5, 0),
        col = "orange")

##### PIK3CA
# Find proportions
allProps <- c()
for(subtype in interesting){
  totalPatients <- sum(table(patientData$PIK3CA_Status[patientData$PAM50toCoreCL == subtype]))
  pos <- table(patientData$PIK3CA_Status[patientData$PAM50toCoreCL == subtype])[["Mut"]]
  proportion <- pos/totalPatients
  allProps <- c(allProps, proportion)
}
names(allProps) <- interesting

# Make plots
barplot(allProps, 
        col = c(subtypeColorsSingle[1],
                subtypeColorsSingle[1],
                subtypeColorsSingle[1]),
        ylab = "Proportion with mutation",
        ylim = c(0, 0.25),
        main = "PIK3CA",
        las = 2)


barplot(allProps,
        add = TRUE,
        las = 2,
        angle = 60,
        density = c(20, 0, 0),
        col = "yellow")

barplot(allProps,
        add = TRUE,
        las = 2,
        angle = 60,
        density = c(0, 5, 0),
        col = "orange")
dev.off()

#####################################################
## StromalScore and ImmuneScore

pdf("./Output/Figures/METABRIC_CoreClaudinLow_StromalImmuneScore.pdf", width = 4, height = 3, pointsize = 8)
par(mar=c(10,4,4,2), mfrow = c(1, 2))

##### ImmuneScore
# Subset relevant data
currentData <- data.frame(patientID = patientData$PATIENT_ID, Score = patientData$ImmuneScore, ClaudinLow = patientData$PAM50toCoreCL, PAM50 = patientData$PAM50)

# Remove LumB, LumA, HER2-enriched and normal-like tumors
currentData <- currentData[!currentData$PAM50 == "LumB", ]
currentData <- currentData[!currentData$PAM50 == "LumA", ]
currentData <- currentData[!currentData$PAM50 == "Her2", ]
currentData <- currentData[!currentData$PAM50 == "Normal", ]

# Drop unused levels
currentData <- droplevels(currentData)

# Factor so that subtypes come in the wanted order
currentData$ClaudinLow <- factor(currentData$ClaudinLow , levels = levels(as.factor(currentData$ClaudinLow))[c(2, 3, 1)])


bp <- boxplot(Score ~ ClaudinLow,
              data = currentData,
              ylab = "Immune score",
              ylim = c(-1000, 4000),
              xlab = "",
              main = "Immune score",
              outline = FALSE,
              las = 2,
              col = c(subtypeColorsSingle[1],
                      subtypeColorsSingle[1],
                      subtypeColorsSingle[1]),
              frame = FALSE)

rect(xleft = c(0.6),
     ybottom = bp$stats[2, c(1)],
     xright = c(1.4),
     ytop = bp$stats[4, c(1)],
     angle = 60,
     density = 20,
     col = "yellow",
     border = "black")

rect(xleft = c(1.6),
     ybottom = bp$stats[2, c(2)],
     xright = c(2.4),
     ytop = bp$stats[4, c(2)],
     angle = 60,
     density = 5,
     col = "orange",
     border = "black")

rect(xleft = c(0.6, 1.6),
     ybottom = bp$stats[3, c(1, 2)],
     xright = c(1.4, 2.4),
     ytop = bp$stats[3, c(1, 2)],
     col = "black",
     lwd = 2)

##### StromalScore
# Subset relevant data
currentData <- data.frame(patientID = patientData$PATIENT_ID, Score = patientData$StromalScore, ClaudinLow = patientData$PAM50toCoreCL, PAM50 = patientData$PAM50)

# Remove LumB, LumA, HER2-enriched and normal-like tumors
currentData <- currentData[!currentData$PAM50 == "LumB", ]
currentData <- currentData[!currentData$PAM50 == "LumA", ]
currentData <- currentData[!currentData$PAM50 == "Her2", ]
currentData <- currentData[!currentData$PAM50 == "Normal", ]

# Drop unused levels
currentData <- droplevels(currentData)

# Factor so that subtypes come in the wanted order
currentData$ClaudinLow <- factor(currentData$ClaudinLow , levels = levels(as.factor(currentData$ClaudinLow))[c(2, 3, 1)])

# Make plot
bp <- boxplot(Score ~ ClaudinLow,
              data = currentData,
              ylab = "Stromal score",
              main = "Stromal score",
              ylim = c(-1000, 2000),
              xlab = "",
              outline = FALSE,
              las = 2,
              col = c(subtypeColorsSingle[1],
                      subtypeColorsSingle[1],
                      subtypeColorsSingle[1]),
              frame = FALSE)

rect(xleft = c(0.6),
     ybottom = bp$stats[2, c(1)],
     xright = c(1.4),
     ytop = bp$stats[4, c(1)],
     angle = 60,
     density = 20,
     col = "yellow",
     border = "black")

rect(xleft = c(1.6),
     ybottom = bp$stats[2, c(2)],
     xright = c(2.4),
     ytop = bp$stats[4, c(2)],
     angle = 60,
     density = 5,
     col = "orange",
     border = "black")

rect(xleft = c(0.6, 1.6),
     ybottom = bp$stats[3, c(1, 2)],
     xright = c(1.4, 2.4),
     ytop = bp$stats[3, c(1, 2)],
     col = "black",
     lwd = 2)

dev.off()

print("Finished generating core claudin-low figures")