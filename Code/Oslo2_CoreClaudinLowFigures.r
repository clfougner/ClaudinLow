##################################################
## Figures for the Oslo2 cohort
##
## Requires the output from ./Code/Oslo2_CoreClaudinLow.r (./Output/Oslo2_patientData.txt). Make sure that the path is set to /path/to/ClaudinLow/
##################################################

print("Generating other figures for the Oslo2 cohort")

# Read patientData
patientData <- read.table(file = "./Output/Oslo2_patientData.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Subtype colors, basal-like, ERBB2+, LumA, LumB, normal-like
subtypeColorsSingle <- c("#DC362A",
                         "#ED1E78",
                         "#333A8E",
                         "#2F6DB9",
                         "#267255")

claudinLowColor <- "yellow"

subtypeColors <- c(subtypeColorsSingle[1], subtypeColorsSingle[1],
                   subtypeColorsSingle[2], subtypeColorsSingle[2],
                   subtypeColorsSingle[3], subtypeColorsSingle[3],
                   subtypeColorsSingle[4], subtypeColorsSingle[4],
                   subtypeColorsSingle[5], subtypeColorsSingle[5])


##################################################
## Distribution of tumors by subtype (all)

# Get subtype distribution
tumorDistribution <- table(patientData$PAM50toCL)

# Format subtype distribution correctly
tumorDistributionCounts <- data.frame(Subtype = c(tumorDistribution[["Basal"]],
                                                  tumorDistribution[["BasalClaudinLow"]],
                                                  tumorDistribution[["Her2"]],
                                                  tumorDistribution[["Her2ClaudinLow"]],
                                                  tumorDistribution[["LumA"]],
                                                  tumorDistribution[["LumAClaudinLow"]],
                                                  tumorDistribution[["LumB"]],
                                                  tumorDistribution[["LumBClaudinLow"]],
                                                  tumorDistribution[["Normal"]],
                                                  tumorDistribution[["NormalClaudinLow"]]), 
                                      row.names = c("Basal",
                                                    "BasalClaudinLow",
                                                    "Her2",
                                                    "Her2ClaudinLow",
                                                    "LumA",
                                                    "LumAClaudinLow",
                                                    "LumB",
                                                    "LumBClaudinLow",
                                                    "Normal",
                                                    "NormalClaudinLow"))

tumorDistributionCounts <- as.matrix(tumorDistributionCounts)

# Matrix for only claudin-low tumors
claudinCounts <- as.matrix(tumorDistributionCounts[grep(x = rownames(tumorDistributionCounts), pattern = "ClaudinLow"), ], dimnames(x = "Subtypes"))

# Get n
pts <- nrow(patientData)

# Add 2 because two claudin-low tumors have unknown PAM50 subtype
clds <- sum(claudinCounts[, 1]) + 2

# Make plot
pdf(file = "./Output/Figures/Oslo2_SubtypeDistributions.pdf", height = 5, width = 8)
par(mfrow = c(2,1), xpd = NA)

# All tumors top bar
barplot(tumorDistributionCounts,
        col = subtypeColors,
        horiz = TRUE,
        xlab = "",
        axisnames = FALSE,
        axes = FALSE,
        legend = TRUE,
        args.legend = list(x = "top",
                           bty = "n",
                           ncol = 5,
                           inset = c(0, -1),
                           cex = 0.75,
                           xpd = NA))
title(ylab=paste("Subtype distribution \n (n = ", pts, ")", sep = ""), mgp = c(1,1,0), cex.lab = 1.2)      

barplot(tumorDistributionCounts,
        add = TRUE,
        horiz = TRUE,
        xlab = "",
        axisnames = FALSE,
        axes = FALSE,
        legend = TRUE,
        args.legend = list(x = "top",
                           bty = "n",
                           ncol = 5,
                           inset = c(0, -1),
                           cex = 0.75,
                           xpd = NA),
        angle = 60,
        density = c(0, 20),
        col = claudinLowColor)

# Claudin-low bottom bar
barplot(claudinCounts,
        col = subtypeColorsSingle,
        horiz = TRUE,
        xlab = "",
        axisnames = FALSE,
        axes = FALSE)
title(ylab=paste("Claudin-low \n (n = ", clds, ")", sep = ""), mgp = c(1,1,0), cex.lab = 1.2)    

barplot(claudinCounts,
        add = TRUE,
        horiz = TRUE,
        xlab = "",
        axisnames = FALSE,
        axes = FALSE,
        angle = c(60, 120),
        density = 20,
        col = claudinLowColor)
dev.off()


##################################################
## Distribution of subtypes in core and other claudin-low tumors

# Find distribution of core claudin-low tumors
coreCLDistribution <- table(patientData$pam50[patientData$CoreCLOtherCL == "CoreClaudinLow"])
coreCLDistributionCounts <- data.frame(CoreClaudinLow = c(coreCLDistribution[["Basal"]],
                                                          0, #coreCLDistribution[["Her2"]], # No HER2E CoreCL
                                                          coreCLDistribution[["LumA"]],
                                                          coreCLDistribution[["LumB"]], 
                                                          coreCLDistribution[["Normal"]]), 
                                       row.names = c("Basal",
                                                     "Her2",
                                                     "LumA",
                                                     "LumB",
                                                     "Normal"))

# Find distribution of other claudin-low tumors
otherCLDistribution <- table(patientData$pam50[patientData$CoreCLOtherCL == "OtherClaudinLow"])
otherCLDistributionCounts <- data.frame(OtherClaudinLow = c(otherCLDistribution[["Basal"]],
                                                            otherCLDistribution[["Her2"]],
                                                            otherCLDistribution[["LumA"]],
                                                            otherCLDistribution[["LumB"]],
                                                            0), #otherCLDistribution[["Normal"]]), # No Normal-like OtherCL
                                        row.names = c("Basal",
                                                      "Her2",
                                                      "LumA",
                                                      "LumB",
                                                      "Normal"))

# Gather the above two in a single matrix
clDistribution <- as.matrix(cbind(coreCLDistributionCounts, otherCLDistributionCounts))

# Make plot
pdf("./Output/Figures/Oslo2_CoreClaudinLow_Distribution.pdf", width = 4, height = 3)
par(mar=c(5,8,4,2))
barplot(clDistribution,
        col = subtypeColorsSingle,
        horiz = TRUE,
        xlab = "Number of tumors",
        xlim = c(0, 30),
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
dev.off()

print("Finished generating plots for the Oslo2 cohort")
