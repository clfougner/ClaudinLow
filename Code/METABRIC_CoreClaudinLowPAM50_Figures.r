##################################################
## Figures for CoreCL tumors in the context of their PAM50 subtype (METABRIC cohort)
##
## Requires the output from ./Code/METABRIC_CoreClaudinLow.r (./Output/METABRIC_patientData_CoreClaudinLow.txt). Make sure that the working directory is set to /path/to/ClaudinLow/
##################################################

# Load the required packages
library(survival)
library(survminer)
library(ggsci)

print("METABRIC: Generating figures for CoreCL tumors in the context of their intrinsic subtype")

##################################################
## Get relevant data

# Read patient data output from ./Code/METABRIC_OutputAllDataPerSample.r
patientData <- read.table(file = "./Output/METABRIC_patientData_CoreClaudinLow.txt", header = TRUE, sep = "\t")

# Set only CoreCL as claudin-low, consider OtherCL as non-claudin-low
patientData$PAM50ClaudinLow[patientData$PAM50toCoreCL == "Basal"] <- "Basal"
patientData$PAM50ClaudinLow[patientData$PAM50toCoreCL == "BasalOtherClaudinLow"] <- "Basal"
patientData$PAM50ClaudinLow[patientData$PAM50toCoreCL == "BasalCoreClaudinLow"] <- "BasalClaudinLow"

patientData$PAM50ClaudinLow[patientData$PAM50toCoreCL == "Her2"] <- "Her2"
patientData$PAM50ClaudinLow[patientData$PAM50toCoreCL == "Her2OtherClaudinLow"] <- "Her2"
patientData$PAM50ClaudinLow[patientData$PAM50toCoreCL == "Her2CoreClaudinLow"] <- "Her2ClaudinLow"

patientData$PAM50ClaudinLow[patientData$PAM50toCoreCL == "LumA"] <- "LumA"
patientData$PAM50ClaudinLow[patientData$PAM50toCoreCL == "LumAOtherClaudinLow"] <- "LumA"
patientData$PAM50ClaudinLow[patientData$PAM50toCoreCL == "LumACoreClaudinLow"] <- "LumAClaudinLow"

patientData$PAM50ClaudinLow[patientData$PAM50toCoreCL == "LumB"] <- "LumB"
patientData$PAM50ClaudinLow[patientData$PAM50toCoreCL == "LumBOtherClaudinLow"] <- "LumB"
patientData$PAM50ClaudinLow[patientData$PAM50toCoreCL == "LumBCoreClaudinLow"] <- "LumBClaudinLow"

patientData$PAM50ClaudinLow[patientData$PAM50toCoreCL == "Normal"] <- "Normal"
patientData$PAM50ClaudinLow[patientData$PAM50toCoreCL == "NormalOtherClaudinLow"] <- "Normal"
patientData$PAM50ClaudinLow[patientData$PAM50toCoreCL == "NormalCoreClaudinLow"] <- "NormalClaudinLow"

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
## Distribution of tumors by subtype

# Tumor number of tumors from each subtype
tumorDistribution <- table(patientData$PAM50ClaudinLow)

# Format the subtype distribution
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
clds <- sum(claudinCounts[, 1])

# Make plot
pdf(file = "./Output/Figures/METABRIC_CoreClaudinLow_SubtypeDistribution.pdf", height = 5, width = 8, pointsize = 10)
  par(mfrow = c(2,1), xpd = NA)
  
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
  
  # Claudin-low tumors only
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
## Remove HER2-enriched and LumB due to low sample numbers

# Remove Her2 and LumB
patientData <- patientData[-grep(x = patientData$PAM50ClaudinLow, pattern = "Her2"), ]
patientData <- patientData[-grep(x = patientData$PAM50ClaudinLow, pattern = "LumB"), ]
patientData <- droplevels(patientData)

# Subtype colors for only basal-like, LumA, normal-like
subtypeColors <- c(subtypeColorsSingle[1], subtypeColorsSingle[1],
                   subtypeColorsSingle[3], subtypeColorsSingle[3],
                   subtypeColorsSingle[5], subtypeColorsSingle[5])

subtypeColorsSingle <- c(subtypeColorsSingle[1],
                         subtypeColorsSingle[3],
                         subtypeColorsSingle[5])


##################################################
## Get sample numbers and proportions for basal-like, LumA and normal-like tumors

tumorDistribution <- table(patientData$PAM50ClaudinLow)

# Format subtype distributions
tumorDistributionCounts <- data.frame(Subtype = c(tumorDistribution[["Basal"]],
                                                  tumorDistribution[["BasalClaudinLow"]],
                                                  tumorDistribution[["LumA"]],
                                                  tumorDistribution[["LumAClaudinLow"]],
                                                  tumorDistribution[["Normal"]],
                                                  tumorDistribution[["NormalClaudinLow"]]), 
                                      row.names = c("Basal",
                                                    "BasalClaudinLow",
                                                    "LumA",
                                                    "LumAClaudinLow",
                                                    "Normal",
                                                    "NormalClaudinLow"))

tumorDistributionCounts <- as.matrix(tumorDistributionCounts)


##################################################
## Distribution of ER
allProps <- c()
for(subtype in rownames(tumorDistributionCounts)){
  totalPatients <- sum(table(patientData$ER_IHC[patientData$PAM50ClaudinLow == subtype]))
  pos <- table(patientData$ER_IHC[patientData$PAM50ClaudinLow == subtype])[["Positve"]]
  proportion <- pos/totalPatients
  allProps <- c(allProps, proportion)
}
names(allProps) <- rownames(tumorDistributionCounts)

pdf(file = "./Output/Figures/METABRIC_CoreClaudinLow_ER_distribution.pdf", height = 3, width = 2.5, pointsize = 8)
  par(mar=c(10,4,4,2))
  barplot(allProps, 
          col = subtypeColors,
          ylab = "Proportion ER+",
          main = "Estrogen receptor",
          las = 2,
          space=c(0.2, 0.07, 0.2, 0.07))
  
  barplot(allProps,
          add = TRUE,
          las = 2,
          space=c(0.2, 0.07, 0.2, 0.07),
          angle = 60,
          density = c(0, 20),
          col = claudinLowColor)
dev.off()


##################################################
## Distribution of TP53 mutations
allProps <- c()
for(subtype in rownames(tumorDistributionCounts)){
  totalPatients <- sum(table(patientData$TP53_Status[patientData$PAM50ClaudinLow == subtype]))
  pos <- table(patientData$TP53_Status[patientData$PAM50ClaudinLow == subtype])[["Mut"]]
  proportion <- pos/totalPatients
  allProps <- c(allProps, proportion)
}
names(allProps) <- rownames(tumorDistributionCounts)

pdf(file = "./Output/Figures/METABRIC_CoreClaudinLow_TP53_distribution.pdf", height = 3, width = 2.5, pointsize = 8)
  par(mar=c(10,4,4,2))
  barplot(allProps, 
          col = subtypeColors,
          ylab = "Proportion with mutation",
          main = "TP53 mutation",
          ylim = c(0, 1),
          las = 2,
          space=c(0.2, 0.07, 0.2, 0.07))
  
  barplot(allProps,
          add = TRUE,
          las = 2,
          space=c(0.2, 0.07, 0.2, 0.07),
          ylim = c(0, 1),
          angle = 60,
          density = c(0, 20),
          col = claudinLowColor)
dev.off()


##################################################
## Distribution of PIK3CA mutations
allProps <- c()
for(subtype in rownames(tumorDistributionCounts)){
  totalPatients <- sum(table(patientData$PIK3CA_Status[patientData$PAM50ClaudinLow == subtype]))
  pos <- table(patientData$PIK3CA_Status[patientData$PAM50ClaudinLow == subtype])[["Mut"]]
  proportion <- pos/totalPatients
  allProps <- c(allProps, proportion)
}
names(allProps) <- rownames(tumorDistributionCounts)

pdf(file = "./Output/Figures/METABRIC_CoreClaudinLow_PIK3CA_distribution.pdf", height = 3, width = 2.5, pointsize = 8)
  par(mar=c(10,4,4,2))
  barplot(allProps, 
          col = subtypeColors,
          ylab = "Proportion with mutation",
          main = "PIK3CA mutation",
          ylim = c(0, 0.6),
          las = 2,
          space=c(0.2, 0.07, 0.2, 0.07))
  
  barplot(allProps,
          add = TRUE,
          las = 2,
          space=c(0.2, 0.07, 0.2, 0.07),
          ylim = c(0, 0.6),
          angle = 60,
          density = c(0, 20),
          col = claudinLowColor)
dev.off()


##################################################
## Distribution of MYC gain
allProps <- c()
for(subtype in rownames(tumorDistributionCounts)){
  totalPatients <- sum(table(patientData$MYC_SIMPLE[patientData$PAM50ClaudinLow == subtype]))
  pos <- table(patientData$MYC_SIMPLE[patientData$PAM50ClaudinLow == subtype])[["AMP"]]
  proportion <- pos/totalPatients
  allProps <- c(allProps, proportion)
}
names(allProps) <- rownames(tumorDistributionCounts)

pdf(file = "./Output/Figures/METABRIC_CoreClaudinLow_MYC_distribution.pdf", height = 3, width = 2.5, pointsize = 8)
  par(mar=c(10,4,4,2))
  barplot(allProps, 
          col = subtypeColors,
          ylab = "Proportion with gain",
          main = "MYC gain",
          ylim = c(0, 0.8),
          las = 2,
          space=c(0.2, 0.07, 0.2, 0.07))
  
  barplot(allProps,
          add = TRUE,
          las = 2,
          space=c(0.2, 0.07, 0.2, 0.07),
          ylim = c(0, 0.8),
          angle = 60,
          density = c(0, 20),
          col = claudinLowColor)
dev.off()


##################################################
## Distribution of MDM4 gain
allProps <- c()
for(subtype in rownames(tumorDistributionCounts)){
  totalPatients <- sum(table(patientData$MDM4_SIMPLE[patientData$PAM50ClaudinLow == subtype]))
  pos <- table(patientData$MDM4_SIMPLE[patientData$PAM50ClaudinLow == subtype])[["AMP"]]
  proportion <- pos/totalPatients
  allProps <- c(allProps, proportion)
}
names(allProps) <- rownames(tumorDistributionCounts)

pdf(file = "./Output/Figures/METABRIC_CoreClaudinLow_MDM4_distribution.pdf", height = 3, width = 2.5, pointsize = 8)
  par(mar=c(10,4,4,2))
  barplot(allProps, 
          col = subtypeColors,
          ylab = "Proportion with gain",
          main = "MDM4 gain", 
          ylim = c(0, 0.8),
          las = 2,
          space=c(0.2, 0.07, 0.2, 0.07))
  
  barplot(allProps,
          add = TRUE,
          las = 2,
          space=c(0.2, 0.07, 0.2, 0.07),
          ylim = c(0, 0.8),
          angle = 60,
          density = c(0, 20),
          col = claudinLowColor)
dev.off()

##################################################
## Distribution of IntClust4ER+
allProps <- c()
for(subtype in rownames(tumorDistributionCounts)){
  totalPatients <- sum(table(patientData$INTCLUST[patientData$PAM50ClaudinLow == subtype]))
  pos <- table(patientData$INTCLUST[patientData$PAM50ClaudinLow == subtype])[["4ER+"]]
  proportion <- pos/totalPatients
  allProps <- c(allProps, proportion)
}
names(allProps) <- rownames(tumorDistributionCounts)

pdf(file = "./Output/Figures/METABRIC_CoreClaudinLow_IC4ER_distribution.pdf", height = 3, width = 5, pointsize = 8)
par(mfrow = c(1,2), xpd = NA, mar = c(10,4,4,2))

barplot(allProps, 
        col = subtypeColors,
        ylab = "Proportion",
        main = "IntClust4ER+",
        ylim = c(0, 1),
        las = 2,
        space=c(0.2, 0.07, 0.2, 0.07))

barplot(allProps,
        add = TRUE,
        las = 2,
        space=c(0.2, 0.07, 0.2, 0.07),
        angle = 60,
        density = c(0, 20),
        col = claudinLowColor)

##################################################
## Distribution of IntClust4ER-
allProps <- c()
for(subtype in rownames(tumorDistributionCounts)){
  totalPatients <- sum(table(patientData$INTCLUST[patientData$PAM50ClaudinLow == subtype]))
  pos <- table(patientData$INTCLUST[patientData$PAM50ClaudinLow == subtype])[["4ER-"]]
  proportion <- pos/totalPatients
  allProps <- c(allProps, proportion)
}
names(allProps) <- rownames(tumorDistributionCounts)

barplot(allProps, 
        col = subtypeColors,
        ylab = "Proportion",
        main = "IntClust4ER-",
        ylim = c(0, 1),
        las = 2,
        space=c(0.2, 0.07, 0.2, 0.07))

barplot(allProps,
        add = TRUE,
        las = 2,
        space=c(0.2, 0.07, 0.2, 0.07),
        angle = 60,
        density = c(0, 20),
        col = claudinLowColor)
dev.off()

##################################################
## Distribution of IntClust4
allProps <- c()
for(subtype in rownames(tumorDistributionCounts)){
  totalPatients <- sum(table(patientData$INTCLUST[patientData$PAM50ClaudinLow == subtype]))
  pos <- table(patientData$INTCLUST[patientData$PAM50ClaudinLow == subtype])[["4ER+"]] + table(patientData$INTCLUST[patientData$PAM50ClaudinLow == subtype])[["4ER-"]]
    proportion <- pos/totalPatients
  allProps <- c(allProps, proportion)
}
names(allProps) <- rownames(tumorDistributionCounts)

pdf(file = "./Output/Figures/METABRIC_CoreClaudinLow_IC4_distribution.pdf", height = 3, width = 2.5, pointsize = 8)
par(mar=c(10,4,4,2))
barplot(allProps, 
        col = subtypeColors,
        ylab = "Proportion",
        main = "IntClust4",
        ylim = c(0, 1),
        las = 2,
        space=c(0.2, 0.07, 0.2, 0.07))

barplot(allProps,
        add = TRUE,
        las = 2,
        space=c(0.2, 0.07, 0.2, 0.07),
        angle = 60,
        density = c(0, 20),
        col = claudinLowColor)
dev.off()

##################################################
## Distribution of GII
pdf("./Output/Figures/METABRIC_CoreClaudinLow_GII_distribution.pdf", height = 3, width = 2.5, pointsize = 8)
  par(mar=c(10,4,4,2))
  giiBP <- boxplot(GII_PloidyCorrected ~ PAM50ClaudinLow,
                   data = patientData,
                   ylab = "GII",
                   main = "Genomic instability index",
                   ylim = c(0, 1),
                   xlab = "",
                   at = c(1, 1.85, 3, 3.85, 5, 5.85),
                   las = 2,
                   frame = FALSE,
                   outline = FALSE,
                   col = subtypeColors)
  
  rect(xleft = c(1.45, 3.45, 5.45),
       ybottom = giiBP$stats[2, c(2, 4, 6)],
       xright = c(2.25, 4.25, 6.25),
       ytop = giiBP$stats[4, c(2, 4, 6)],
       angle = 60,
       density = 20,
       col = claudinLowColor,
       border = "black")
  
  rect(xleft = c(1.45, 3.45, 5.45),
       ybottom = giiBP$stats[3, c(2, 4, 6)],
       xright = c(2.25, 4.25, 6.25),
       ytop = giiBP$stats[3, c(2, 4, 6)],
       lwd = 2)

dev.off()


##################################################
## Distribution of mutation counts
pdf("./Output/Figures/METABRIC_CoreClaudinLow_Mutation_distribution.pdf", height = 3, width = 2.5, pointsize = 8)
  par(mar = c(10,4,4,2))
  mutBP <- boxplot(log2(MutationCounts + 1) ~ PAM50ClaudinLow,
                   data = patientData,
                   ylab = "Log2(Mutations + 1)",
                   main = "Mutations", 
                   xlab = "",
                   ylim = c(0, 6),
                   at = c(1, 1.85, 3, 3.85, 5, 5.85),
                   las = 2,
                   frame = FALSE,
                   outline = FALSE,
                   col = subtypeColors)
  
  rect(xleft = c(1.45, 3.45, 5.45),
       ybottom = mutBP$stats[2, c(2, 4, 6)],
       xright = c(2.25, 4.25, 6.25),
       ytop = mutBP$stats[4, c(2, 4, 6)],
       angle = 60,
       density = 20,
       col = claudinLowColor,
       border = "black")
  
  rect(xleft = c(1.45, 3.45, 5.45),
       ybottom = mutBP$stats[3, c(2, 4, 6)],
       xright = c(2.25, 4.25, 6.25),
       ytop = mutBP$stats[3, c(2, 4, 6)],
       lwd = 2)
  
dev.off()

##################################################
## Distribution of MKI67
pdf("./Output/Figures/METABRIC_CoreClaudinLow_Ki67_distribution.pdf", height = 3, width = 2.5, pointsize = 8)
par(mar = c(10,4,4,2))
mki67BP <- boxplot(MKI67 ~ PAM50ClaudinLow,
                 data = patientData,
                 ylab = "MKI67 expression",
                 main = "MKI67", 
                 xlab = "",
                 at = c(1, 1.85, 3, 3.85, 5, 5.85),
                 las = 2,
                 frame = FALSE,
                 outline = FALSE,
                 col = subtypeColors)

rect(xleft = c(1.45, 3.45, 5.45),
     ybottom = mki67BP$stats[2, c(2, 4, 6)],
     xright = c(2.25, 4.25, 6.25),
     ytop = mki67BP$stats[4, c(2, 4, 6)],
     angle = 60,
     density = 20,
     col = claudinLowColor,
     border = "black")

rect(xleft = c(1.45, 3.45, 5.45),
     ybottom = mki67BP$stats[3, c(2, 4, 6)],
     xright = c(2.25, 4.25, 6.25),
     ytop = mki67BP$stats[3, c(2, 4, 6)],
     lwd = 2)

dev.off()


##################################################
## Distribution of ImmuneScore and StromalScore
pdf("./Output/Figures/METABRIC_CoreClaudinLow_ImmuneScore_distribution.pdf", height = 3, width = 2.5, pointsize = 8)

  ## ImmuneScore
  par(mar = c(10,4,4,2))
  immuneBP <- boxplot(ImmuneScore ~ PAM50ClaudinLow,
                   data = patientData,
                   ylab = "Immune score",
                   main = "Immune score",
                   xlab = "",
                   ylim = c(-1000, 4000),
                   at = c(1, 1.85, 3, 3.85, 5, 5.85),
                   las = 2,
                   frame = FALSE,
                   outline = FALSE,
                   col = subtypeColors)
  
  rect(xleft = c(1.45, 3.45, 5.45),
       ybottom = immuneBP$stats[2, c(2, 4, 6)],
       xright = c(2.25, 4.25, 6.25),
       ytop = immuneBP$stats[4, c(2, 4, 6)],
       angle = 60,
       density = 20,
       col = claudinLowColor,
       border = "black")
  
  rect(xleft = c(1.45, 3.45, 5.45),
       ybottom = immuneBP$stats[3, c(2, 4, 6)],
       xright = c(2.25, 4.25, 6.25),
       ytop = immuneBP$stats[3, c(2, 4, 6)],
       lwd = 2)
  
dev.off()
  
  ## StromalScore
pdf("./Output/Figures/METABRIC_CoreClaudinLow_StromalScore_distribution.pdf", height = 3, width = 2.5, pointsize = 8)
  par(mar = c(10,4,4,2))
  stromalBP <- boxplot(StromalScore ~ PAM50ClaudinLow,
                   data = patientData,
                   ylab = "Stromal score",
                   main = "Stromal score",
                   xlab = "",
                   ylim = c(-1000, 2000),
                   at = c(1, 1.85, 3, 3.85, 5, 5.85),
                   las = 2,
                   frame = FALSE,
                   outline = FALSE,
                   col = subtypeColors)
  
  rect(xleft = c(1.45, 3.45, 5.45),
       ybottom = stromalBP$stats[2, c(2, 4, 6)],
       xright = c(2.25, 4.25, 6.25),
       ytop = stromalBP$stats[4, c(2, 4, 6)],
       angle = 60,
       density = 20,
       col = claudinLowColor,
       border = "black")
  
  rect(xleft = c(1.45, 3.45, 5.45),
       ybottom = stromalBP$stats[3, c(2, 4, 6)],
       xright = c(2.25, 4.25, 6.25),
       ytop = stromalBP$stats[3, c(2, 4, 6)],
       lwd = 2)
dev.off()

print("METABRIC: Finished generating figures for CoreCL tumors in the context of their intrinsic subtype")