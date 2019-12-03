##################################################
## Figures for claudin-low tumors in the context of their intrinsic subtype (METABRIC cohort)
##
## Requires the output from ./Code/METABRIC_patientData.r (./Output/METABRIC_patientData.txt). Ensure that the working directory is set to /path/to/ClaudinLow/
##################################################

print("METABRIC: generating figures for claudin-low tumors stratified by intrinsic subtype")
print("Ignore warnings generated here")

# Load the required packages
library(survival)
library(survminer)
library(ggsci)

##################################################
## Read input data data

# Read patient data output from ./Code/METABRIC_patientData.r
patientData <- read.table(file = "./Output/METABRIC_patientData.txt", header = TRUE, sep = "\t")

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
pdf(file = "./Output/Figures/METABRIC_SubtypeDistribution.pdf", height = 5, width = 8, pointsize = 10)
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
          density = c(0, 10),
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
          density = 10,
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

pdf(file = "./Output/Figures/METABRIC_ER_distribution.pdf", height = 3, width = 2.5, pointsize = 8)
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

pdf(file = "./Output/Figures/METABRIC_TP53_distribution.pdf", height = 3, width = 2.5, pointsize = 8)
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

pdf(file = "./Output/Figures/METABRIC_PIK3CA_distribution.pdf", height = 3, width = 2.5, pointsize = 8)
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

pdf(file = "./Output/Figures/METABRIC_MYC_distribution.pdf", height = 3, width = 2.5, pointsize = 8)
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

pdf(file = "./Output/Figures/METABRIC_MDM4_distribution.pdf", height = 3, width = 2.5, pointsize = 8)
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

pdf(file = "./Output/Figures/METABRIC_IC4ER_distribution.pdf", height = 3, width = 5, pointsize = 8)
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
        ylim = c(0, 0.4),
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

pdf(file = "./Output/Figures/METABRIC_IC4_distribution.pdf", height = 3, width = 2.5, pointsize = 8)
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
pdf("./Output/Figures/METABRIC_GII_distribution.pdf", height = 3, width = 2.5, pointsize = 8)
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
pdf("./Output/Figures/METABRIC_Mutation_distribution.pdf", height = 3, width = 2.5, pointsize = 8)
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
pdf("./Output/Figures/METABRIC_MKI67_distribution.pdf", height = 3, width = 2.5, pointsize = 8)
par(mar=c(10,4,4,2))
mki67BP <- boxplot(MKI67 ~ PAM50ClaudinLow,
                 data = patientData,
                 ylab = "MKI67 gene expression",
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
## Distribution of Age at diagnosis
pdf("./Output/Figures/METABRIC_AgeAtDiagnosis_distribution.pdf", height = 3, width = 2.5, pointsize = 8)
par(mar=c(10,4,4,2))
ageBP <- boxplot(AGE_AT_DIAGNOSIS ~ PAM50ClaudinLow,
                 data = patientData,
                 ylab = "Age at diagnosis",
                 main = "Age at diagnosis",
                 ylim = c(min(pretty(patientData$AGE_AT_DIAGNOSIS)), max(pretty(patientData$AGE_AT_DIAGNOSIS))),
                 xlab = "",
                 at = c(1, 1.85, 3, 3.85, 5, 5.85),
                 las = 2,
                 frame = FALSE,
                 outline = FALSE,
                 col = subtypeColors)

rect(xleft = c(1.45, 3.45, 5.45),
     ybottom = ageBP$stats[2, c(2, 4, 6)],
     xright = c(2.25, 4.25, 6.25),
     ytop = ageBP$stats[4, c(2, 4, 6)],
     angle = 60,
     density = 20,
     col = claudinLowColor,
     border = "black")

rect(xleft = c(1.45, 3.45, 5.45),
     ybottom = ageBP$stats[3, c(2, 4, 6)],
     xright = c(2.25, 4.25, 6.25),
     ytop = ageBP$stats[3, c(2, 4, 6)],
     lwd = 2)

dev.off()


##################################################
## Distribution of ImmuneScore and StromalScore
pdf("./Output/Figures/METABRIC_Immune_Stromal_Score_distribution.pdf", height = 3, width = 5, pointsize = 8)

  ## ImmuneScore
  par(mfrow = c(1,2), xpd = NA, mar = c(10,4,4,2))
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

  ## StromalScore
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


##################################################
## Kaplan-Meier plots

# Create a data frame for the relevant data
survivalDF_root <- data.frame(sample = patientData$PATIENT_ID,
                              subtype = patientData$PAM50ClaudinLow,
                              survival = patientData$OS_MONTHS,
                              status = patientData$OS_STATUS,
                              DiseaseSpecificSurvival = patientData$VITAL_STATUS)


################################################
## All claudin lows
survivalDF <- survivalDF_root[grep(x = survivalDF_root$subtype, pattern = "ClaudinLow"), ]
survivalDF$SurvObj <- Surv(survivalDF$survival, event = survivalDF$DiseaseSpecificSurvival== "Died of Disease")
survivalByType <- survfit(SurvObj ~ subtype, data = survivalDF)

survivalPlotClaudinLow <- ggsurvplot(survivalByType,
                                     data = survivalDF,
                                     palette = c(subtypeColorsSingle[1],
                                                 subtypeColorsSingle[2],
                                                 subtypeColorsSingle[3]),
                                     pval = TRUE,
                                     legend = "none",
                                     censor = FALSE,
                                     xlim = c(0, 200),
                                     break.x.by = 50,
                                     xlab = "Months",
                                     title = "All claudin-low")

################################################
## Basal-like + ClaudinLow
survivalDF <- survivalDF_root[grep(x = survivalDF_root$subtype, pattern = "Basal"), ]
survivalDF$SurvObj <- Surv(survivalDF$survival, event = survivalDF$DiseaseSpecificSurvival== "Died of Disease")
survivalByType <- survfit(SurvObj ~ subtype, data = survivalDF)

survivalPlotBasal <- ggsurvplot(survivalByType,
                                data = survivalDF,
                                palette = c(subtypeColorsSingle[1],
                                            claudinLowColor),
                                pval = TRUE,
                                legend = "none",
                                censor = FALSE,
                                xlim = c(0, 200),
                                break.x.by = 50,
                                xlab = "Months",
                                ylab = "",
                                title = "Basal-like")

################################################
## LumA + ClaudinLow
survivalDF <- survivalDF_root[grep(x = survivalDF_root$subtype, pattern = "LumA"), ]
survivalDF$SurvObj <- Surv(survivalDF$survival, event = survivalDF$DiseaseSpecificSurvival== "Died of Disease")
survivalByType <- survfit(SurvObj ~ subtype, data = survivalDF)

survivalPlotLumA <- ggsurvplot(survivalByType,
                               data = survivalDF,
                               palette = c(subtypeColorsSingle[2],
                                           claudinLowColor),
                               pval = TRUE,
                               legend = "none",
                               censor = FALSE,
                               xlim = c(0, 200),
                               break.x.by = 50,
                               xlab = "Months",
                               ylab = "",
                               title = "LumA")

################################################
## Normal-like + ClaudinLow
survivalDF <- survivalDF_root[grep(x = survivalDF_root$subtype, pattern = "Normal"), ]
survivalDF$SurvObj <- Surv(survivalDF$survival, event = survivalDF$DiseaseSpecificSurvival== "Died of Disease")
survivalByType <- survfit(SurvObj ~ subtype, data = survivalDF)

survivalPlotNormal <- ggsurvplot(survivalByType,
                                 data = survivalDF,
                                 palette = c(subtypeColorsSingle[3],
                                             claudinLowColor),
                                 pval = TRUE,
                                 legend = "none",
                                 censor = FALSE,
                                 xlim = c(0, 200),
                                 break.x.by = 50,
                                 xlab = "Months",
                                 ylab = "",
                                 title = "Normal-like")


## Write Kaplan-Meier plots to file
pdf(file = "./Output/Figures/METABRIC_KaplanMeier_ClaudinLow.pdf", height = 2.5, width = 10, onefile=FALSE)
  arrange_ggsurvplots(list(survivalPlotClaudinLow,
                           survivalPlotBasal,
                           survivalPlotLumA,
                           survivalPlotNormal),
                      nrow = 1, ncol = 4)
dev.off()


################################################
## Kaplan-Meier plot for all subtypes, claudin-low tumors not stratified by subtype
survivalDF <- survivalDF_root
survivalDF$subtype <- as.character(survivalDF$subtype)
survivalDF[grep(x = survivalDF$subtype, pattern = "ClaudinLow"), "subtype"] <- "ClaudinLow"
survivalDF$SurvObj <- Surv(survivalDF$survival, event = survivalDF$DiseaseSpecificSurvival== "Died of Disease")
survivalByType <- survfit(SurvObj ~ subtype, data = survivalDF)

survivalPlotClaudinLow <- ggsurvplot(survivalByType,
                                     data = survivalDF,
                                     palette = c(subtypeColorsSingle[1],
                                                 claudinLowColor,
                                                 subtypeColorsSingle[2],
                                                 subtypeColorsSingle[3]),
                                     pval = TRUE,
                                     legend = "right",
                                     censor = FALSE,
                                     xlim = c(0, 200),
                                     break.x.by = 50,
                                     xlab = "Months",
                                     title = "Claudin-low not stratified")

pdf("./Output/Figures/METABRIC_Survival_ClaudinLowNotStratified.pdf", width = 6, height = 4, onefile = FALSE)
  print(survivalPlotClaudinLow)
dev.off()


################################################
## GII vs. distance to CL-centroid
patientData <- read.table(file = "./Output/METABRIC_patientData.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Linear regression of distance to claudin-low centroid vs. GII
fit <- lm(patientData$CL_distance ~ patientData$GII_PloidyCorrected)

# Function to get p-value (source: https://stackoverflow.com/questions/5587676/pull-out-p-values-and-r-squared-from-a-linear-regression)
lmp <- function (modelobject) {
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1], f[2], f[3], lower.tail = F)
  attributes(p) <- NULL
  return(p)
}

# Plot
patientData$INTCLUST <- factor(patientData$INTCLUST, levels = levels(as.factor(patientData$INTCLUST))[c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 2)])
pdf(file = "./Output/Figures/METABRIC_CLDist_vs_GII.pdf", width = 7, height = 6, useDingbats = FALSE)
  palette(pal_ucscgb(alpha = 1)(11))
  plot(patientData$CL_distance ~ patientData$GII_PloidyCorrected,
       col = patientData$INTCLUST,
       pch = 19,
       ylim = c(10, 70),
       cex=1,
       bty="n",
       xlab = "Genomic instability index",
       main = "Correlation between GII and distance \n to claudin-low centroid",
       ylab = "Distance to claudin-low centroid")
  
  abline(fit, cex = 10)
  text(x = 0.9, y = max(patientData$CL_distance), labels = paste("R^2 =", round(summary(fit)$r.squared, 2)))
  # The following p value has been checked manually, hard coded for formatting
  text(x = 0.9, y = 0.95*max(patientData$CL_distance), labels = "P < 0.001")
  
  legend("bottomleft", inset=.02, legend = levels(patientData$INTCLUST), fill = 1:length(levels(patientData$INTCLUST)), horiz=FALSE, title = "IntClust", ncol = 6)
  palette("default")
dev.off()



################################################
## Estimate score vs. distance to CL-centroid

subtypeColorsSingle <- c("#DC362A",
                         "#ED1E78",
                         "#333A8E",
                         "#2F6DB9",
                         "#267255")

claudinLowColor <- "yellow"

patientData$PAM50 <- factor(patientData$PAM50, levels = levels(as.factor(patientData$PAM50)))

## CLDist ~ Stromal
pdf(file = "./Output/Figures/METABRIC_CLDist_vs_ESTIMATE_Stromal_9CLFill.pdf",
    width = 3.5,
    height = 3.5,
    pointsize = 7, useDingbats = FALSE)

  # Linear regression of distance to claudin-low centroid vs. GII
  fit <- lm(patientData$CL_distance ~ patientData$StromalScore)
  palette(subtypeColorsSingle)
  plot(patientData$CL_distance[patientData$ClaudinLow != "ClaudinLow"] ~ patientData$StromalScore[patientData$ClaudinLow != "ClaudinLow"],
       bg = patientData$PAM50[patientData$ClaudinLow != "ClaudinLow"],
       lwd = 0.25,
       pch = 21,
       ylim = c(10, 70),
       xlim = c(min(patientData$StromalScore), max(patientData$StromalScore)),
       font.main = 1,
       cex = 1,
       bty="n",
       xlab = "Stromal score",
       main = "Correlation between stromal score and distance \n to claudin-low centroid",
       ylab = "Distance to claudin-low centroid")
  
  for(i in (1:nrow(patientData))[patientData$ClaudinLow == "ClaudinLow"]){
    
    points(patientData$CL_distance[i] ~ patientData$StromalScore[i],
           bg = patientData$PAM50[i],
           lwd = 0.25,
           pch = 21,
           cex = 1,
           bty="n")
    
    points(patientData$CL_distance[i] ~ patientData$StromalScore[i],
           col = claudinLowColor,
           pch = 21,
           cex = 1,
           lwd = 1.5)
    
    points(patientData$CL_distance[i] ~ patientData$StromalScore[i],
           col = "black",
           pch = 21,
           cex = 1.3,
           lwd = 0.25)
    
  }
  
  abline(fit, cex = 10)
  text(x = 1500, y = max(patientData$CL_distance), labels = paste("R^2 =", round(summary(fit)$r.squared, 2)))
  # The following p value has been checked manually, hard coded for formatting
  text(x = 1500, y = 0.95*max(patientData$CL_distance), labels = "P < 0.001")
dev.off()


## OtherDist ~ Immune
pdf(file = "./Output/Figures/METABRIC_OtherDist_vs_ESTIMATE_Stromal_9CLFill.pdf",
    width = 3.5,
    height = 3.5,
    pointsize = 7, useDingbats = FALSE)

  # Linear regression of distance to claudin-low centroid vs. GII
  fit <- lm(patientData$Other_distance ~ patientData$StromalScore)
  palette(subtypeColorsSingle)
  plot(patientData$Other_distance[patientData$ClaudinLow != "ClaudinLow"] ~ patientData$StromalScore[patientData$ClaudinLow != "ClaudinLow"],
       bg = patientData$PAM50[patientData$ClaudinLow != "ClaudinLow"],
       lwd = 0.25,
       pch = 21,
       ylim = c(10, 70),
       xlim = c(min(patientData$StromalScore), max(patientData$StromalScore)),
       font.main = 1,
       cex = 1,
       bty="n",
       xlab = "Stromal score",
       main = "Correlation between stromal score and distance \n to other centroid",
       ylab = "Distance to other centroid")
  
  for(i in (1:nrow(patientData))[patientData$ClaudinLow == "ClaudinLow"]){
    
    points(patientData$Other_distance[i] ~ patientData$StromalScore[i],
           bg = patientData$PAM50[i],
           lwd = 0.25,
           pch = 21,
           cex = 1,
           bty="n")
    
    points(patientData$Other_distance[i] ~ patientData$StromalScore[i],
           col = claudinLowColor,
           pch = 21,
           cex = 1,
           lwd = 1.5)
    
    points(patientData$Other_distance[i] ~ patientData$StromalScore[i],
           col = "black",
           pch = 21,
           cex = 1.3,
           lwd = 0.25)
    
  }
  
  abline(fit, cex = 10)
  text(x = 0, y = max(patientData$Other_distance), labels = paste("R^2 =", round(summary(fit)$r.squared, 2)))
  # The following p value has been checked manually, hard coded for formatting
  text(x = 0, y = 0.95*max(patientData$Other_distance), labels = "P < 0.001")
dev.off()


## CLDist ~ Immune
pdf(file = "./Output/Figures/METABRIC_CLDist_vs_ESTIMATE_Immune_9CLFill.pdf",
    width = 3.5,
    height = 3.5,
    pointsize = 7, useDingbats = FALSE)

  # Linear regression of distance to claudin-low centroid vs. GII
  fit <- lm(patientData$CL_distance ~ patientData$ImmuneScore)
  palette(subtypeColorsSingle)
  plot(patientData$CL_distance[patientData$ClaudinLow != "ClaudinLow"] ~ patientData$ImmuneScore[patientData$ClaudinLow != "ClaudinLow"],
       bg = patientData$PAM50[patientData$ClaudinLow != "ClaudinLow"],
       lwd = 0.25,
       pch = 21,
       ylim = c(10, 70),
       xlim = c(min(patientData$ImmuneScore), max(patientData$ImmuneScore)),
       font.main = 1,
       cex = 1,
       bty="n",
       xlab = "Immune score",
       main = "Correlation between immune score and distance \n to claudin-low centroid",
       ylab = "Distance to claudin-low centroid")
  
  for(i in (1:nrow(patientData))[patientData$ClaudinLow == "ClaudinLow"]){
    
    points(patientData$CL_distance[i] ~ patientData$ImmuneScore[i],
           bg = patientData$PAM50[i],
           lwd = 0.25,
           pch = 21,
           cex = 1,
           bty="n")
    
    points(patientData$CL_distance[i] ~ patientData$ImmuneScore[i],
           col = claudinLowColor,
           pch = 21,
           cex = 1,
           lwd = 1.5)
    
    points(patientData$CL_distance[i] ~ patientData$ImmuneScore[i],
           col = "black",
           pch = 21,
           cex = 1.3,
           lwd = 0.25)
    
  }
  
  abline(fit, cex = 10)
  text(x = 2500, y = max(patientData$CL_distance), labels = paste("R^2 =", round(summary(fit)$r.squared, 2)))
  # The following p value has been checked manually, hard coded for formatting
  text(x = 2500, y = 0.95*max(patientData$CL_distance), labels = "P < 0.001")
dev.off()


## OtherDist ~ Immune
pdf(file = "./Output/Figures/METABRIC_OtherDist_vs_ESTIMATE_Immune_9CLFill.pdf",
    width = 3.5,
    height = 3.5,
    pointsize = 7, useDingbats = FALSE)

  # Linear regression of distance to claudin-low centroid vs. GII
  fit <- lm(patientData$Other_distance ~ patientData$ImmuneScore)
  palette(subtypeColorsSingle)
  plot(patientData$Other_distance[patientData$ClaudinLow != "ClaudinLow"] ~ patientData$ImmuneScore[patientData$ClaudinLow != "ClaudinLow"],
       bg = patientData$PAM50[patientData$ClaudinLow != "ClaudinLow"],
       lwd = 0.25,
       pch = 21,
       ylim = c(10, 70),
       xlim = c(min(patientData$ImmuneScore), max(patientData$ImmuneScore)),
       font.main = 1,
       cex = 1,
       bty="n",
       xlab = "Immune score",
       main = "Correlation between immune score and distance \n to other centroid",
       ylab = "Distance to other centroid")
  
  for(i in (1:nrow(patientData))[patientData$ClaudinLow == "ClaudinLow"]){
    
    points(patientData$Other_distance[i] ~ patientData$ImmuneScore[i],
         bg = patientData$PAM50[i],
         lwd = 0.25,
         pch = 21,
         cex = 1,
         bty="n")
    
    points(patientData$Other_distance[i] ~ patientData$ImmuneScore[i],
           col = claudinLowColor,
           pch = 21,
           cex = 1,
           lwd = 1.5)
    
    points(patientData$Other_distance[i] ~ patientData$ImmuneScore[i],
           col = "black",
           pch = 21,
           cex = 1.3,
           lwd = 0.25)
    
  }
  
  abline(fit, cex = 10)
  text(x = 2500, y = max(patientData$Other_distance), labels = paste("R^2 =", round(summary(fit)$r.squared, 2)))
  # The following p value has been checked manually, hard coded for formatting
  text(x = 2500, y = 0.95*max(patientData$Other_distance), labels = "P < 0.001")
dev.off()

