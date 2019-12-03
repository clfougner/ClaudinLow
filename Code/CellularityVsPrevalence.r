##################################################
## Claudin-low prevalence vs cellularity cut-off in all cohorts
##
## Requires the output from ./Code/METABRIC_patientData.r, ./Code/Oslo2_CoreClaudinLow.r and ./Code/TCGA.r
## Make sure that working directory is set to /path/to/ClaudinLow/
##################################################

print("Plotting claudin-low prevalence in a cohort against tumor cellularity cutoff")

##################################################
## TCGA

# Read patientData from TCGA-BRCA
tcga <- read.table(file = "./Output/TCGA-BRCA_patientData.txt", header = TRUE, sep = "\t")

# Find number of claudin-low tumors
tcgaClaudinLow <- table(tcga$ClaudinLow)[["ClaudinLow"]]

# Find proportion of all tumors which are claudin-low
tcgaProp <- tcgaClaudinLow/nrow(tcga)


##################################################
## METABRIC
# Read patientData from METABRIC
mb <- read.table(file = "./Output/METABRIC_patientData_CoreClaudinLow.txt", header = TRUE, sep = "\t")

# Read original annotation files
mb_disc <- read.table(file = "./Data/table_S2_revised.txt", sep = "\t", header = TRUE)
mb_val <- read.table(file = "./Data/table_S3_revised.txt", sep = "\t", header = TRUE)

# Find number of claudin-low tumors - discovery
mbClaudinLow_disc <- table(mb$ClaudinLow[mb$PATIENT_ID %in% mb_disc$METABRIC_ID])[["ClaudinLow"]]

# Find proportion of all tumors which are claudin-low - discovery
mbProp_disc <- mbClaudinLow_disc/table(mb$PATIENT_ID %in% mb_disc$METABRIC_ID)[["TRUE"]]

# Find number of claudin-low tumors - validation
mbClaudinLow_val <- table(mb$ClaudinLow[mb$PATIENT_ID %in% mb_val$METABRIC_ID])[["ClaudinLow"]]
# Find proportion of all tumors which are claudin-low - validation
mbProp_val <- mbClaudinLow_val/table(mb$PATIENT_ID %in% mb_val$METABRIC_ID)[["TRUE"]]


## Update patientData with cohort
mb$Cohort[mb$PATIENT_ID %in% mb_disc$METABRIC_ID] <- "Discovery"
mb$Cohort[mb$PATIENT_ID %in% mb_val$METABRIC_ID] <- "Validation"

write.table(mb, file = "./Output/METABRIC_patientData_CoreClaudinLow.txt", sep = "\t", row.names = FALSE, quote = FALSE)

##################################################
## Oslo2
# Read patientData from Oslo2
oslo2 <- read.table(file = "./Output/Oslo2_patientData.txt", sep = "\t", header = TRUE)

# Find number of claudin-low tumors
osloClaudinLow <- table(oslo2$ClaudinLow)[["ClaudinLow"]]

# Find proportion of all tumors which are claudin-low
oslo2prop <- osloClaudinLow/nrow(oslo2)


##################################################
## Gather data and plot

# Cellularity cut-offs in the cohorts are 0.6 (TCGA), 0.4 (METABRIC) and 0 (Oslo2). Data gathered from the respective cohorts' associated publications
prevalences <- matrix(data = c(0.6, 0.4, 0, 0,
                               tcgaProp, mbProp_disc, mbProp_val, oslo2prop),
                      nrow = 2, byrow = TRUE)

prevalences <- t(prevalences)
row.names(prevalences) <- c("TCGA-BRCA", "MB-Disc.", "MB-Val.", "Oslo2")

pdf("./Output/Figures/ClaudinLowPrevalences.pdf", height = 4, width = 4)
  plot(prevalences,
       col="lightblue",
       pch = 3,
       cex = 2,
       bty="n",
       ylim = c(0, 0.1),
       xlab = "Cellularity cut-off",
       ylab = "Claudin-low proportion")
  
  text(x = c(prevalences[1:2, 1] * c(0.85, 1), 0.06, 0.03),
       y = prevalences[, 2] * c(0.7, 0.8, 0.85, 0.85),
       labels=rownames(prevalences),
       cex = 0.9, font = 2)
dev.off()

print("Finished all analyses.")