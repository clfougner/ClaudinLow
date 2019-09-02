##################################################
## Check mutations and CNAs in claudin-low and non-claudin-low tumors
##
## Requires the output from ./Code/METABRIC_CoreClaudinLow.r (./Output/METABRIC_patientData_CoreClaudinLow.txt. Make sure that path is set to /path/to/ClaudinLow/
##################################################

print("METABRIC: Checking for significantly over-represented mutations and CNAs in claudin-low tumors stratified by intrinsic subtype")

########################################
## Read input data

## PatientData
# Read patient data output from ./Code/METABRIC_OutputAllDataPerSample.r
patientData <- read.table(file = "./Output/METABRIC_patientData_CoreClaudinLow.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

## Mutations
# Read mutations
mutations <- read.table(file = "./Data/brca_metabric/data_mutations_extended.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE, quote = "")

# Remove mutations found in patients without all required data
mutations <- mutations[mutations$Tumor_Sample_Barcode %in% patientData$PATIENT_ID, ]


## CNAs
# Read CNAs
cna <- read.table(file = "./Data/brca_metabric/data_CNA.txt", sep = "\t", header = TRUE)

# Set gene symbols as row names
cna <- data.frame(cna[, 3:ncol(cna)], row.names = cna$Hugo_Symbol)

# Rename samples to match with other data types
colnames(cna) <- gsub("\\.", "-", colnames(cna))

# Remove samples without all other required data
cna <- cna[, colnames(cna) %in% patientData$PATIENT_ID]

# Re-order samples to match other data
cna <- cna[, order(colnames(cna))]

# Read genes associated with CNA in cancer (https://cancer.sanger.ac.uk/census)
cnaGenes <- read.table(file = "./ReferenceFiles/CNAGenes.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Separate into gains and losses
cnaGain <- cna[cnaGenes$GeneSymbol[cnaGenes$Class == "Amplification"], ]
cnaLosses <- cna[cnaGenes$GeneSymbol[cnaGenes$Class == "Deletion"], ]

# Set all gains = 1, and remove losses in genes associated with gains in cancer
cnaGain[cnaGain == 2] <- 1
cnaGain[cnaGain == -2] <- 0
cnaGain[cnaGain == -1] <- 0

# Set all losses = 1, and remove gains in genes associated with losses in cancer
cnaLosses[cnaLosses == 1] <- 0
cnaLosses[cnaLosses == 2] <- 0
cnaLosses[cnaLosses == -2] <- 1
cnaLosses[cnaLosses == -1] <- 1

# Gather gains and losses in a data frame
allCNA <- rbind(cnaGain, cnaLosses)


########################################
## Mutations: Basal-like

# Create an empty data frame for comparisons
basalMutations <- data.frame(Gene = unique(mutations$Hugo_Symbol), BasalPos = NA, BasalNeg = NA, BasalProportion = NA, BasalClaudinLowPos = NA, BasalClaudinLowNeg = NA, BasalClaudinLowProportion = NA, CL_to_NonCL_Ratio = NA, p = NA, p.adj = NA, row.names = unique(mutations$Hugo_Symbol))

# Find number of basal-like and basal-like claudin-low tumors
basalTumorNumber <- table(patientData$PAM50ClaudinLow)[["Basal"]]
basalClaudinLowTumorNumber <- table(patientData$PAM50ClaudinLow)[["BasalClaudinLow"]]

# Find sample IDs for basal-like and basal-like claudin-low tumors
basalSamples <- as.character(patientData$PATIENT_ID[patientData$PAM50ClaudinLow == "Basal"])
basalClaudinLowSamples <- as.character(patientData$PATIENT_ID[patientData$PAM50ClaudinLow == "BasalClaudinLow"])

for(i in row.names(basalMutations)){

  # Subset mutations in the given gene for basal-like tumors
  currentBasalMutations <- mutations[mutations$Hugo_Symbol == i & mutations$Tumor_Sample_Barcode %in% basalSamples, ]

  # Find number of basal-like tumors with mutations in the given gene
  basalMutations[i, "BasalPos"] <- nrow(currentBasalMutations)

  # Find number of basal-like tumors without mutations in the given gene
  basalMutations[i, "BasalNeg"] <- basalTumorNumber - nrow(currentBasalMutations)

  # Find the proportion of basal-like tumors with mutations in the given gene
  basalMutations[i, "BasalProportion"] <- basalMutations[i, "BasalPos"]/basalTumorNumber


  # Subset mutations in the given gene for basal-like claudin-low tumors
  currentBasalClaudinLowMutations <- mutations[mutations$Hugo_Symbol == i & mutations$Tumor_Sample_Barcode %in% basalClaudinLowSamples, ]

  # Find number of basal-like claudin-low tumors with mutations in the given gene
  basalMutations[i, "BasalClaudinLowPos"] <- nrow(currentBasalClaudinLowMutations)

  # Find number of basal-like claudin-low tumors without mutations in the given gene
  basalMutations[i, "BasalClaudinLowNeg"] <- basalClaudinLowTumorNumber - nrow(currentBasalClaudinLowMutations)

  # Find the proportion of basal-like claudin-low tumors with mutations in the given gene
  basalMutations[i, "BasalClaudinLowProportion"] <- basalMutations[i, "BasalClaudinLowPos"]/basalClaudinLowTumorNumber


  # Find the ratio of mutations in basal-like claudin-low tumors to mutations in basal-like non-claudin-low tumors
  basalMutations[i, "CL_to_NonCL_Ratio"] <- basalMutations[i, "BasalClaudinLowProportion"]/basalMutations[i, "BasalProportion"]

  # Test significance of the difference between the two
  test <- fisher.test(x = matrix(data = c(basalMutations[i, "BasalPos"], basalMutations[i, "BasalNeg"],
                                          basalMutations[i, "BasalClaudinLowPos"], basalMutations[i, "BasalClaudinLowNeg"]),
                                 nrow = 2, byrow = FALSE))

  # Insert the p value into the data frame
  basalMutations[i, "p"] <- test$p.value
}

# Adjust for multiple hypothesis testing
basalMutations[, "p.adj"] <- p.adjust(p = basalMutations[, "p"], method = "bonferroni")

# Write to file
write.table(x = basalMutations, file = "./Output/METABRIC_ClaudinLowBasal_vs_Basal_mutations.txt", sep = "\t", quote = FALSE)

########################################
## Mutations: LumA

# Create an empty data frame for comparisons
lumAMutations <- data.frame(Gene = unique(mutations$Hugo_Symbol), LumAPos = NA, LumANeg = NA, LumAProportion = NA, LumAClaudinLowPos = NA, LumAClaudinLowNeg = NA, LumAClaudinLowProportion = NA, CL_to_NonCL_Ratio = NA, p = NA, p.adj = NA, row.names = unique(mutations$Hugo_Symbol))

# Find number of LumA and LumA claudin-low tumors
lumATumorNumber <- table(patientData$PAM50ClaudinLow)[["LumA"]]
lumAClaudinLowTumorNumber <- table(patientData$PAM50ClaudinLow)[["LumAClaudinLow"]]

# Find sample IDs for LumA and LumA claudin-low tumors
lumASamples <- as.character(patientData$PATIENT_ID[patientData$PAM50ClaudinLow == "LumA"])
lumAClaudinLowSamples <- as.character(patientData$PATIENT_ID[patientData$PAM50ClaudinLow == "LumAClaudinLow"])

for(i in row.names(lumAMutations)){

  # Subset mutations in the given gene for LumA tumors
  currentLumAMutations <- mutations[mutations$Hugo_Symbol == i & mutations$Tumor_Sample_Barcode %in% lumASamples, ]

  # Find number of LumA tumors with mutations in the given gene
  lumAMutations[i, "LumAPos"] <- nrow(currentLumAMutations)

  # Find number of LumA tumors without mutations in the given gene
  lumAMutations[i, "LumANeg"] <- lumATumorNumber - nrow(currentLumAMutations)

  # Find the proportion of LumA tumors with mutations in the given gene
  lumAMutations[i, "LumAProportion"] <- lumAMutations[i, "LumAPos"]/lumATumorNumber


  # Subset mutations in the given gene for LumA claudin-low tumors
  currentLumAClaudinLowMutations <- mutations[mutations$Hugo_Symbol == i & mutations$Tumor_Sample_Barcode %in% lumAClaudinLowSamples, ]

  # Find number of LumA claudin-low tumors with mutations in the given gene
  lumAMutations[i, "LumAClaudinLowPos"] <- nrow(currentLumAClaudinLowMutations)

  # Find number of LumA claudin-low tumors without mutations in the given gene
  lumAMutations[i, "LumAClaudinLowNeg"] <- lumAClaudinLowTumorNumber - nrow(currentLumAClaudinLowMutations)

  # Find the proportion of LumA claudin-low tumors with mutations in the given gene
  lumAMutations[i, "LumAClaudinLowProportion"] <- lumAMutations[i, "LumAClaudinLowPos"]/lumAClaudinLowTumorNumber

  # Find the ratio of mutations in LumA claudin-low tumors to mutations in LumA non-claudin-low tumors
  lumAMutations[i, "CL_to_NonCL_Ratio"] <- lumAMutations[i, "LumAClaudinLowProportion"]/lumAMutations[i, "LumAProportion"]

  # Test significance of the difference between the two
  test <- fisher.test(x = matrix(data = c(lumAMutations[i, "LumAPos"], lumAMutations[i, "LumANeg"],
                                          lumAMutations[i, "LumAClaudinLowPos"], lumAMutations[i, "LumAClaudinLowNeg"]),
                                 nrow = 2, byrow = FALSE))

  # Insert the p value into the data frame
  lumAMutations[i, "p"] <- test$p.value
}

# Adjust for multiple hypothesis testing
lumAMutations[, "p.adj"] <- p.adjust(p = lumAMutations[, "p"], method = "bonferroni")

# Write to file
write.table(x = lumAMutations, file = "./Output/METABRIC_ClaudinLowLumA_vs_LumA_mutations.txt", sep = "\t", quote = FALSE)


########################################
## Mutations: Normal-like

# Create an empty data frame for comparisons
normalMutations <- data.frame(Gene = unique(mutations$Hugo_Symbol), NormalPos = NA, NormalNeg = NA, NormalProportion = NA, NormalClaudinLowPos = NA, NormalClaudinLowNeg = NA, NormalClaudinLowProportion = NA, CL_to_NonCL_Ratio = NA, p = NA, p.adj = NA, row.names = unique(mutations$Hugo_Symbol))

# Find number of normal-like and normal-like claudin-low tumors
normalTumorNumber <- table(patientData$PAM50ClaudinLow)[["Normal"]]
normalClaudinLowTumorNumber <- table(patientData$PAM50ClaudinLow)[["NormalClaudinLow"]]

# Find sample IDs for normal-like and normal-like claudin-low tumors
normalSamples <- as.character(patientData$PATIENT_ID[patientData$PAM50ClaudinLow == "Normal"])
normalClaudinLowSamples <- as.character(patientData$PATIENT_ID[patientData$PAM50ClaudinLow == "NormalClaudinLow"])

for(i in row.names(normalMutations)){

  # Subset mutations in the given gene for normal-like tumors
  currentNormalMutations <- mutations[mutations$Hugo_Symbol == i & mutations$Tumor_Sample_Barcode %in% normalSamples, ]

  # Find number of normal-like tumors with mutations in the given gene
  normalMutations[i, "NormalPos"] <- nrow(currentNormalMutations)

  # Find number of normal-like tumors without mutations in the given gene
  normalMutations[i, "NormalNeg"] <- normalTumorNumber - nrow(currentNormalMutations)

  # Find the proportion of normal-like tumors with mutations in the given gene
  normalMutations[i, "NormalProportion"] <- normalMutations[i, "NormalPos"]/normalTumorNumber

  # Subset mutations in the given gene for normal-like claudin-low tumors
  currentNormalClaudinLowMutations <- mutations[mutations$Hugo_Symbol == i & mutations$Tumor_Sample_Barcode %in% normalClaudinLowSamples, ]

  # Find number of normal-like claudin-low tumors with mutations in the given gene
  normalMutations[i, "NormalClaudinLowPos"] <- nrow(currentNormalClaudinLowMutations)

  # Find number of normal-like claudin-low tumors without mutations in the given gene
  normalMutations[i, "NormalClaudinLowNeg"] <- normalClaudinLowTumorNumber - nrow(currentNormalClaudinLowMutations)

  # Find the proportion of normal-like claudin-low tumors with mutations in the given gene
  normalMutations[i, "NormalClaudinLowProportion"] <- normalMutations[i, "NormalClaudinLowPos"]/normalClaudinLowTumorNumber


  # Find the ratio of mutations in normal-like claudin-low tumors to mutations in normal-like non-claudin-low tumors
  normalMutations[i, "CL_to_NonCL_Ratio"] <- normalMutations[i, "NormalClaudinLowProportion"]/normalMutations[i, "NormalProportion"]

  # Test significance of the difference between the two
  test <- fisher.test(x = matrix(data = c(normalMutations[i, "NormalPos"], normalMutations[i, "NormalNeg"],
                                          normalMutations[i, "NormalClaudinLowPos"], normalMutations[i, "NormalClaudinLowNeg"]),
                                 nrow = 2, byrow = FALSE))

  # Insert the p value into the data frame
  normalMutations[i, "p"] <- test$p.value
}

# Adjust for multiple hypothesis testing
normalMutations[, "p.adj"] <- p.adjust(p = normalMutations[, "p"], method = "bonferroni")

# Write to file
write.table(x = normalMutations, file = "./Output/METABRIC_ClaudinLowNormal_vs_Normal_mutations.txt", sep = "\t", quote = FALSE)


########################################
## Mutations: Basal CoreCL

# Create an empty data frame for comparisons
basalMutations <- data.frame(Gene = unique(mutations$Hugo_Symbol), BasalPos = NA, BasalNeg = NA, BasalProportion = NA, BasalClaudinLowPos = NA, BasalClaudinLowNeg = NA, BasalClaudinLowProportion = NA, CL_to_NonCL_Ratio = NA, p = NA, p.adj = NA, row.names = unique(mutations$Hugo_Symbol))

# Find number of basal-like+basalOtherCL and basal-like CoreCL tumors
basalTumorNumber <- table(patientData$PAM50toCoreCL)[["Basal"]] + table(patientData$PAM50toCoreCL)[["BasalOtherClaudinLow"]]
basalClaudinLowTumorNumber <- table(patientData$PAM50toCoreCL)[["BasalCoreClaudinLow"]]

# Find sample IDs for basal-like+basalOtherCL and basal-like CoreCL claudin-low tumors
basalSamples <- as.character(patientData$PATIENT_ID[patientData$PAM50toCoreCL == "Basal" | patientData$PAM50toCoreCL == "BasalOtherClaudinLow"])
basalClaudinLowSamples <- as.character(patientData$PATIENT_ID[patientData$PAM50toCoreCL == "BasalCoreClaudinLow"])

for(i in row.names(basalMutations)){

  # Subset mutations in the given gene for basal-like+basalOtherCL tumors
  currentBasalMutations <- mutations[mutations$Hugo_Symbol == i & mutations$Tumor_Sample_Barcode %in% basalSamples, ]

  # Find number of basal-like+basalOtherCL tumors with mutations in the given gene
  basalMutations[i, "BasalPos"] <- nrow(currentBasalMutations)

  # Find number of basal-like+basalOtherCL tumors without mutations in the given gene
  basalMutations[i, "BasalNeg"] <- basalTumorNumber - nrow(currentBasalMutations)

  # Find the proportion of basal-like+basalOtherCL tumors with mutations in the given gene
  basalMutations[i, "BasalProportion"] <- basalMutations[i, "BasalPos"]/basalTumorNumber

  # Subset mutations in the given gene for basal-like CoreCL tumors
  currentBasalClaudinLowMutations <- mutations[mutations$Hugo_Symbol == i & mutations$Tumor_Sample_Barcode %in% basalClaudinLowSamples, ]

  # Find number of basal-like CoreCL tumors with mutations in the given gene
  basalMutations[i, "BasalClaudinLowPos"] <- nrow(currentBasalClaudinLowMutations)

  # Find number of basal-like CoreCL tumors without mutations in the given gene
  basalMutations[i, "BasalClaudinLowNeg"] <- basalClaudinLowTumorNumber - nrow(currentBasalClaudinLowMutations)

  # Find the proportion of basal-like CoreCL tumors with mutations in the given gene
  basalMutations[i, "BasalClaudinLowProportion"] <- basalMutations[i, "BasalClaudinLowPos"]/basalClaudinLowTumorNumber


  # Find the ratio of mutations in basal-like CoreCL tumors to mutations in basal-like+basalOtherCL tumors
  basalMutations[i, "CL_to_NonCL_Ratio"] <- basalMutations[i, "BasalClaudinLowProportion"]/basalMutations[i, "BasalProportion"]

  # Test significance of the difference between the two
  test <- fisher.test(x = matrix(data = c(basalMutations[i, "BasalPos"], basalMutations[i, "BasalNeg"],
                                          basalMutations[i, "BasalClaudinLowPos"], basalMutations[i, "BasalClaudinLowNeg"]),
                                 nrow = 2, byrow = FALSE))

  # Insert the p value into the data frame
  basalMutations[i, "p"] <- test$p.value
}

# Adjust for multiple hypothesis testing
basalMutations[, "p.adj"] <- p.adjust(p = basalMutations[, "p"], method = "bonferroni")

# Write to file
write.table(x = basalMutations, file = "./Output/METABRIC_CoreClaudinLowBasal_vs_BasalNonCLOtherCL_mutations.txt", sep = "\t", quote = FALSE)


########################################
## CNA: Basal-like

# Create an empty data frame for comparisons
basalCNA <- data.frame(Gene = row.names(allCNA), Class = c(rep("Gain", times = nrow(cnaGain)), rep("Loss", times = nrow(cnaLosses))), BasalPos = NA, BasalNeg = NA, BasalProportion = NA, BasalClaudinLowPos = NA, BasalClaudinLowNeg = NA, BasalClaudinLowProportion = NA, CL_to_NonCL_Ratio = NA, p = NA, p.adj = NA, row.names = row.names(allCNA))

for(i in row.names(basalCNA)){

  # Subset CNAs in the given gene for basal-like tumors
  currentBasalCNA <- allCNA[i, patientData$PATIENT_ID[patientData$PAM50ClaudinLow == "Basal"]]

  # Find number of basal-like tumors with CNA in the given gene
  basalCNA[i, "BasalPos"] <- sum(currentBasalCNA)

  # Find number of basal-like tumors without CNA in the given gene
  basalCNA[i, "BasalNeg"] <- length(currentBasalCNA) - sum(currentBasalCNA)

  # Find the proportion of basal-like tumors with CNA in the given gene
  basalCNA[i, "BasalProportion"] <- basalCNA[i, "BasalPos"]/length(currentBasalCNA)


  # Subset CNAs in the given gene for basal-like claudin-low tumors
  currentBasalClaudinLowCNA <- allCNA[i, patientData$PATIENT_ID[patientData$PAM50ClaudinLow == "BasalClaudinLow"]]

  # Find number of basal-like claudin-low tumors with CNA in the given gene
  basalCNA[i, "BasalClaudinLowPos"] <- sum(currentBasalClaudinLowCNA)

  # Find number of basal-like claudin-low tumors without CNA in the given gene
  basalCNA[i, "BasalClaudinLowNeg"] <- length(currentBasalClaudinLowCNA) - sum(currentBasalClaudinLowCNA)

  # Find the proportion of basal-like claudin-low tumors with CNA in the given gene
  basalCNA[i, "BasalClaudinLowProportion"] <- basalCNA[i, "BasalClaudinLowPos"]/length(currentBasalClaudinLowCNA)


  # Find the ratio of CNA in basal-like claudin-low tumors to CNA in basal-like non-claudin-low tumors
  basalCNA[i, "CL_to_NonCL_Ratio"] <- basalCNA[i, "BasalClaudinLowProportion"]/basalCNA[i, "BasalProportion"]

  # Test significance of the difference between the two
  test <- fisher.test(x = matrix(data = c(basalCNA[i, "BasalPos"], basalCNA[i, "BasalNeg"],
                                          basalCNA[i, "BasalClaudinLowPos"], basalCNA[i, "BasalClaudinLowNeg"]),
                                 nrow = 2, byrow = FALSE))

  # Insert the p value into the data frame
  basalCNA[i, "p"] <- test$p.value
}

# Adjust for multiple hypothesis testing
basalCNA[, "p.adj"] <- p.adjust(p = basalCNA[, "p"], method = "bonferroni")

# Write to file
write.table(x = basalCNA, file = "./Output/METABRIC_ClaudinLowBasal_vs_Basal_CNA.txt", sep = "\t", quote = FALSE)


########################################
## CNA: LumA

# Create an empty data frame for comparisons
lumACNA <- data.frame(Gene = row.names(allCNA), Class = c(rep("Gain", times = nrow(cnaGain)), rep("Loss", times = nrow(cnaLosses))), LumAPos = NA, LumANeg = NA, LumAProportion = NA, LumAClaudinLowPos = NA, LumAClaudinLowNeg = NA, LumAClaudinLowProportion = NA, CL_to_NonCL_Ratio = NA, p = NA, p.adj = NA, row.names = row.names(allCNA))

for(i in row.names(lumACNA)){

  # Subset CNAs in the given gene for LumA tumors
  currentLumACNA <- allCNA[i, patientData$PATIENT_ID[patientData$PAM50ClaudinLow == "LumA"]]

  # Find number of LumA tumors with CNA in the given gene
  lumACNA[i, "LumAPos"] <- sum(currentLumACNA)

  # Find number of LumA tumors without CNA in the given gene
  lumACNA[i, "LumANeg"] <- length(currentLumACNA) - sum(currentLumACNA)

  # Find the proportion of LumA tumors with CNA in the given gene
  lumACNA[i, "LumAProportion"] <- lumACNA[i, "LumAPos"]/length(currentLumACNA)

  # Subset CNAs in the given gene for LumA claudin-low tumors
  currentLumAClaudinLowCNA <- allCNA[i, patientData$PATIENT_ID[patientData$PAM50ClaudinLow == "LumAClaudinLow"]]

  # Find number of LumA claudin-low tumors with CNA in the given gene
  lumACNA[i, "LumAClaudinLowPos"] <- sum(currentLumAClaudinLowCNA)

  # Find number of LumA claudin-low tumors without CNA in the given
  lumACNA[i, "LumAClaudinLowNeg"] <- length(currentLumAClaudinLowCNA) - sum(currentLumAClaudinLowCNA)

  # Find the proportion of LumA claudin-low tumors with CNA in the given gene
  lumACNA[i, "LumAClaudinLowProportion"] <- lumACNA[i, "LumAClaudinLowPos"]/length(currentLumAClaudinLowCNA)

  # Find the ratio of CNA in LumA claudin-low tumors to CNA in LumA non-claudin-low tumors
  lumACNA[i, "CL_to_NonCL_Ratio"] <- lumACNA[i, "LumAClaudinLowProportion"]/lumACNA[i, "LumAProportion"]

  # Test significance of the difference between the two
  test <- fisher.test(x = matrix(data = c(lumACNA[i, "LumAPos"], lumACNA[i, "LumANeg"],
                                          lumACNA[i, "LumAClaudinLowPos"], lumACNA[i, "LumAClaudinLowNeg"]),
                                 nrow = 2, byrow = FALSE))

  # Insert the p value into the data frame
  lumACNA[i, "p"] <- test$p.value
}

# Adjust for multiple hypothesis testing
lumACNA[, "p.adj"] <- p.adjust(p = lumACNA[, "p"], method = "bonferroni")

# Write to file
write.table(x = lumACNA, file = "./Output/METABRIC_ClaudinLowLumA_vs_LumA_CNA.txt", sep = "\t", quote = FALSE)


########################################
## CNA: Normal-like

# Create an empty data frame for comparisons
normalCNA <- data.frame(Gene = row.names(allCNA), Class = c(rep("Gain", times = nrow(cnaGain)), rep("Loss", times = nrow(cnaLosses))), NormalPos = NA, NormalNeg = NA, NormalProportion = NA, NormalClaudinLowPos = NA, NormalClaudinLowNeg = NA, NormalClaudinLowProportion = NA, CL_to_NonCL_Ratio = NA, p = NA, p.adj = NA, row.names = row.names(allCNA))

for(i in row.names(normalCNA)){

  # Subset CNAs in the given gene for normal-like tumors
  currentNormalCNA <- allCNA[i, patientData$PATIENT_ID[patientData$PAM50ClaudinLow == "Normal"]]

  # Find number of normal-like tumors with CNA in the given gene
  normalCNA[i, "NormalPos"] <- sum(currentNormalCNA)

  # Find number of normal-like tumors without CNA in the given gene
  normalCNA[i, "NormalNeg"] <- length(currentNormalCNA) - sum(currentNormalCNA)

  # Find the proportion of normal-like tumors with CNA in the given gene
  normalCNA[i, "NormalProportion"] <- normalCNA[i, "NormalPos"]/length(currentNormalCNA)

  # Subset CNAs in the given gene for normal-like claudin-low tumors
  currentNormalClaudinLowCNA <- allCNA[i, patientData$PATIENT_ID[patientData$PAM50ClaudinLow == "NormalClaudinLow"]]

  # Find number of normal-like claudin-low tumors with CNA in the given gene
  normalCNA[i, "NormalClaudinLowPos"] <- sum(currentNormalClaudinLowCNA)

  # Find number of normal-like claudin-low tumors without CNA in the given gene
  normalCNA[i, "NormalClaudinLowNeg"] <- length(currentNormalClaudinLowCNA) - sum(currentNormalClaudinLowCNA)

  # Find the proportion of normal-like claudin-low tumors with CNA in the given gene
  normalCNA[i, "NormalClaudinLowProportion"] <- normalCNA[i, "NormalClaudinLowPos"]/length(currentNormalClaudinLowCNA)

  # Find the ratio of CNA in normal-like claudin-low tumors to CNA in normal-like non-claudin-low tumors
  normalCNA[i, "CL_to_NonCL_Ratio"] <- normalCNA[i, "NormalClaudinLowProportion"]/normalCNA[i, "NormalProportion"]

  # Test significance of the difference between the two
  test <- fisher.test(x = matrix(data = c(normalCNA[i, "NormalPos"], normalCNA[i, "NormalNeg"],
                                          normalCNA[i, "NormalClaudinLowPos"], normalCNA[i, "NormalClaudinLowNeg"]),
                                 nrow = 2, byrow = FALSE))

  # Insert the p value into the data frame
  normalCNA[i, "p"] <- test$p.value
}

# Adjust for multiple hypothesis testing
normalCNA[, "p.adj"] <- p.adjust(p = normalCNA[, "p"], method = "bonferroni")

# Write to file
write.table(x = normalCNA, file = "./Output/METABRIC_ClaudinLowNormal_vs_Normal_CNA.txt", sep = "\t", quote = FALSE)


########################################
## CNA: Basal CoreCL

# Create an empty data frame for comparisons
basalCNA <- data.frame(Gene = row.names(allCNA), Class = c(rep("Gain", times = nrow(cnaGain)), rep("Loss", times = nrow(cnaLosses))), BasalPos = NA, BasalNeg = NA, BasalProportion = NA, BasalClaudinLowPos = NA, BasalClaudinLowNeg = NA, BasalClaudinLowProportion = NA, CL_to_NonCL_Ratio = NA, p = NA, p.adj = NA, row.names = row.names(allCNA))

for(i in row.names(basalCNA)){

  # Subset CNAs in the given gene for basal-like+basalOtherCL tumors
  currentBasalCNA <- allCNA[i, patientData$PATIENT_ID[patientData$PAM50toCoreCL == "Basal" | patientData$PAM50toCoreCL == "BasalOtherClaudinLow"]]

  # Find number of basal-like+basalOtherCL tumors with CNA in the given gene
  basalCNA[i, "BasalPos"] <- sum(currentBasalCNA)

  # Find number of basal-like+basalOtherCL tumors without CNA in the given gene
  basalCNA[i, "BasalNeg"] <- length(currentBasalCNA) - sum(currentBasalCNA)

  # Find the proportion of basal-like+basalOtherCL tumors with CNA in the given gene
  basalCNA[i, "BasalProportion"] <- basalCNA[i, "BasalPos"]/length(currentBasalCNA)


  # Subset CNAs in the given gene for basal-like CoreCL tumors
  currentBasalClaudinLowCNA <- allCNA[i, patientData$PATIENT_ID[patientData$PAM50toCoreCL == "BasalCoreClaudinLow"]]

  # Find number of basal-like CoreCL tumors with CNA in the given gene
  basalCNA[i, "BasalClaudinLowPos"] <- sum(currentBasalClaudinLowCNA)

  # Find number of basal-like CoreCL tumors without CNA in the given gene
  basalCNA[i, "BasalClaudinLowNeg"] <- length(currentBasalClaudinLowCNA) - sum(currentBasalClaudinLowCNA)

  # Find the proportion of basal-like CoreCL tumors with CNA in the given gene
  basalCNA[i, "BasalClaudinLowProportion"] <- basalCNA[i, "BasalClaudinLowPos"]/length(currentBasalClaudinLowCNA)

  # Find the ratio of CNA in basal-like CoreCL tumors to CNA in basal-like non-claudin-low tumors
  basalCNA[i, "CL_to_NonCL_Ratio"] <- basalCNA[i, "BasalClaudinLowProportion"]/basalCNA[i, "BasalProportion"]

  # Test significance of the difference between the two
  test <- fisher.test(x = matrix(data = c(basalCNA[i, "BasalPos"], basalCNA[i, "BasalNeg"],
                                          basalCNA[i, "BasalClaudinLowPos"], basalCNA[i, "BasalClaudinLowNeg"]),
                                 nrow = 2, byrow = FALSE))

  # Insert the p value into the data frame
  basalCNA[i, "p"] <- test$p.value
}

# Adjust for multiple hypothesis testing
basalCNA[, "p.adj"] <- p.adjust(p = basalCNA[, "p"], method = "bonferroni")

# Write to file
write.table(x = basalCNA, file = "./Output/METABRIC_CoreClaudinLowBasal_vs_BasalNonCLOtherCL_CNA.txt", sep = "\t", quote = FALSE)

print("METABRIC: Finished checking for significantly over-represented mutations and CNAs in claudin-low tumors stratified by intrinsic subtype")
