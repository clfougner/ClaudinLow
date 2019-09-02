# 1) Clone the repository from github.com/clfougner/ClaudinLow
# 2) Download the required data as described in ".README.md" and in "./Data/README.md"
# 3) Install the packages listed in ".README.md"
# 4) Set working directory to "/path/to/ClaudinLow/"
# 5) Run the following scripts:
# Note that the scripts are written so that they can be run independently (given that the required files are available), so there is some redundancy in the scripts.

# METABRIC
source("./Code/METABRIC_patientData.r") # Prepare required data
source("./Code/METABRIC_ClaudinLowPAM50_Figures.r") # Make figures for claudin-low tumors stratified by intrinsic subtype
source("./Code/METABRIC_CoreClaudinLow.r") # Core claudin-low heatmap; add CoreCL to patientData
source("./Code/METABRIC_CoreClaudinLowFigures.r") # Make figures related to core claudin-low tumors
source("./Code/METABRIC_CoreClaudinLowPAM50_Figures.r") # Make figures analogous to METABRIC_ClaudinLowPAM50_Figures.r for CoreCL
source("./Code/METABRIC_ClaudinLow_Mutations_CNAs.r") # Make files with overview of all mutations and CNAs, with significance testing

# OSLO2
source("./Code/Oslo2_CoreClaudinLow.r") # Prepared required data, identify claudin-low and CoreCL tumors
source("./Code/Oslo2_CoreClaudinLowFigures.r") # Make figures of CoreCL tumors

# TCGA
source("./Code/TCGA.r") # Prepared required data, identify claudin-low and CoreCL tumors

# Other
source("./Code/CellularityVsPrevalence.r") # Plot the prevalence of claudin-low tumors versus their prevalence in a cohort 
