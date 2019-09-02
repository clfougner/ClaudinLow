# Re-definition of claudin-low as a breast cancer phenotype
Repository containing all analyses in Fougner et al. 2019.

## Data
Here, we use three datasets, METABRIC, OSLO2 and TCGA-BRCA. See `./Data/README.md` for details.

### METABRIC
The METABRIC dataset was first described by [Curtis et al. 2012](https://www.nature.com/articles/nature10983) and later by [Pereira et al. 2016](https://www.nature.com/articles/ncomms11479). The following files are used:

* Processed data was downloaded from [cBioportal](http://www.cbioportal.org/study?id=brca_metabric) and saved to `./Data/brca_metabric/`.

* Supplementary files 2 and 3 from [Curtis et al. 2012](https://www.nature.com/articles/nature10983) are saved to `./Data/table_S2_revised.txt` and `./Data/table_S3_revised.txt`.

* Segmented copy number data from [Pereira et al. 2016](https://www.nature.com/articles/ncomms11479) [(available here)](https://raw.githubusercontent.com/cclab-brca/mutationalProfiles/master/Data/ascatSegments_withoutCNVs.txt) is saved to `./Data/ascatSegments_withoutCNVs.txt`; IDs are matched with [this file](https://raw.githubusercontent.com/cclab-brca/mutationalProfiles/master/Data/tumorIdMap.txt), saved to `./Data/tumorIdMap.txt`. These files are re-hosted in this repository.

### Oslo2
Gene expression data, hormone receptor status, and intrinsic (PAM50) subtype from the [Oslo2 cohort](https://breast-cancer-research.biomedcentral.com/articles/10.1186/s13058-017-0812-y) is downloaded from [GEO accession GSE80999](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80999) using the getGEO command in R in the `./Code/Oslo2_CoreClaudinLow.r` script; this is done automatically and the data does not need to be downloaded from the linked repository.

### TCGA-BRCA
Gene expression data from the [TCGA cohort](https://www.nature.com/articles/nature11412) is downloaded from [cBioportal](http://www.cbioportal.org/study?id=brca_tcga_pan_can_atlas_2018/) and saved to `./Data/brca_tcga_pan_can_atlas_2018`.

## Reference Files
All reference files are included in this repository, and will be downloaded by cloning it.

## Analyses
Scripts for all analyses are saved to `./Code/`. All scripts should be run with this directory set as the working directory (`./`).

All analyses were performed in [R](https://www.r-project.org) version 3.6.0 on macOS version 10.14.3, using the following packages:

* [gtools](https://cran.r-project.org/web/packages/gtools/index.html)
* [survival](https://github.com/therneau/survival)
* [survminer](https://rpkgs.datanovia.com/survminer/index.html)
* [circlize](https://cran.r-project.org/package=circlize)[(publication)](https://academic.oup.com/bioinformatics/article/30/19/2811/2422259)
* [ggplot2](https://ggplot2.tidyverse.org)
* [ggsignif](https://github.com/const-ae/ggsignif)
* [ggsci](https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html#discussion)
* [SigClust](https://cran.r-project.org/package=sigclust) [(publication)](https://www.tandfonline.com/doi/abs/10.1198/016214508000000454)
* [genefu](https://bioconductor.org/packages/release/bioc/html/genefu.html) [(publication)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6410906/)
* [ESTIMATE](https://bioinformatics.mdanderson.org/public-software/estimate/) [(publication)](https://www.nature.com/articles/ncomms3612)
* [ComplexHeatmap](http://bioconductor.org/packages/ComplexHeatmap/) [(publication)](https://academic.oup.com/bioinformatics/article/32/18/2847/1743594)
* [Biobase](https://bioconductor.org/packages/Biobase/) [(publication)](https://www.nature.com/articles/nmeth.3252)
* [GEOquery](https://www.bioconductor.org/packages/GEOquery/)[(publication)](https://academic.oup.com/bioinformatics/article/23/14/1846/190290)

These packages can be installed in R by running the following:
```r
install.packages("gtools")
install.packages("survival")
install.packages("survminer")
install.packages("circlize")
install.packages("ggplot2")
install.packages("ggsignif")
install.packages("ggsci")
install.packages("sigclust")

source("https://bioconductor.org/biocLite.R")
biocLite("genefu")

#library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap", version = "3.8")
BiocManager::install("Biobase", version = "3.8")
BiocManager::install("GEOquery", version = "3.8")
```

### To run the analyses:
1) Clone this repository
2) Download all the required data as described here and in `./Data/README.md`
3) Install the packages listed above
4) Run the following in R:

```r
setwd("/path/to/ClaudinLow/")
source("./Code/EntirePipeline.r")
```

## Output
Output from all folders is saved to the `./Output/` directory, with figures saved to `./Output/Figures/`.
