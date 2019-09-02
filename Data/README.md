# METABRIC
* Processed data was downloaded from [cBioportal](http://www.cbioportal.org/study?id=brca_metabric) [(download link - large file)](http://download.cbioportal.org/brca_metabric.tar.gz) and saved to `./Data/brca_metabric/`. Files used in the analyses were downloaded on 28.03.2019.

* Supplementary files 2 and 3 from [Curtis et al. 2012](https://www.nature.com/articles/nature10983) are saved to `./Data/table_S2_revised.txt` and `./Data/table_S3_revised.txt`. These files can not be hosted here and must be downloaded from their source. Files used in the analyses were downloaded on 28.03.2019.

* The copy number segment file from [Pereira et al. 2016](https://www.nature.com/articles/ncomms11479) is downloaded from [here](https://raw.githubusercontent.com/cclab-brca/mutationalProfiles/master/Data/ascatSegments_withoutCNVs.txt) and saved to `./Data/ascatSegments_withoutCNVs.txt`. The copy number segment file uses different sample IDs, which can be mapped with [this file](https://raw.githubusercontent.com/cclab-brca/mutationalProfiles/master/Data/tumorIdMap.txt), saved to `./Data/tumorIdMap.txt`. These files are licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) and can therefore be re-hosted here. These file will be included if the entire repository is cloned and does not need to be dowloaded separately. The files used in the analyses were downloaded on 31.03.2019.

The folder structure for METABRIC should therefore be:
```
./Data/
├──brca_metabric/
│ 	├──case_lists/
│ 	│	├──cases_all.txt
│ 	│	├──cases_cna.txt
│ 	│	├──cases_cnaseq.txt
│ 	│	├──cases_complete.txt
│ 	│	├──cases_nat_comm_2016.txt
│ 	│	├──cases_nature_2012.txt
│ 	│	├──cases_RNA_Seq_mRNA.txt
│ 	│	└──cases_sequenced.txt
│	├──data_clinical_patient.txt
│	├──data_clinical_sample.txt
│	├──data_CNA.txt
│	├──data_expression_median.txt
│	├──data_gene_matrix.txt
│	├──data_gene_panel_metabric_173.txt
│	├──data_mRNA_median_Zscores.txt
│	├──data_mutations_extended.txt
│	├──data_mutations_mskcc.txt
│	├──LICENSE
│	├──meta_clinical_patient.txt
│	├──meta_clinical_sample.txt
│	├──meta_CNA.txt
│	├──meta_expression_median.txt
│	├──meta_gene_matrix.txt
│	├──meta_mRNA_median_Zscores.txt
│	├──meta_mutations_extended.txt
│	└──meta_study.txt
├──ascatSegments_withoutCNVs.txt
├──table_S2_revised.txt
├──table_S3_revised.txt
└──tumorIdMap.txt
```

# Oslo2
Gene expression data, hormone receptor status, and intrinsic (PAM50) subtype from the [Oslo2 cohort](https://breast-cancer-research.biomedcentral.com/articles/10.1186/s13058-017-0812-y) is downloaded from [GEO accession GSE80999](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80999) using the getGEO command in R in the `./Code/Oslo2_CoreClaudinLow.r` script; this is done automatically and the data does not need to be downloaded from the linked repository.

# TCGA
Gene expression data from [TCGA-BRCA](https://www.nature.com/articles/nature11412) is downloaded from [cBioportal](http://www.cbioportal.org/study?id=brca_tcga_pan_can_atlas_2018) - [download link - large file](http://download.cbioportal.org/brca_tcga_pan_can_atlas_2018.tar.gz). The downloaded file is saved in `./Data/` and unzipped so that all files are located under `./Data/brca_tcga_pan_can_atlas_2018`.
