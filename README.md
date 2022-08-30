SARS-CoV-2-IHMV
===============

In this repository we provide all the code to replicate the analyses presented in https://www.cell.com/iscience/fulltext/S2589-0042(21)00084-5 

Specifically, the code is organized in 5 directories: 

*step0_data_preprocessing*, to perform variant calling and data preprocessing; 

*step1_signatures_discovery_training_set*, to perform signatures discovery on the training dataset from NCBI BioProject PRJNA645906; 

*step2_signatures_validation_testing_set*, to perform signatures assignments on the testing datasets from NCBI BioProjects PRJNA625551, PRJNA633948, PRJNA636748 and PRJNA647529; 

*step3_signatures_bootstrap*, to perform assessment of signatures significance by bootstrap both on training dataset and testing datasets; 

*step4_VERSO_step1_training_set*, to perform phylogenetic inference by VERSO step 1 on training dataset; 

*step5_MrBayes_training_set*, to perform phylogenetic inference by MrBayes on training dataset. 

*figures*, to generate the figures included in the manuscript. Move to this directory and run or source the 'generate_plots.R' script.

### License

All the software contained in this repository is distributed under Apache License. 
