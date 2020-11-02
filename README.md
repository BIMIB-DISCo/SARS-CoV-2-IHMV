SARS-CoV-2-IHMV
===============

The repository includes all the code to replicate the analyses presented in https://www.biorxiv.org/content/10.1101/2020.07.06.189944v1 

Specifically, we provide:

step0_data_preprocessing, to perform variant calling and data preprocessing; 
step1_signatures_discovery_training_set, to perform signatures discovery on the training dataset from NCBI BioProject PRJNA645906; 
step2_signatures_validation_testing_set, to perform signatures assignments on the testing datasets from NCBI BioProjects PRJNA625551, PRJNA633948, PRJNA636748 and PRJNA647529; 
step3_signatures_bootstrap, to perform assessment of signatures significance by bootstrap both on training dataset and testing datasets; 
step4_MrBayes_training_set, to perform phylogenetic inference by MyBayes on training dataset. 

### License

All the software contained in this repository is distributed under an Apache License. 
