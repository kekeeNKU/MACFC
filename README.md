# Multi-variable AUC for Sifting Complementary Features and Its Biomedical Application

Data
1. UCI data and multi-class data downloaded from https://archive.ics.uci.edu/ml/index.php
2. TCGA data downloaded from https://xenabrowser.net/datapages/ in the “gene expression RNAseq -IlluminaHiSeq pancan normalized” version. The file is too large to upload on Github, please see details on the provided website.

Code
1. data_process.py is code for data preprocessing.
2. feature_selection.py is code for sifing an optimal feature set.
3. feature_ranking.py is code for ranking features from high to low, and also sorting features with its frequency.
 
Contributions for this paper:
1. Globally evaluate the complementarity among features.
2. Screen discriminative combination of features that are complementary to each other from global view by using MACFC.

Note: features means genes when apply to gene expression datasets.
