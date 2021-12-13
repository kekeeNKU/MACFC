# Multi-variable AUC for Sifting Complementary Features and Its Biomedical Application

### Data
1. UCI data and multi-class data were downloaded from [UC Irvine Machine Learning Repository](https://archive.ics.uci.edu/ml/index.php)
2. TCGA datasets are preprocessed data from the Xena platform in the “gene expression RNAseq -IlluminaHiSeq pancan normalized” version. The file is too large to upload on Github. Please see details and download datasets on the provided web links.
* [BRCA (TCGA Breast Cancer)](https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap%2FHiSeqV2_PANCAN&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
* [KIRC (TCGA Kidney Clear Cell Carcinoma)](https://xenabrowser.net/datapages/?dataset=TCGA.KIRC.sampleMap%2FHiSeqV2_PANCAN&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
* [LIHC (TCGA Liver Cancer)](https://xenabrowser.net/datapages/?dataset=TCGA.LIHC.sampleMap%2FHiSeqV2_PANCAN&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
* [LUAD (TCGA Lung Adenocarcinoma)](https://xenabrowser.net/datapages/?dataset=TCGA.LUAD.sampleMap%2FHiSeqV2_PANCAN&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
* [PRAD (TCGA Prostate Cancer)](https://xenabrowser.net/datapages/?dataset=TCGA.PRAD.sampleMap%2FHiSeqV2_PANCAN&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)
* [THCA (TCGA Thyroid Cancer)](https://xenabrowser.net/datapages/?dataset=TCGA.THCA.sampleMap%2FHiSeqV2_PANCAN&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443)

### Code
1. data_process.py is code for data preprocessing.
2. feature_selection.py is code for sifting an optimal feature set.
3. feature_ranking.py is code for ranking features from high to low and also sorting features with their frequency.
 
### Contributions of this algorithm:
1. Globally evaluate the complementarity among features.
2. Screen discriminative combination of features that are complementary to each other from a global view.

Note: features mean genes when applied to gene expression datasets.
