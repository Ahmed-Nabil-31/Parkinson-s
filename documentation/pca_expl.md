Yes. In fact, you can write this in a way that is both scientifically accurate and defensible. I would avoid claiming that the lack of clustering is "good" or "bad." Instead, explain what PCA is designed to show and interpret your findings objectively.

---

## Results

**Principal Component Analysis (PCA)**

Principal component analysis (PCA) was performed on variance-stabilized gene expression data to assess the overall transcriptomic variation among samples and to identify potential outliers or technical batch effects prior to downstream differential expression analysis. Across all three Parkinson's disease (PD) datasets (GSE133101, GSE136666, and GSE168496), the first two principal components explained between **21–31%** and **7–15%** of the total variance, respectively. Although individual samples exhibited considerable dispersion, PD and control samples did not form distinct clusters along the first two principal components.

No obvious clustering attributable to technical artifacts or sample processing was observed, and no extreme outliers were identified. The overlap between PD and control samples suggests that disease status was not the dominant contributor to global transcriptomic variance within these postmortem substantia nigra samples.

---

## Discussion / Interpretation

The absence of complete separation between PD and control samples in PCA is not unexpected for bulk RNA-seq analyses of human postmortem brain tissue. Principal component analysis is an **unsupervised dimensionality reduction technique** that identifies the largest sources of variance across the entire transcriptome rather than maximizing separation between predefined experimental groups. Consequently, disease-associated transcriptional changes may not necessarily dominate the first principal components if other biological sources of variation contribute more strongly to overall variance. ([arXiv][1])

Several factors likely contribute to the observed overlap between PD and control samples. Human substantia nigra tissue exhibits substantial biological heterogeneity arising from inter-individual genetic variation, age, sex, disease duration, medication history, postmortem interval (PMI), RNA integrity, and cellular composition. These factors frequently contribute more to global transcriptomic variability than disease status alone, particularly in neurodegenerative disorders where only a subset of genes is consistently dysregulated. ([PubMed Central (PMC)][2])

Furthermore, PCA summarizes variance across all expressed genes simultaneously. While Parkinson's disease is characterized by differential expression of specific molecular pathways—including mitochondrial dysfunction, protein folding, neuroinflammation, and dopaminergic signaling—these disease-associated genes represent only a fraction of the transcriptome. Consequently, significant differential expression can be detected by statistical models such as DESeq2 even when the overall transcriptomic profiles of PD and control samples partially overlap in PCA space. ([ScienceDirect][3])

Because no obvious clustering indicative of technical batch effects was observed, the PCA results suggest that the major sources of variation are primarily biological rather than technical. Therefore, each dataset was analyzed independently using DESeq2 to identify differentially expressed genes, followed by cross-study comparison to identify recurrent PD-associated transcriptional signatures.

---

## Why this interpretation is valid

There are three key reasons this explanation is supported by the literature:

1. **PCA is a quality control tool, not a statistical test for differential expression.** Its primary purpose in RNA-seq workflows is to visualize major sources of variance, detect outliers, and identify possible batch effects rather than to determine whether two biological groups are significantly different. ([PubMed Central (PMC)][2])

2. **Human Parkinson's disease is highly heterogeneous.** Even studies with more than 100 postmortem RNA-seq samples report substantial molecular heterogeneity within PD, to the extent that patients can be partitioned into distinct molecular subtypes. This indicates that variability among PD patients themselves is considerable and can obscure simple PD-versus-control separation in PCA. ([PubMed Central (PMC)][4])

3. **Differential expression does not require PCA separation.** DESeq2 performs gene-wise statistical modeling, testing each gene individually while accounting for within-group variability. Therefore, hundreds of significant DEGs may still be identified despite overlapping PCA distributions, provided that consistent gene-specific expression differences exist between conditions. ([PubMed Central (PMC)][2])

---

### One suggestion for your manuscript

Since you **did not observe evidence of a strong technical batch effect**, I would **remove the post-correction PCA figures** and simply present **one PCA per dataset** (before any SVA). Because you ultimately decided that SVA was unnecessary, showing only the standard PCA avoids implying that batch correction materially changed the data when it did not. This is also the approach most RNA-seq papers take unless batch correction substantially alters the sample distribution.

[1]: https://arxiv.org/abs/1404.1100?utm_source=chatgpt.com "A Tutorial on Principal Component Analysis"
[2]: https://pmc.ncbi.nlm.nih.gov/articles/PMC9284682/?utm_source=chatgpt.com "Batch effect detection and correction in RNA-seq data using machine-learning-based automated assessment of quality - PMC"
[3]: https://www.sciencedirect.com/science/article/abs/pii/S0304394026001060?utm_source=chatgpt.com "Dysregulation of mRNAs and hub genes in Parkinson’s disease within post mortem substantia Nigra: using three methods differential expression genes analysis - ScienceDirect"
[4]: https://pmc.ncbi.nlm.nih.gov/articles/PMC8910747/?utm_source=chatgpt.com "A Hybrid Machine Learning and Network Analysis Approach Reveals Two Parkinson’s Disease Subtypes from 115 RNA-Seq Post-Mortem Brain Samples - PMC"
