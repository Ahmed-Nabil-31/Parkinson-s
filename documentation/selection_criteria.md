Here's a table you can include in your proposal or thesis. I've focused on the information that justifies **why these datasets were selected**, rather than every experimental detail.

| GEO Accession | Study Summary                                                                                                                                                                                                                                                                                                         | Samples Used in This Study                                                    | Tissue Used                 | Experimental Design | Why Selected                                                                                                                                                                                |
| ------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------- | --------------------------- | ------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **GSE136666** | Bulk RNA-seq analysis of postmortem human substantia nigra (SN) and putamen from Parkinson's disease (PD) patients and age-matched controls to identify differentially expressed genes associated with PD pathology. SN samples were pooled from two individuals; putamen samples were pooled from three individuals. | **Only SN samples:** PD vs age-matched healthy controls (putamen excluded)    | Postmortem substantia nigra | Bulk RNA-seq        | Human study, bulk RNA-seq, PD and control groups, same brain region (SN), gene-level count data available, directly relevant to PD transcriptomics.                                         |
| **GSE168496** | Bulk RNA-seq study investigating sex-dependent transcriptomic differences in Parkinson's disease. Although the original objective was to study sex differences, the dataset contains PD and non-demented control substantia nigra samples suitable for conventional differential expression analysis.                 | PD vs non-demented controls                                                   | Postmortem substantia nigra | Bulk RNA-seq        | Human study, bulk RNA-seq, includes non-demented controls, same brain region (SN), suitable for PD vs control analysis after ignoring sex stratification (or adjusting for sex if desired). |
| **GSE133101** | Transcriptomic profiling of multiple brain regions, including substantia nigra, middle temporal gyrus (MTG), and amygdala, from Parkinson's disease patients and controls.                                                                                                                                            | **Only substantia nigra samples:** PD vs controls (MTG and amygdala excluded) | Postmortem substantia nigra | Bulk RNA-seq        | Human study, bulk RNA-seq, contains both PD and control SN samples. Only SN samples were selected to maintain tissue consistency across datasets.                                           |

---

# Dataset Selection Criteria

The following inclusion criteria were used to ensure comparability among datasets:

| Criterion              | Requirement                                                |
| ---------------------- | ---------------------------------------------------------- |
| Species                | **Homo sapiens**                                           |
| Disease                | Parkinson's disease                                        |
| Control group          | Healthy or non-demented controls                           |
| Tissue                 | **Substantia nigra only**                                  |
| Sequencing method      | Bulk RNA sequencing                                        |
| Data availability      | Gene-level count matrix available in GREIN                 |
| Sample type            | Postmortem brain tissue                                    |
| Analysis compatibility | Suitable for differential expression analysis using DESeq2 |
| Experimental design    | Case-control comparison (PD vs control)                    |

---

# Exclusion Criteria

| Criterion                                     | Reason                                                                         |
| --------------------------------------------- | ------------------------------------------------------------------------------ |
| Single-cell or single-nucleus RNA-seq         | Not comparable with bulk RNA-seq workflow                                      |
| Brain regions other than substantia nigra     | Different baseline transcriptomes; excluded to reduce biological heterogeneity |
| Studies without healthy/non-demented controls | Cannot perform PD vs control differential expression                           |
| Non-human datasets                            | Outside the scope of this study                                                |
| Microarray datasets                           | Different technology; not directly comparable with RNA-seq                     |

---

## Why these three datasets are comparable

| Feature                             | GSE136666 | GSE168496 | GSE133101       |
| ----------------------------------- | --------- | --------- | --------------- |
| Human                               | ✓         | ✓         | ✓               |
| Parkinson's disease                 | ✓         | ✓         | ✓               |
| Healthy/non-demented controls       | ✓         | ✓         | ✓               |
| Bulk RNA-seq                        | ✓         | ✓         | ✓               |
| Gene-level counts                   | ✓         | ✓         | ✓               |
| Postmortem tissue                   | ✓         | ✓         | ✓               |
| **Substantia nigra analyzed**       | ✓         | ✓         | ✓ (subset only) |
| Suitable for recurrent DEG analysis | ✓         | ✓         | ✓               |

### Minor considerations to mention in your Methods

To be transparent, you can include a short note:

* **GSE136666:** Only substantia nigra samples were included; putamen samples were excluded. Substantia nigra samples were generated from pooled tissue of two individuals.
* **GSE168496:** Although the original study focused on sex-specific transcriptomic differences, only the Parkinson's disease versus non-demented control comparison was analyzed in this study.
* **GSE133101:** Only substantia nigra samples were included. Samples from the middle temporal gyrus and amygdala were excluded to maintain tissue consistency across datasets.

This level of documentation is typically sufficient to justify dataset selection in a thesis or manuscript.
