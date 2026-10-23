That's a really good question, and it's one of the most common misconceptions about PCA in RNA-seq.

The important point is this:

> **The PCA was *colored* by PD and control, but it was *not calculated* using PD and control.**

These are two different things.

Here's what happens internally:

```r
vsd <- vst(dds)
pca <- plotPCA(vsd, intgroup = "condition")
```

`plotPCA()` first computes PCA **using only the expression matrix** (all genes × all samples).

It does **not** use the `condition` variable when calculating the principal components.

Only **after** PCA is computed does it use

```r
intgroup = "condition"
```

to color the points red and blue.

So mathematically, PCA is asking:

> "What are the directions of greatest variance among all samples?"

It is **not** asking:

> "Can I best separate PD from controls?"

That second question would require a **supervised** method such as Linear Discriminant Analysis (LDA), PLS-DA, or a classifier.

---

### Why don't PD and controls separate?

For them to separate automatically in PCA, **disease status must explain one of the largest sources of variance**.

Imagine the total variance comes from

* age
* sex
* neuronal loss
* glial infiltration
* postmortem interval
* RNA integrity
* medication
* individual genetics
* **PD**

If PD contributes only ~8% of total variance while age contributes 20%, then PC1 will mostly represent age, not PD.

Yet DESeq2 can still identify genes consistently altered by PD.

---

### A classic example

The original DESeq2 paper (Love et al., 2014) explicitly uses PCA as a **quality assessment tool**, not as evidence of differential expression. They note that PCA is intended to detect outliers and dominant sources of variation, while differential expression is determined by gene-wise statistical modeling rather than overall sample clustering.

Similarly, the RNA-seq workflow by Love, Anders, Huber, and colleagues explains that PCA is an **unsupervised exploratory analysis** whose principal components represent the directions of greatest variance, regardless of the experimental groups.

---

### What would be concerning?

Suppose your PCA looked like this:

```
PC1

Dataset 1: ●●●●●●

Dataset 2: ▲▲▲▲▲▲

Dataset 3: ■■■■■■
```

instead of

```
PC1

PD     ● ● ● ● ●
Ctrl   ○ ○ ○ ○ ○
```

Then you'd conclude that **dataset/batch**, not biology, dominates the variation.

Your plots do **not** show that. Instead, they show intermingled PD and control samples without a clear technical grouping.

---

### So in your manuscript, avoid saying:

> "PCA showed no biological differences between PD and controls."

That would be incorrect.

Instead, write:

> "PCA did not reveal complete separation between PD and control samples along the first two principal components. Because PCA is an unsupervised method that captures the dominant sources of overall transcriptomic variation rather than maximizing separation between predefined experimental groups, this overlap indicates that disease status was not the principal determinant of global variance. Nevertheless, disease-associated transcriptional differences can still be identified through gene-wise statistical modeling using DESeq2."

That interpretation is fully consistent with how PCA is intended to be used in bulk RNA-seq analyses.
