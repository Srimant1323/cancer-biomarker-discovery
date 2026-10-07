# Cancer Biomarker Discovery from Gene Expression Data (Lung Cancer)

[!\[Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YOUR_USERNAME/cancer-biomarker-discovery/blob/main/notebooks/cancer_biomarker_discovery.ipynb)

Identification and independent validation of candidate lung cancer biomarkers from public microarray data, using Python (pandas, scikit-learn, seaborn, scipy, statsmodels).

## Key results

* **Data:** 60 paired tumor / adjacent-normal lung samples (120 arrays) from non-smoking women, GEO accession [GSE19804](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19804). 16,241 genes after quality control.
* **Differential expression (paired t-test, FDR < 0.05, |log2FC| > 1):** 817 genes, of which **271 are upregulated** and **546 are downregulated** in tumors.
* **Biology check:** reproduced the downregulation of **SEMA5A**, the biomarker reported by the original study (log2FC = -1.56, FDR = 3.6e-19).
* **Classifier:** a 20-gene logistic regression reached a **cross-validated AUC of 0.983** (patient-grouped 5-fold CV, gene selection inside each fold).
* **Independent validation:** trained on GSE19804 and tested without retraining on [GSE18842](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18842) (91 samples): **AUC = 0.9995**.
* **Replication:** 98.0% of the 817 significant genes changed in the same direction in the independent cohort.

## Dataset

|||
|-|-|
|Discovery cohort|GSE19804, 60 patients x 2 tissues (tumor, adjacent normal), Affymetrix GPL570|
|Validation cohort|GSE18842, 46 tumor and 45 control lung samples, Affymetrix GPL570|
|Data type|Microarray, already log2-scaled and normalized (values 3.04 to 14.89 in GSE19804)|
|Original study|Lu et al., Cancer Epidemiol Biomarkers Prev 2010 (PMID 20802022)|

## Methods

1. **Quality control.** No missing values. Per-sample boxplots showed identical distributions, consistent with quantile normalization (RMA), so no further normalization was applied.
2. **Cleaning (probes to genes).** 54,675 probes, then 42,986 probes mapped to a single gene symbol, then 21,655 unique genes after averaging duplicate probes, then 16,241 genes after removing the lowest-variance quartile (a label-free filter).
3. **PCA** on standardized expression of the 16,241 genes.
4. **Differential expression.** Paired t-test per gene (each tumor compared with the same patient's normal tissue), Benjamini-Hochberg FDR correction. Because the data are already log2, log2FC is the mean paired difference. Thresholds: FDR < 0.05 and |log2FC| > 1 (at least 2-fold).
5. **Classifier.** Standardization, top 20 genes by ANOVA F-test, logistic regression, all inside a scikit-learn pipeline so gene selection is repeated within each training fold. `GroupKFold` keeps both samples of a patient in the same fold to prevent leakage.
6. **External validation.** The classifier was fit on GSE19804 only and applied once to GSE18842 using the 16,241 shared genes, with and without per-gene z-scoring within each dataset (no labels used) to correct dataset-level offsets.
7. **Pathway enrichment.** GO Biological Process enrichment (Enrichr via `gseapy`) for the up and down gene lists.

## Results

### Sample structure (PCA)

PC1 and PC2 explain 18.9% and 12.5% of the variance. Tumor and normal samples separate partially along a consistent diagonal direction; connecting each patient's two samples shows that almost every patient shifts the same way from normal to tumor.

!\[PCA](results/pca\_plot.png)
!\[PCA paired](results/pca\_paired.png)

Eleven samples lie far from the main cloud. They come largely from the same patients (131, 132, 135, 137), with both tumor and normal samples shifted together, which suggests a patient-specific or batch effect. The GEO metadata (stage, gender, protocols) did not explain it. No samples were removed; the paired test compares each tumor with its own normal tissue.

### Differential expression

!\[Volcano](results/volcano\_plot.png)
!\[MA plot](results/ma\_plot.png)

|Status|Genes|
|-|-|
|Upregulated in tumor|271|
|Downregulated in tumor|546|
|Not significant|15,424|

### Top 10 upregulated genes (ranked by log2FC among significant genes)

|Gene|log2FC|FDR|
|-|-|-|
|COL10A1|4.03|1.7e-20|
|SPINK1|3.34|3.5e-10|
|CTHRC1|3.19|1.4e-18|
|MMP12|3.18|2.4e-13|
|MMP1|3.04|2.2e-10|
|CST1|2.99|2.7e-15|
|COL11A1|2.97|1.3e-13|
|GREM1|2.87|4.8e-12|
|SPP1|2.78|9.5e-18|
|HS6ST2|2.70|1.9e-15|

### Top 10 downregulated genes (ranked by log2FC among significant genes)

|Gene|log2FC|FDR|
|-|-|-|
|AGER|-3.93|1.3e-22|
|WIF1|-3.77|3.7e-13|
|TMEM100|-3.57|1.1e-15|
|FCN3|-3.41|4.0e-16|
|SOSTDC1|-3.39|1.0e-16|
|SCGB1A1|-3.38|2.6e-10|
|CLDN18|-3.25|4.4e-15|
|GKN2|-3.25|1.2e-13|
|CA4|-3.19|4.0e-23|
|SFTPC|-3.14|1.7e-11|

A second table ranked by statistical significance instead of effect size is in `results/biomarker\_table\_by\_FDR.csv`. Four genes appear in the top 10 of both rankings (AGER, CA4, COL10A1, CTHRC1). The full results for all 16,241 genes are in `results/full\_DE\_results.csv`.

!\[Heatmap](results/heatmap.png)

The heatmap of the top 20 genes is illustrative: those genes were selected because they differ between groups, so the clean separation is partly circular. The classifier and the external cohort below are the independent evidence.

### Known biomarker: SEMA5A

The original study of this cohort reported downregulation of SEMA5A in tumors. This analysis reproduces it (log2FC = -1.56, FDR = 3.6e-19). SEMA5A ranks 61st of 16,241 genes by FDR but 200th of 546 downregulated genes by fold change, so it does not appear in the top-10 effect-size table. It also replicates in the independent cohort (log2FC = -1.51, FDR = 2.7e-12).

### Classifier and validation

|Setting|AUC|Accuracy|
|-|-|-|
|GSE19804, 5-fold patient-grouped CV|0.983|95.8% (115/120)|
|GSE18842, no harmonization|0.9995|93.4% (85/91)|
|GSE18842, per-gene z-score|0.9995|97.8% (89/91)|

Within each patient, the tumor scored higher than the paired normal in 59 of 60 patients. All three tumors missed in cross-validation were stage 1A, which suggests early tumors are harder to detect, though this rests on only three samples.

!\[ROC CV](results/roc\_curve.png)
!\[ROC validation](results/roc\_validation.png)

The ranking quality (AUC) is unaffected by the offset between the two datasets, but accuracy at a fixed 0.5 cutoff improved with per-gene standardization.

### Replication of gene-level results

Among the 817 significant discovery genes, 98.0% changed in the same direction in GSE18842. The Spearman correlation of log2FC across all genes was 0.73.

### Pathway enrichment (GO Biological Process)

!\[Pathways](results/pathway\_enrichment.png)

* **Upregulated genes** are enriched for extracellular matrix organization and for mitotic chromosome segregation and spindle processes.
* **Downregulated genes** are enriched for inflammatory response, regulation of angiogenesis, and regulation of the ERK/MAPK cascade.

## Limitations

* Tumor versus adjacent normal lung is a high-contrast comparison. A near-perfect AUC shows the signature is robust across two studies, not that it is a clinical diagnostic test.
* Both cohorts used the same platform (GPL570), so cross-platform robustness was not tested. Both are from tissue, not blood, so this says nothing about non-invasive biomarkers.
* Many downregulated genes (alveolar, vascular, immune markers) likely reflect differences in cell composition between tumor and normal tissue as well as gene regulation within cells.
* The classifier's coefficients should not be read as per-gene effect directions because the selected genes are correlated (for example CA4 has a positive coefficient despite being downregulated in tumors). Gene direction is taken from the differential expression results.
* Validation differential expression used an unpaired test (GSE18842 is not fully paired), and the per-gene z-score used the test cohort's own distribution, which would need a different normalization approach for single new patients.
* Enrichment used Enrichr's default gene background rather than the 16,241 measured genes. Many top GO terms are redundant.
* The cohort is single-sex (female), non-smoking, from Taiwan, and mostly early stage (35 of 60 tumors are stage I).
* The cause of the PCA outliers could not be identified from the GEO metadata.

## Repository structure

```
cancer-biomarker-discovery/
├── README.md
├── requirements.txt
├── .gitignore
├── notebooks/
│   └── cancer\_biomarker\_discovery.ipynb
├── data/
│   └── sample\_labels.csv
└── results/
    ├── qc\_boxplot.png
    ├── pca\_plot.png
    ├── pca\_paired.png
    ├── volcano\_plot.png
    ├── ma\_plot.png
    ├── heatmap.png
    ├── roc\_curve.png
    ├── roc\_validation.png
    ├── pathway\_enrichment.png
    ├── biomarker\_table.csv
    ├── biomarker\_table\_by\_FDR.csv
    ├── classifier\_genes.csv
    └── full\_DE\_results.csv
```

## How to run

1. Click the **Open in Colab** badge above.
2. Choose **Runtime > Run all**. The notebook downloads both datasets from GEO automatically (internet access required).

To run locally: `pip install -r requirements.txt`, then open the notebook with Jupyter.

## References

* Lu TP, Tsai MH, Lee JM, et al. Identification of a novel biomarker, SEMA5A, for non-small cell lung carcinoma in nonsmoking women. *Cancer Epidemiol Biomarkers Prev* 2010;19(10):2590-7. PMID 20802022.
* GEO series GSE19804 (discovery) and GSE18842 (validation).

## Author

Srimant Bhardwaj 



