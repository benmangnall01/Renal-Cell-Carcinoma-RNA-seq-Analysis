# Renal-Cell-Carcinoma-RNA-seq-Analysis

## Overview

This project analyses bulk RNA-seq data from papillary renal cell carcinoma (pRCC) and matched adjacent-normal kidney tissue.

The main research question was:

Which genes and biological pathways are consistently altered between pRCC tumour and matched normal tissue?

The study contains 106 RNA-seq samples from 53 matched tumour/normal patient pairs (GSE180777).

The analysis was performed primarily in R and covers the workflow from raw sequencing reads through differential expression and pathway analysis.

**Methods** 

**Quality control**

Raw and trimmed FASTQ files were assessed with FastQC.

QC included:

per-base sequence quality \
sequence composition\
GC content\
duplication\
adapter content\
read-length distribution

**Gene filtering**

Genes were retained when they had at least 10 reads in at least 10 samples.

This reduced the dataset from 60,675 genes to 41,671 genes for downstream analysis.

**Differential expression**

Differential expression was performed with DESeq2 as Tumour vs Normal cells

**Pathway analysis**

Biological interpretation was performed using:

Gene Ontology (GO) Biological Process enrichment\
GO Biological Process GSEA

The GSEA used the complete ranked differential-expression result rather than only statistically significant genes.

**Results**\
**Differential expression**

The analysis identified widespread transcriptional differences between tumour and matched normal tissue.

Among the most statistically significant genes were:

**Increased in tumour**\
MKI67\
BUB1\
E2F8\
TPX2\
TOP2A\
CENPE\
ASPM\
DTL

These genes are associated with cell proliferation, DNA replication and mitotic processes.

**Decreased in tumour**\
ALDOB\
F11\
UMOD\
AGMAT\
PCK1\
SLC9A3\
G6PC\
PLPPR1

These changes contribute to a strong metabolic/renal-tissue expression signature.

The strongest individual changes included large positive and negative log2 fold changes, with many genes reaching extremely small adjusted p-values.

**Pathway analysis**

GO enrichment of tumour-upregulated genes was strongly dominated by immune-related biological processes, including:

leukocyte-mediated immunity\
lymphocyte-mediated immunity\
regulation of lymphocyte activation\
leukocyte proliferation\
lymphocyte proliferation\
adaptive immune response\
leukocyte migration\
regulation of T-cell activation

These results indicate a strong tumour-associated immune expression signal.

GO enrichment of downregulated genes was dominated by metabolic and renal processes, particularly small-molecule, organic-acid and amino-acid metabolism together with renal transport processes.

GSEA supported these findings.

Strong negative enrichment was observed for metabolic processes including:

organic acid catabolic process\
carboxylic acid catabolic process\
monocarboxylic acid catabolic process\
fatty acid catabolic process

while positive enrichment was observed for proliferative processes including:

sister chromatid segregation\
chromosome segregation\
cell-cycle checkpoint signalling\
mitotic cell-cycle processes

For example, the GSEA normalized enrichment score (NES) was approximately -3.25 for organic acid catabolism and +3.06 for sister chromatid segregation.

**Biological interpretation**

Taken together, the results indicate three broad transcriptional patterns distinguishing pRCC tumour from adjacent-normal tissue:

Increased cell proliferation and cell-cycle activity\
Strong immune-associated transcriptional activity\
Reduced renal/metabolic and small-molecule catabolic programs

Because this is bulk RNA-seq, the immune-related signal may reflect differences in the cellular composition of tumour and normal tissue as well as changes in gene expression within individual cell types.

**Key figures**

The main outputs include:

Library-size QC\
Detected-gene QC\
PCA\
MA plot\
Volcano plot\
GO enrichment of upregulated genes\
GO enrichment of downregulated genes\
GO GSEA

**Conclusions**

This analysis demonstrates an end-to-end RNA-seq workflow from raw sequencing data to biological interpretation.

The pRCC samples show coordinated transcriptional changes involving tumour proliferation, immune-associated processes and loss of renal/metabolic expression programs.

The project also demonstrates independent processing of raw RNA-seq reads through quality control, trimming, alignment and gene-level quantification rather than relying exclusively on a deposited processed dataset.

**Dataset**

GEO accession: GSE180777
