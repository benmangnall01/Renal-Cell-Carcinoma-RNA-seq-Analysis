# Renal-Cell-Carcinoma-RNA-seq-Analysis

## Overview

This project analyses bulk RNA-seq data from papillary renal cell carcinoma (pRCC) and matched normal kidney tissue.

The main research aim was to identify which genes and biological pathways are consistently altered between pRCC tumor and matched normal tissue.

The study contains 106 RNA-seq samples from 53 matched tumor/normal patient pairs (GSE180777).

## Methods 

**Quality control**

Raw and trimmed FASTQ files were assessed with FastQC.

QC included:

Per-base sequence quality \
Sequence composition\
GC content\
Duplication\
Adapter content\
Read-length distribution

**Gene filtering**

Genes were retained when they had at least 10 reads in at least 10 samples. This reduced the dataset from 60,675 genes to 41,671 genes for downstream analysis.

**Differential expression**

Differential expression was performed with DESeq2 as tumor vs Normal cells

**Pathway analysis**

Biological interpretation was performed using Gene Ontology (GO) Biological Process enrichment and GO Biological Process GSEA. The GSEA used the complete ranked differential-expression result rather than only statistically significant genes.

## Results
**Differential expression**

The analysis identified widespread transcriptional differences between tumor and matched normal tissue. Among the most statistically significant genes were:

**Increased in tumor:**\
MKI67\
BUB1\
E2F8\
TPX2\
TOP2A\
CENPE\
ASPM\
DTL

These genes are associated with cell proliferation, DNA replication and mitotic processes.

**Decreased in tumor:**\
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

GO enrichment of tumor-upregulated genes was strongly dominated by immune-related biological processes, including:

Leukocyte-mediated immunity\
Lymphocyte-mediated immunity\
Regulation of lymphocyte activation\
Leukocyte proliferation\
Lymphocyte proliferation\
Adaptive immune response\
Leukocyte migration\
Regulation of T-cell activation

These results indicate a strong tumor-associated immune expression signal.

GO enrichment of downregulated genes was dominated by metabolic and renal processes, particularly small-molecule, organic-acid and amino-acid metabolism together with renal transport processes. GSEA supported these findings.

Strong negative enrichment was observed for metabolic processes including:

Organic acid catabolic process\
Carboxylic acid catabolic process\
Monocarboxylic acid catabolic process\
Fatty acid catabolic process

Positive enrichment was observed for proliferative processes including:

Pister chromatid segregation\
Chromosome segregation\
Cell-cycle checkpoint signalling\
Mitotic cell-cycle processes

For example, the GSEA normalized enrichment score (NES) was approximately -3.25 for organic acid catabolism and +3.06 for sister chromatid segregation.

**Biological interpretation**

Taken together, the results indicate three broad transcriptional patterns distinguishing pRCC tumor from adjacent-normal tissue:

Increased cell proliferation and cell-cycle activity\
Strong immune-associated transcriptional activity\
Reduced renal/metabolic and small-molecule catabolic programs

Because this is bulk RNA-seq, the immune-related signal may reflect differences in the cellular composition of tumor and normal tissue as well as changes in gene expression within individual cell types.

**Key figures**

Library-size QC\
Detected-gene QC\
PCA\
MA plot\
Volcano plot\
GO enrichment of upregulated genes\
GO enrichment of downregulated genes\
GO GSEA

## Conclusions

This analysis demonstrates an end-to-end RNA-seq workflow from raw sequencing data to biological interpretation.

The pRCC samples show coordinated transcriptional changes involving tumor proliferation, immune-associated processes and loss of renal/metabolic expression.

## Dataset

GEO accession: GSE180777
