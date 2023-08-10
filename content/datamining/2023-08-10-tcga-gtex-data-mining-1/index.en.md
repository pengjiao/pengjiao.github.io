---
title: TCGA GTEx data mining 1 1
author: Pengjiao
date: '2023-08-09'
slug: []
categories:
  - R
tags:
  - R Markdown
weight: 5
output: 
  html_document:
    toc: yes
    toc_float: yes
    highlight: tango
---




The FOLR1 gene and the MUC16 gene are significant genes involved in ovarian cancer, potentially serving as therapeutic targets.

In this blog post, I will investigate the expression levels of the MUC16 gene in ovarian cancer patients who exhibit low expression of the FOLR1 gene. The analysis will involve comparing gene exression data from ovarian cancer tissue with that from normal ovary tissue. However, given the lack of normal ovary control in The Cancer Genome Atlas (TCGA) data, I have incorporated normal ovary data from the Genotype-Tissue Expression (GTEx) database.

The RECOUNT2 and UCSC Xena projects, which provide gene expression data from TCGA and GTEx database, utilize consistent preprocessing pipelines to reduce batch effects and enhance the reliability of comparisons across diverse studies. Yet, these two projects may employ different preprocessing and normalization methods. To ensure the robustness of results, I have treated the two projects separately in this analysis.

For data acquisition, the R packages TCGAbiolinks (using the TCGAquery_recount2 function) and UCSCXenaTools were used to query and download TCGA-OV data and normal ovarian tissue data from GTEx, available in both the RECOUNT2 and UCSC Xena projects.


## 01. Load library

We first need to load the necessary R packages that will be utilized throughout our study.


```r
library(SummarizedExperiment)
library(TCGAbiolinks)
library(recount)
library(biomaRt)
library(org.Hs.eg.db)
library(dplyr)
library(tibble)
library(CCSBUtils)
library(ggpubr)
library(stringr)
```


## 02. Data from RECOUNT2 project

### Data prepare

Next, we obtain and prepare our gene expression data from the RECOUNT2 project using R package TCGAbiolinks (TCGAquery_recount2). 

TCGAbiolinks is an R package that provides an integrated interface to access, analyze, and visualize TCGA data. With its 'TCGAquery_recount2' function, it allows us to access the preprocessed gene expression data available on the RECOUNT2 project with ease.

In the following steps, we'll focus on tidying up this data, ensuring it's ready for the subsequent analysis stages. We'll select our genes of interest (FOLR1 and MUC16), and create subsets of the data corresponding to ovarian cancer samples and normal ovary tissue samples.

After we obtain the data, we extract the gene expression data for the ovary tissue from both projects into SummarizedExperiment (SE) objects. These objects are then saved in an RData file for future reference and reanalysis.


```r
## query from Recount2 platform: Ovary Carcinoma
recount.gtex <- TCGAquery_recount2(project="GTEX", tissue="ovary")
recount.tcga <- TCGAquery_recount2(project="TCGA", tissue="ovary")

##  save data for reanalysis
# save(recount.gtex, recount.tcga, file = "recount_gtex_tcga_ovary.RData")
#load("recount_gtex_tcga_ovary.RData")

## to get the SE object
SE.recount.gtex <- recount.gtex$GTEX_ovary
SE.recount.tcga <- recount.tcga$TCGA_ovary
```


The Recount2 project uses Rail-RNA as its primary computational pipeline to process and align RNA sequencing (RNA-seq) data. Rail-RNA is an alignment and transcriptome reconstruction tool. There is difference between coverage counts from Recount2 project using Rail-RNA and typical read count matrices. Rail-RNA employs a "binning" procedure to the aligned reads, resulting in what are sometimes called "read coverage" or "coverage counts". These coverage counts are the number of reads that overlap with each base pair of a gene, rather than just counting each read once as is done in traditional "read counts". However, most methods were developed for read typical count matrices. Therefore, the Recount2 data must be scaled in order to compare samples within and across studies. Scaling ensures that differences in sequencing depth across samples do not distort the counts, making them directly comparable. In the context of your analysis, scaling the Recount2 data ensures that the gene expression levels you observe are reflective of true biological differences, rather than technical differences from sequencing.


There are two methods to scale the count data: either auc or mapped_reads. If set to auc it will scale the counts by the total coverage of the sample. That is, the area under the curve (AUC) of the coverage. If set to mapped_reads it will scale the counts by the number of mapped reads, whether the library was paired-end or not, and the desired read length (L).

\[
\sum_{i} \frac{coverage_i}{AUC} * target = \text{{scaled read counts}} \tag{1}
\]

\[
\sum_{i} \frac{coverage_i}{\text{{Read Length}} * target}{mapped} = \text{{scaled read counts}} \tag{2}
\]



```r
## preparing/scaling Recount2 data because it was sequenced using Rail-RNA
eset.gtex <- assays(scale_counts(recount.gtex$GTEX_ovary, round = TRUE))$counts
eset.tcga <- assays(scale_counts(recount.tcga$TCGA_ovary, round = TRUE))$counts

## Check that the number of reads is less than or equal to 40 million
rse_scaled <- scale_counts(recount.gtex$GTEX_ovary, round = TRUE)
summary(colSums(assays(rse_scaled)$counts)) / 1e6
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   31.33   34.77   35.37   35.52   36.42   38.67
```


The next step is to segregate the TCGA data into primary tumors and normal samples. Note that TCGA does not typically include normal tissue samples. We then merge the GTEx and TCGA cancer data, joining them by gene names. This creates a new dataframe where each row represents a gene and the columns contain expression data from the two different sources. Duplicates, which might have same gene ID are removed. Lastly, we normalize and filter the data using TCGAanalyze functions from the TCGAbiolinks package. This includes GC content normalization, which adjusts for the varying guanine and cytosine content across different genes. Then, we apply a quantile filter, removing genes that have expression values in the lowest 25% quantile, as they may represent noise.

This completes the data preparation process. With this, our gene expression data from RECOUNT2 is now ready for further analysis.


```r
## replacing UUIDs with TCGA barcodes:
colnames(eset.tcga) <- colData(recount.tcga$TCGA_ovary)$gdc_cases.samples.portions.analytes.aliquots.submitter_id

## segregate between primary tumors and normal samples
eset.tcga.cancer <- eset.tcga[,which(colData(recount.tcga$TCGA_ovary)$gdc_cases.samples.sample_type=="Primary Tumor")]
eset.tcga.normal <- eset.tcga[,which(colData(recount.tcga$TCGA_ovary)$gdc_cases.samples.sample_type=="Solid Tissue Normal")] # there is no tissue normal

## merging data by row names
dataPrep.ov <- merge(as.data.frame(eset.gtex), as.data.frame(eset.tcga.cancer), by=0, all=TRUE)

## rename gene id and remove duplicated rows with same gene id
dataPrep.ov <- column_to_rownames(dataPrep.ov, var = "Row.names")
newgene <- str_split(rownames(dataPrep.ov), "\\.", simplify = T)[,1]
dataPrep.fil <- dataPrep.ov[which(!duplicated(newgene)),]
rownames(dataPrep.fil) <- newgene[which(!duplicated(newgene))]

## normalization and filtering data
dataNorm.ov <- TCGAanalyze_Normalization(tabDF = dataPrep.fil,
                                          geneInfo = geneInfoHT,
                                          method = "gcContent")

dataFilt.ov <- TCGAanalyze_Filtering(tabDF = dataNorm.ov,
                                      method = "quantile", 
                                      qnt.cut =  0.25)
```

### Differential analysis

Then, we perform differential expression analysis between normal and tumor samples using the TCGAanalyze_DEA function from the TCGAbiolinks package. This analysis is performed on a set of 108 normal and 430 tumor samples across 42755 mRNA or genes.

The 'method' argument is set to 'glmLRT', which specifies the use of the likelihood ratio test for model fitting in the limma package.

After the differential expression analysis, the Ensembl gene IDs in the result are converted to HGNC symbols using the BioMart package. This makes the results more interpretable as HGNC symbols are more commonly used in the literature. The 'getBM' function from the BioMart package is used for this conversion. The resulting table contains both the Ensembl gene IDs and corresponding HGNC symbols for each differentially expressed gene.



```r
## set "metadata" argument to true when dealing with TCGA data
#in order to extract batch correction data
#processing 108 normal and 430 tumor samples on 42755 mRNA or genes
DEG.ov <- TCGAanalyze_DEA( mat1 = dataFilt.ov[,colnames(eset.gtex)],
                            mat2 = dataFilt.ov[,colnames(eset.tcga.cancer)],
                            metadata =FALSE,
                            pipeline="limma",
                            voom = TRUE,
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.01 ,
                            logFC.cut = 1,
                            method = "glmLRT")
```

<img src="{{< blogdown/postref >}}index.en_files/figure-html/unnamed-chunk-5-1.png" width="672" />

```r
## converting ensenmbl gene ids to huugo gymbols using Biomart package
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), values=rownames(DEG.ov), mart= mart)

DEG.ov <- rownames_to_column(DEG.ov, var = "ensembl_gene_id")
deg.ov <- merge(DEG.ov, G_list, by = "ensembl_gene_id")
```

### Examine the expression levels

We define the low expression of the FOLR1 gene as expression level of FOLR1 gene was lower than median expression level. First we tranform counts into TPMs using getTPM function.

\[
RPKM =\frac{Number\ of\ Reads}{(\frac{Gene\ Length}{1000})*(\frac{Total\ reads}{10^6})}
\]

\[
TPM = \frac{{RPKM}}{{\sum{RPKM}}} \times 10^6
\]



```r
## transform counts into TPMs
gtex.tpm <- getTPM(scale_counts(recount.gtex$GTEX_ovary, round = TRUE)) %>% data.frame()
tcga.tpm <- getTPM(scale_counts(recount.tcga$TCGA_ovary, round = TRUE)) %>% data.frame()

## replacing UUIDs with TCGA barcodes:
colnames(tcga.tpm) <- colData(recount.tcga$TCGA_ovary)$gdc_cases.samples.portions.analytes.aliquots.submitter_id

## rename the gene id and remove duplicated rows with same gene id from TPM data
tpmgenes <- str_split(rownames(tcga.tpm), "\\.", simplify = T)[,1]
tcga.tpm.df <- tcga.tpm[which(!duplicated(tpmgenes)),]
rownames(tcga.tpm.df) <- tpmgenes[which(!duplicated(tpmgenes))]
gtex.tpm.df <- gtex.tpm[which(!duplicated(tpmgenes)),]
rownames(gtex.tpm.df) <- tpmgenes[which(!duplicated(tpmgenes))]


## converting ensenmbl gene ids to huugo gymbols using Biomart package
G_list.tpm <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"), values=rownames(tcga.tpm.df), mart= mart)

## remove the rows without gene names
### merging data by ensembl_gene_id
tcga.tpm.df$ensembl_gene_id <- rownames(tcga.tpm.df)
tcga.tpm.mer <- merge(G_list.tpm, tcga.tpm.df, by = "ensembl_gene_id")
tcga.tpm.sel <- tcga.tpm.mer[tcga.tpm.mer$hgnc_symbol != "", -1]
tcga.tpm.sel <- tcga.tpm.sel[which(!duplicated(tcga.tpm.sel$hgnc_symbol)),]
rownames(tcga.tpm.sel) <- NULL
tcga.tpm.fil <- column_to_rownames(tcga.tpm.sel, "hgnc_symbol")

gtex.tpm.df$ensembl_gene_id <- row.names.data.frame(gtex.tpm.df)
gtex.tpm.mer <- merge(G_list.tpm, gtex.tpm.df, by = "ensembl_gene_id")
gtex.tpm.sel <- gtex.tpm.mer[gtex.tpm.mer$hgnc_symbol != "",-1]
gtex.tpm.sel <- gtex.tpm.sel[which(!duplicated(gtex.tpm.sel$hgnc_symbol)),]
rownames(gtex.tpm.sel) <- NULL
gtex.tpm.fil <- column_to_rownames(gtex.tpm.sel, "hgnc_symbol")
```


There is no significant difference between low and high expressions of the FOLR1 gene. And MUC16 expression is significantly higher in patients with both low and high FOLR1 gene expression compared to normal ovary tissues.


```r
## calculate the mean tpm per gene
tcga.tpm.mean <- apply(tcga.tpm.fil, 1, mean, na.rm = T)
tcga.tpm.mean[names(tcga.tpm.mean) == "FOLR1"]
```

```
##    FOLR1 
## 419.3596
```

```r
tcga.tpm.mean[names(tcga.tpm.mean) == "MUC16"]
```

```
##    MUC16 
## 42.91637
```

```r
## extract patients with FOLR1 expression level lower than median value
lower.sample <- colnames(tcga.tpm.fil)[which(tcga.tpm.fil[rownames(tcga.tpm.fil) == "FOLR1",] < tcga.tpm.mean[names(tcga.tpm.mean) == "FOLR1"])] # 260
upper.sample <- colnames(tcga.tpm.fil)[which(tcga.tpm.fil[rownames(tcga.tpm.fil) == "FOLR1",] >= tcga.tpm.mean[names(tcga.tpm.mean) == "FOLR1"])] # 170

## extract tpm data of folr1 and muc16
lower.tcga.tpm <- tcga.tpm.fil[rownames(tcga.tpm.fil) %in% c("FOLR1", "MUC16"),colnames(tcga.tpm.fil) %in% lower.sample]
upper.tcga.tpm <- tcga.tpm.fil[rownames(tcga.tpm.fil) %in% c("FOLR1", "MUC16"),colnames(tcga.tpm.fil) %in% upper.sample]
lower.gtex.tpm <- gtex.tpm.fil[rownames(gtex.tpm.fil) %in% c("FOLR1", "MUC16"),]

## data used for boxplot
ov.lower.df <- data.frame(t(lower.tcga.tpm))
ov.lower.df$group <- "Low_FOLR1"
ov.upper.df <- data.frame(t(upper.tcga.tpm))
ov.upper.df$group <- "High_FOLR1"
nor.df <- data.frame(t(lower.gtex.tpm))
nor.df$group <- "Normal"
box.df <- rbind(ov.lower.df, ov.upper.df, nor.df)
box.df$FOLR1 <- log2(box.df$FOLR1 + 1)
box.df$MUC16 <- log2(box.df$MUC16 + 1)
box.df$group <- factor(box.df$group, levels = c("Normal", "Low_FOLR1", "High_FOLR1"))



my_comparisons <- list( c("Normal", "Low_FOLR1"), c("Normal", "High_FOLR1"), c("Low_FOLR1", "High_FOLR1") )

ggboxplot(box.df, x = "group", y = "MUC16",
          color = "group", palette = "lancet",
          add = "jitter") +
  stat_compare_means(comparisons = my_comparisons)
```

<img src="{{< blogdown/postref >}}index.en_files/figure-html/unnamed-chunk-7-1.png" width="672" />

```r
ggboxplot(box.df, x = "group", y = "FOLR1",
          color = "group", palette = "lancet",
          add = "jitter") +
  stat_compare_means(comparisons = my_comparisons)
```

<img src="{{< blogdown/postref >}}index.en_files/figure-html/unnamed-chunk-7-2.png" width="672" />
