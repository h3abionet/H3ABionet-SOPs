---
title: 16s data processing and microbiome analysis
keywords: 16s rRNA, microbiome
tags: [genomics_analysis]
last_updated: September 17, 2014
summary: "This document briefly outlines the processing of 16s rRNA data, and provides guidance to key concepts and terms used"
sidebar: 16Smicro_sidebar
permalink: 16s-rRNA_4.html
folder: genomics_analysis/16s-rRNA/
author_profile: true
authors:
 - Gerrit_Botha
---



## **Phase 3: Measure diversity and other statistical analysis** {#phase3}



OTU information (number of OTUs, abundance of OTUs) and the phylogenetic tree generated from the phase 2 is utilized to estimate diversity within and between samples. Additional statistical analysis to test the significance of the diversities can also be done.

_Input_: classified OTU table with abundance, phylogenetic tree and sample metadata

_Output_: alpha and beta diversity metrics, distance matrix, results from statistical tests, rarefaction plots, PCoA plots, heatmaps

### **_Determine alpha diversity_** {#alpha_diversity}

Alpha diversity is a measure of diversity within a sample. It gives an indication of richness and/or evenness of species present in a sample. The accuracy of the measured diversity is mainly affected by the sequencing depth (number of reads per sample). Sequencing depth must be high enough to capture the true diversity within a sample. Samples with higher number of reads would show higher diversity than samples with lower number of reads. Rarefaction analysis is therefore required to understand the actual diversity within a sample and to determine if your sequencing effort is sufficient and if the total diversity within the sample has been captured. Mothur as well as QIIME have tools to generate multiple rarefactions and then measure alpha diversity on the rarefied OTU tables. Several popular alpha diversity measures are available both in Mothur and QIIME: Shannon index, chao1, observed species, and phylogenetic diversity whole tree.

<span style="text-decoration:underline;">Software</span>:  mothur, QIIME

### **_Determine beta diversity_** {#beta_diversity}

Beta diversity is a measure of diversity between samples. One of the most commonly used metrics is the Unifrac distance that compares samples using phylogenetic information. An all-by all or pairwise matrix of the beta diversity metrics between all the samples in the study is generated and can be visualized in different ways such as a tree, graph, network etc. Mothur and QIIME have several tools to generate distance metrics, phylogenetic trees and PCoA plots.

<span style="text-decoration:underline;">Software</span>:  UNIFRAC for distance metrics, mothur, QIIME

### **_Other statistical analysis_** {#other}

Additional statistical tests between samples or groups of samples can be done in QIIME. For alpha diversity a parametric or non-parametric t-test can be performed on a rarefied number of sequences. For beta diversity the Mantel, partial Mantel and Mantel correlogram matrix correlation can be used to compare distance matrices. Multivariate analyses are also available for testing significance between the distance matrix and other factors. The statistical methods available are: adonis, ANOSIM, BEST, Moran's I, MRPP, PERMANOVA, PERMDISP, and db-RDA. Native methods in R and other R packages such as phyloseq and ade4 can also be considered for these types of analyses.

<span style="text-decoration:underline;">Software</span>:  QIIME, R packages (phyloseq, ade4)





## **References** {#references}

[^1]: Bokulich, Nicholas A., et al. ["Quality-filtering vastly improves diversity estimates from Illumina amplicon sequencing."](https://www.nature.com/articles/nmeth.2276) Nature methods 10.1 (2013): 57.

[^2]:  Edgar, Robert C., et al. ["UCHIME improves sensitivity and speed of chimera detection."](https://academic.oup.com/bioinformatics/article/27/16/2194/255262) Bioinformatics 27.16 (2011): 2194-2200.

[^3]: Quince, Christopher, et al. ["Accurate determination of microbial diversity from 454 pyrosequencing data."](https://www.nature.com/articles/nmeth.1361) Nature methods 6.9 (2009): 639.


[//]: <> (Below are the common abbreviations in the page.)
*[SOPs]: Standard Operating Procedures
*[16S rRNA]: 16S ribosomal RNA
*[OTU]: Operational Taxonomic Unit
*[Q score]: Phred Quality score, a measure of sequencing accuracy.
*[PCR]:
