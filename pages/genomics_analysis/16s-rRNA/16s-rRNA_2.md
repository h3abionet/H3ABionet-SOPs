---
title: 16s data processing and microbiome analysis
keywords: 16s rRNA, microbiome
tags: [genomics_analysis]
last_updated: September 17, 2014
summary: "This document briefly outlines the processing of 16s rRNA data, and provides guidance to key concepts and terms used"
sidebar: 16Smicro_sidebar
permalink: 16s-rRNA_2.html
folder: genomics_analysis/16s-rRNA/
author_profile: true
authors:
 - Gerrit_Botha
---



## Phase 1: Preprocessing of reads {#preprocessing}

It is essential to preprocess raw reads before subjecting them for downstream analysis. The preprocessing includes the removal of low quality bases, ambiguous bases and adapter sequences, the stitching together of paired reads, and the detection of chimeric reads. Sequencing errors, reads with ambiguous bases and chimeras can all cause the appearance of spurious OTUs if they are not removed.

_Input_: raw reads (multiplexed or demultiplexed)

_Output_: high quality reads ready for OTU picking

### **_QC plots and stats_** {#qc_plots}

The first step in the data preprocessing is to check the quality of bases in all the reads. Once we understand the quality spectrum of the reads, we can decide on the parameters for trimming low quality bases. If the raw reads are not demultiplexed, demultiplexing should be performed before proceeding to the next step. Illumina's CASAVA or QIIME tools can perform the demultiplexing task.

<span style="text-decoration:underline;">Software</span>: FASTQC, PRINSEQ, SolexaQA

### **_Trim and Filter reads_** {#trim}

At the 3' end of reads there are often adaptor sequences left from library preparation. These adaptor bases need to be removed, and low quality bases need to be trimmed off. Any of the indicated programs can be used for this. Bokulich et al. 2013 [^1] recommend a minimum phred quality score of 3 to trim low quality bases at the ends of the reads. Jeraldo et al. 2014 (in review) recommend trimming the 3' end of the reads with a moving average score of 15, with a window size of 4 bases and to removal of any reads shorter than 75% of the original read length. It is also recommended that reads containing ambiguous bases (N) be discarded.

<span style="text-decoration:underline;">Software</span>: Trimmomatic, PRINSEQ, SolexaQA

### **_Paired read stitching_** {#stitching}

When the combined length of reads sequenced from both ends of DNA fragments is longer than the size of the fragment, there is an overlap between the paired reads. The read pairs can be stitched together based on the overlap information, thus generating a single sequence. During the read stitching process, higher quality bases can be selected thus improving the quality of stitched reads. PANDASeq does not do well when the overlap is almost the entire read. PEAR works well for all lengths of overlaps between the paired reads.

<span style="text-decoration:underline;">Software</span>: PEAR, PANDASeq, FLASH, UPARSE merge

### **_Chimera detection_** {#chimera}

Chimeras are artifacts of PCR. These are formed during PCR cycles by the joining of two or more different parent DNA templates. If these chimeras are not removed they may be recognized as novel sequences during the alignment process and therefore mislead interpretation. At present there are no tools that can remove chimeras completely without throwing away non-chimeric sequences. Of the various tools available, UCHIME was found to perform better than ChimeraSlayer, which was the best program to detect chimeras before UCHIME was developed[^2]. To remove chimeras from 454 sequences Perseus can also be used[^3].

<span style="text-decoration:underline;">Software</span>: UCHIME, ChimeraSlayer, Perseus





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
