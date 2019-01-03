---
title: 16s data processing and microbiome analysis
keywords: 16s rRNA, microbiome
tags: [genomics_analysis]
last_updated: September 17, 2014
summary: "This document briefly outlines the processing of 16s rRNA data, and provides guidance to key concepts and terms used"
sidebar: 16Smicro_sidebar
permalink: 16s-rRNA_6.html
folder: genomics_analysis/16s-rRNA/
author_profile: true
authors:
 - Gerrit_Botha
---



## H3ABioNet Assessment exercises
 
### Practice dataset {#practice_dataset}
The input datasets and metadata can be accessed [here](http://h3data.cbio.uct.ac.za/assessments/16SrRNADiversityAnalysis/practice).

### Input data assessment Questions {#input-assessment-questions} 

* Were the number, length and quality of the reads obtained in line with what would be expected for the sequencing platform used?
* Was the input dataset of sufficiently good quality to perform the analysis?
* How did the reads’ quality and GC content affect the way analysis was run?

### Operational assessment questions {#operational-assessment-questions}
* At each step of the workflow, describe which software was used and why:
  * Was the choice affected by the nature and/or quality of the reads?
  * Was the choice made due to the time and cost of the analysis?
  * What are the accuracy and performance considerations for the chosen piece of software?
* For each software, describe which input parameters were chosen, and why:
  * Was the choice affected by the nature and/or quality of the reads?
  * Did the available hardware play a role in the parameter choice?
  * How did the purpose of the study affect the parameter choice?
* For each step of the workflow, how do you know that it completed successfully and that the results are usable for the next step?

### Runtime analysis {#runtime-analysis}

This is useful information for making predictions for the clients and collaborators

* How much time and disk space did each step of the workflow take?
* How did the underlying hardware perform? Was it possible to do other things, or run other analyses on the same computer at the same time?

### Analysis of the results {#results-analysis}
* What percentage of the reads were removed during the quality trimming step?  Did all samples have similar number of reads after the preprocessing of reads steps? What was the median, maximum and minimum read count per sample? How many reads were discarded due to ambiguous bases?
* What percentage of reads could not be stitched? Were unstitched reads retained or discarded?
* How many chimeras were detected?
* How does the trimming or filtering strategy affect the the number of OTUs picked and the  classification and phylogenetic analysis of the OTUs?
* How does the % similarity threshold used during OTU picking affect the number of OTUs identified and the classification and phylogenetic analysis of the OTUs?
* How many OTUs were picked? What percentage of  the OTUs could be classified to the genus and species level? What percentage of OTUs could only be assigned to taxonomic ranks higher than genus?  What is the confidence threshold for the classifications?
* Does the use of a different 16S rRNA database for classification affect the results (e.g. were a lesser or greater number of OTUs classified to lower taxonomic ranks (genus, species))? Were any OTUs classified differently?
* Did the samples have enough sequence depth to capture the diversity? Did the rarefaction curve flatten? Should any samples be excluded because of low read count?
* Were there any differences in the alpha diversity between the samples in the different metadata categories (e.g. higher phylogenetic diversity in treatment 1 vs. treatment 2)?
* When groups of samples were compared (e.g. treatment 1 vs. treatment 2) based on distance metrics, such as unifrac, was there any particular clustering pattern observed?	 	 	
* Were 	any of the OTUs significantly correlated to any of the treatments or other metadata?





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
