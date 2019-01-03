---
title: 16s data processing and microbiome analysis
keywords: 16s rRNA, microbiome
tags: [genomics_analysis]
last_updated: September 17, 2014
summary: "This document briefly outlines the processing of 16s rRNA data, and provides guidance to key concepts and terms used"
sidebar: 16Smicro_sidebar
permalink: 16s-rRNA_3.html
folder: genomics_analysis/16s-rRNA/
author_profile: true
authors:
 - Gerrit_Botha
---



## Phase 2: OTU picking, classification and phylogenetic tree generation {#phase2}

During this phase reads are processed so that comparisons between samples can be made. The first step is to cluster reads based on similarity into OTUs and to select a representative sequence for each OTU. Each OTU is then classified by comparison to a reference database and a phylogeny inference is made based on sequence alignment and the construction of a phylogenetic tree.

_Input_: high quality reads

_Output_: OTUs, representative sequences, OTU table with classification and abundance of each OTU, heatmap, sequence alignment and phylogenetic tree

### **_OTU picking_** {#otu_picking}

OTU picking is the clustering of the preprocessed reads into OTUs. The clusters are formed based on sequence identity. The identity threshold can be defined by the user. Sequences that are more than 97% identical are conventionally assumed to be derived from the same bacterial species/OTU. Other identity percentages can be used, depending on the granularity of the desired clusters and the known divergence in 16S sequences of the OTUs of interest. Three approaches for OTU picking exist. 1) de novo OTU picking groups sequences based on levels of pairwise sequence identity; 2) closed reference OTU picking aligns and groups sequences relative to a reference database, and sequences that are not >97% identical to a known reference are discarded 3) open-reference OTU picking starts with alignment to a reference database, but if the read does not match a known sequence it is not discarded but sent for de novo OTU picking. After the sequences have been clustered into OTUs and counted to estimate OTU abundance, a representative sequence is picked for each OTU. Each OTU is therefore represented by a single sequence and this will speed up downstream analysis. There are multiple choices to select a representative sequence. It can be the first sequence, the longest sequence, the seed sequence used in OTU picking, the most abundant sequence or a random sequence.

<span style="text-decoration:underline;">Software</span>: UPARSE, QIIME

### **_Classification_** {#classification}

Here a taxonomic identity  is assigned to each representative sequence. The taxonomies are pulled from a reference set. There are three main reference databases with aligned, validated and annotated 16S rRNA genes: GreenGenes, Ribosomal Database Project (RDP) and Silva. Each of these databases has strengths and weaknesses that need to be taken into consideration, and all are in common use. Several methods for assigning taxonomy against these reference databases exist, including UCLUST, the RDP classifier and RTAX.

<span style="text-decoration:underline;">Databases</span>: SILVA, GreenGenes, RDP

<span style="text-decoration:underline;">Software</span>:  UCLUST, the RDP classifier, RTAX

### **_Alignment_** {#alignment}

To understand the evolutionary relationships between the sequences in the sample and to perform a diversity analysis, it is necessary to generate a phylogenetic tree of the OTUs. The first step in generating the tree is to generate a multiple alignment of the representative OTU sequences. PyNAST aligns the sequences to a template alignment of reference 16S sequences. Infernal makes use of a Hidden Markov Model that also incorporates secondary structure information.

<span style="text-decoration:underline;">Software</span>:  PyNAST, INFERNAL

### **_Create phylogenetic tree_** {#phylogene}

The phylogenetic tree represents the relationship between the sequences in terms of the evolutionary distance from a common ancestor. In downstream analysis this tree is used for example in calculating the UniFrac distances.

<span style="text-decoration:underline;">Software</span>:  FastTree

An alternate option to most of the steps mentioned in phase 1 and phase 2 is to run IM-TORNADO (Jeraldo et al. 2014 in review). IM-TORNADO is an integrated pipeline that takes demultiplexed reads and trims low quality bases, does paired read stitching, removes chimeras and generates OTU table, phylogenetic tree and assigns taxonomy. Unique feature of IM-TORNADO is that it can analyze paired reads that do not overlap. Non-overlapping paired reads are typically generated from non-overlapping variable regions of 16S rRNA. Such studies try to utilize information in two variable regions instead of one variable region as in any standard 16S rRNA study to define OTUs.





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
