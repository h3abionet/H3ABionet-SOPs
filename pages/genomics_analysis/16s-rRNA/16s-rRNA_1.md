---
title: 16s data processing and microbiome analysis
keywords: 16s rRNA, microbiome
tags: [genomics_analysis]
last_updated: September 17, 2014
summary: "This document briefly outlines the processing of 16s rRNA data, and provides guidance to key concepts and terms used"
sidebar: 16Smicro_sidebar
permalink: 16s-rRNA_1.html
folder: genomics_analysis/16s-rRNA/
author_profile: true
authors:
 - Gerrit_Botha
---



## Introduction {#introduction}

The genes encoding the RNA component of the small subunit of ribosomes, commonly known as the 16S rRNA in bacteria and archaea, are among the most conserved across all kingdoms of life. Nevertheless, they contain regions that are less evolutionarily constrained and whose sequences are indicative of their phylogeny. Amplification of these genomic regions by PCR from an environmental sample and subsequent sequencing of a sufficiently large number of individual amplicons enables the analysis of the diversity of clades in the sample and a rough estimate of their relative abundance. The analytical process is known as "16S rDNA diversity analysis", and is the focus of the present SOP.

The SOP describes the essential steps for processing 16S rRNA gene sequences. The procedure and tools are only recommendations and it is up to the user to evaluate what works best for their needs.

### Glossary of terms and jargon {#glossary}

<table>
  <tr>
   <td>16S rRNA gene
   </td>
   <td>The gene that is responsible for the coding of the 16S ribosomal RNA. The gene is used in constructing phylogenies.
   </td>
  </tr>
  <tr>
   <td>Barcodes
   </td>
   <td>Short nucleotide sequences added onto the ends of the DNA fragments that are to be sequenced. It allows for indexing of samples, so multiple DNA libraries can be mixed together into one sequencing lane.
   </td>
  </tr>
  <tr>
   <td>Variable region
   </td>
   <td>16S rRNA gene sequences contain hypervariable regions that can provide clade-specific signature sequences useful for bacterial identification.
   </td>
  </tr>
  <tr>
   <td>Demultiplex
   </td>
   <td>This is a process of binning reads based on barcodes, primarily used to split them amongst samples.
   </td>
  </tr>
  <tr>
   <td>Operational Taxonomic Unit (OTU)
   </td>
   <td>An operational taxonomic unit is an operational definition of a species or group of species often used when only DNA sequence data is available.
   </td>
  </tr>
  <tr>
   <td>Rarefaction analysis
   </td>
   <td>Rarefaction is a process used to estimate the true diversity of a sample by extracting random subsets of sequences. The analysis estimates diversity from subsets of different sizes and extrapolates the resulting rarefaction curve to an infinite number of sequences. Rarefaction is also used to determine whether the sequencing depth achieved (number of reads per sample) is sufficient to capture the diversity within a sample. A plot is generated showing increase in number of species (or other diversity metrics) as the number of sequence reads increase. A curve that is reaching asymptote indicates that no further diversity would be expected if sequencing depth was increased.
   </td>
  </tr>
  <tr>
   <td>Phred quality score or Q score
   </td>
   <td>Measure of sequencing accuracy. Logarithmically related to the probability that a base is called incorrectly by a sequencer. Example: a Phred score of 30 (Q30) means that the probability of an incorrect call for that base is 1 in 1000 and the base call accuracy is 99.9%. The base call accuracy for a base with Q10 is 90%, Q20 is 99%, Q30 is 99.9%, Q40 is 99.99%, and Q50 = 99.999%.
   </td>
  </tr>
  <tr>
   <td>Adapter
   </td>
   <td>Platform specific nucleotide sequence added to the ends of DNA molecule to facilitate sequencing e.g. in Illumina the adapter facilitates binding to the complementary target sequences mobilised on the flow cell.
   </td>
  </tr>
  <tr>
   <td>Chimera
   </td>
   <td>PCR artefact. Chimeras are potentially formed during PCR when incompletely extended DNA fragments from different templates anneal due to closely related sequences generating recombinants between starting templates. Chimeras can greatly impact estimates of diversity (generally overestimate).
   </td>
  </tr>
  <tr>
   <td>Alpha diversity
   </td>
   <td>Diversity within a single sample. Diversity can be characterised using the number of different species (richness), the abundance of each species (evenness), with indices that combine richness and evenness, and with divergence-based methods (phylogenetic diversity).
   </td>
  </tr>
  <tr>
   <td>Beta diversity
   </td>
   <td>Comparison of diversity between samples.
   </td>
  </tr>
  <tr>
   <td>Unifrac
   </td>
   <td>Beta diversity distance metric based on the phylogenetic distance between the members of communities/samples. Unifrac captures the total amount of evolution that is unique to each sample.
   </td>
  </tr>
</table>

### Schematic workflow of the analysis {#workflow}

| ![16s analysis pipeline]({{ site.baseurl }}/assets/images/16s-overview.png "image_tooltip") |
| :--: |
| Figure 1. Steps in 16s Analysis Workflow |





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
