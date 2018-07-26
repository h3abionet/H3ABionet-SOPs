---
layout: single 
title: RNA-Seq data processing and gene expression analysis
author: "Jenny Zadeh, Radhika Khetani, Jessica Holmes, Chris Fields"
date: "July 4, 2018"
version: "3.0"
#sidebar:
#  nav: "rnatoc"
toc: true
#header:
#  image: /assets/images/H3ABioNet_high_res.png
author_profile: true
authors: 
 - Jenny_Zadeh
 - Radhika_Khetani
 - Jessica_Holmes
 - Chris_Fields
---

## Introduction {#introduction}

This document briefly outlines the essential steps in the process of analyzing gene expression data using RNA sequencing (mRNA, specifically). and recommends commonly used tools and techniques for this purpose. It is assumed in this document that the experimental design is simple and that differential expression is being assessed between 2 experimental conditions, i.e. a simple 1:1 comparison, with some information about analyzing data from complex experimental designs. The focus of the SOP is on single-end strand-specific reads, however special measures to be taken for analysis of paired-end data are also briefly discussed. The recommended coverage for RNA-Seq on human samples is 30-50 million reads (single-end), with a minimum of three replicates per condition, preferably more if one can budget accordingly.  Preference is also generally given for a higher number of replicates with a lower per-sample sequence yield (15-20 million reads) if there is a tradeoff between the number of reads per sample and the total number of replicates.

The procedures outlined below are recommendations to the H3ABioNet groups planning to do differential gene expression analysis on human RNA-Seq data, and are not meant to be prescriptive. Our goal is to help the groups set up their procedures and workflows, and to provide an overview of the main steps involved and the tools that can be used to implement them.


### Glossary of associated terms and jargon  {#glossary-of-associated-terms-and-jargon}

*   FASTQ format & quality scores
    *   FASTQ format is the standard format of raw sequence data
    *   Quality scores assigned in the FASTQ files represent the probability that a certain base was called incorrectly. These [scores are encoded ](http://en.wikipedia.org/wiki/FASTQ_format#Encoding)in various ways and it is important to know the type of encoding for a given FASTQ file.
*   Single-end vs paired-end, read vs fragment
    *   DNA fragments can be sequenced from one end or both ends, single-end (SE) or paired-end (PE) respectively. During data processing PE data will often be represented as **fragments** consisting of 2 reads, instead of 2 separate reads.
*   Strandedness
    *   During library prep RNA can be processed into cDNA such the strand information is maintained. This is important in regions where there are overlapping genes on the two DNA strands
*   TPM, RPKM, FPKM, CPM
    *   TPM - Transcripts Per Million
    *   RPKM - Reads Per Kilobase per Million mapped reads (for SE data)
    *   FPKM - Fragments Per Kilobase per Million mapped reads (for PE data)
    *   CPM - Counts Per Million
*   TMM
    *   This is a type of normalization and is an acronym for "Trimmed Mean of Ms" ([Robinson and Oshlack, Genome Biology, 2010](http://genomebiology.com/2010/11/3/R25/abstract)).




## Procedural steps {#procedural-steps}

[This protocol paper ](http://www.nature.com/nprot/journal/v7/n3/full/nprot.2012.016.html)(Trapnell _et al. _2012) was a very good resource for understanding the procedural steps involved in any RNA-Seq analysis. The datasets they use in that paper are freely available, but the source of RNA was the fruitfly _Drosophila melanogaster_, and not Human tissue. In addition, they exclusively use the "tuxedo" suite developed in their group.

Several papers are now available that describe the steps in greater detail for preparing and analyzing RNA-Seq data, including using more recent statistical tools:



*   (Love _et al. _2016) [https://f1000research.com/articles/4-1070/v2](https://f1000research.com/articles/4-1070/v2)
*   (Law _et al._ 2016) [https://f1000research.com/articles/5-1408/v2](https://f1000research.com/articles/5-1408/v2)

In addition, newer alignment-free methods have also been published and are increasingly being used in analysis (we include a second protocol detailing the use of these):



*   (Bray _et al. _2016) [https://www.nature.com/articles/nbt.3519](https://www.nature.com/articles/nbt.3519)
*   (Patro _et al._ 2017) [https://www.nature.com/articles/nmeth.4197](https://www.nature.com/articles/nmeth.4197)

Tools are suggested in the protocols below.



---
**Figure 1.** Steps in the Workflow.   


![RNA-Seq analysis pipeline](/assets/images/RNA-Seq-Overview.png "image_tooltip")



### _<span style="text-decoration:underline;">Phase 1: Preprocessing of the raw reads</span>_ {#phase-1-preprocessing-of-the-raw-reads}

The following steps prepare reads for analysis and should be always performed prior to alignment.


#### _Step 1.1: Quality check _ {#step-1-1-quality-check}

The overall quality of the sequence information received from the sequencing center will determine how the quality trimming should be set up in Step 1.2. Tools like [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) will enable the collection of this information. Sequencing facilities usually produce read files in FASTQ format, which contain a base sequence and a quality score for each base in a read. FastQC measures several metrics associated with the raw sequence data in the FASTQ file, including read length, average quality score at each sequenced base, GC content, presence of any overrepresented sequences (k-mers), and so on. The key metric to watch for is the graph representing the average quality scores (see figure 2), and the range of scores at each base along the length of the reads (reads are usually the same length at this time, and this length is the X-axis, the Y-axis has the quality scores). Note that for large projects, you may collate all of the FastQC reports by using a tool like [MultiQC](http://multiqc.info). MultiQC will generate an html file that visually summarizes these metrics across all samples, as well as provide tab-delimited files containing all the FastQC stats.

_Note: FastQC has very stringent criteria to assess whether the data "Pass" or "Fail" for a given metric it measures, so even if it looks like your data has "failed" with respect to a given metric, please read carefully about the criteria employed. In most situations a "failed" reading for multiple metrics is not a death sentence for the dataset._



---
**Figure 2.** Graphs generated by FastQC detailing the average quality scores across all reads at each base.  

![FastQC per-base quality scores](/assets/images/RNA-Seq-FASTQC.png "image_tooltip")

---

#### _Step 1.2: Adaptor and Quality trimming + Removal of very short reads_ {#step-1-2-adaptor-and-quality-trimming-removal-of-very-short-reads}

In this step we deal with** 3 major preprocessing steps** that clean up the data and reduce noise in the overall analysis.

1.  Adaptors (glossary term) are artificial pieces of DNA introduced prior to sequencing to ensure that the DNA fragment being sequenced attaches to the sequencing flow cell. Usually these adaptors get sequenced, and have already been removed from the reads. But sometimes bits of adaptors are left behind, anywhere from 90% to 20% of the adaptor length. These need to be removed from the reads. The adaptor sequence for this step will have to be obtained from the same source as the sequence data.

2.  Frequently, the quality of bases sequenced tends to drop off toward one end of the read. A low quality base call means that the nucleotide assigned has a higher probability of being incorrect (see [this link](https://en.wikipedia.org/wiki/Phred_quality_score) for a more in-depth overview of quality scores). It is best to trim off any low quality bases at the ends of reads to ensure the best alignment to the reference. Usually a quality score of <25 is considered as a "poor" quality score.

3.  Once the adaptor remnants and low quality ends have been trimmed, some reads may end up being very short (i.e. <20 bases). These short reads are likely to align to multiple (wrong) locations on the reference, introducing noise. Hence any reads that are shorter than a predetermined cutoff (e.g. 20) need to be removed

One tool that deals with all of these issues at once is [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), though there are there are various alternatives that can perform these 3 clean up steps either combined or one after the other; these are  listed below. For data that are paired ended, it is very important to perform the trimming for both read1 and read2 simultaneously. This is because all downstream applications expect paired information, and if one of the 2 reads is lost because it is too short, then the other read becomes unpaired (orphaned) and cannot be used properly for most applications. Trimmomatic has 2 modes, one for single end data (SE) and another one for paired end data (PE). If using paired end reads, please be sure to use the PE mode with both read1 and read2 FASTQ  files for the same run.

**Alternative Tools:**


    _For adaptor trimming_ - Trim_Galore, BBMap, Flexbar (Dodt _et al._ 2012) and one of the many tools [listed here](https://omictools.com/adapter-trimming-category).


    _For trimming low quality bases from the ends of reads** **_- Trim_Galore, BBMap, FASTX-Toolkit (fastq_quality_filter), PrinSeq, SolexaQA (Cox _et al._ 2010).


    _For removing very short reads_ - PrinSeq, Trim_Galore

**_Step 1.3: Quality recheck _**

Once the trimming step is complete, it is always good practice to make sure that your dataset looks better by rerunning FastQC on the trimmed data. The metrics to compare between trimmed and raw fastq data, in the context of the tool FastQC are listed below:



*   Sequence length and Sequence length distribution (the minimum will be lower in trimmed data)
*   The quality score graph (majority of the bases will have a minimum quality value at or above 30, see Figure 3 below)
*   Per base sequence content (should remain about the same)





---
**Figure 3.** Graphs generated by FastQC detailing the change in average base quality across all reads after trimming in an example dataset.   

![FastQC average quality score before and after trimming](/assets/images/RNA-Seq-FASTQC2.png "image_tooltip")

---


### _<span style="text-decoration:underline;">Phase 2: Determining how many read counts are associated with known genes</span>_ {#phase-2-determining-how-many-read-counts-are-associated-with-known-genes}


#### _Step 2.1: Generation of gene/transcript-level counts_ {#step-2-1-generation-of-gene-transcript-level-counts}

Below are two protocols for generation of gene- and transcript-level counts.  Protocol 1 focuses on the classi alignment-based approach, whereas Protocol 2 uses recent advances in accelerated kmer-based 'pseudo-alignment' approaches for assigning reads to transcripts, which has some distinct advantages to alignment methods but rely on having a comprehensive and reliably annotated transcriptome.


##### **Protocol 1: Alignment-based approach**

Once the data are cleaned up, the next step is alignment to the reference genome. There are various tools available for this step, but it is important that the alignment tool chosen here is a "splice-aware" tool. That means that the tool should have the capability to align reads that contain exonic sequences from 2 exons on either side of one intron (also called intron-spanning reads). STAR, HISAT2, GSNAP, SOAPSplice are some of the many splice-aware aligners available. Note that TopHat used to be a very commonly-used tool from the Tuxedo suite, however this software is no longer supported and has been superseded by HISAT2 (recommended for human data).

When performing alignments it is imperative to set up the parameters properly to ensure the best alignment. Irrespective of the aligner used, it needs the following information:



*   Are the data made up of single-end reads or paired-end reads?
*   Are the data stranded, if so, was the standard dUTP method employed (STAR can detect this automatically)?

The other information that should be provided when setting up the alignment is the gene annotation information, a **GTF **or** GFF3** file that contains the information about the location of all the genes in the context of the reference genome (a **FASTA** file). It is very important to pick the gene annotation file that corresponds to the reference genome, i.e. the same version number and from the same source (Ensembl, UCSC or NCBI).

Note, that most all aligners require an index file to be created from the reference genome (a **FASTA** file). This helps speed up alignment drastically. For STAR, this index can be created using the "genomeGenerate" mode with or without an annotation file.

Once the alignment is complete, the final result will be a file in **SAM** or **BAM** format. For a STAR run, the main alignment output will end in "Aligned.out.sam" by default, but it may also be returned in a sorted or unsorted **BAM **file. In addition, there are two log files that are returned from STAR that report the progress of the run and save a summary of the final results. There is also a splice junctions (SJ) file that details high confidence splice junctions.



Once the alignment is completed, the first step is to check how many reads aligned to the genome. For STAR, all of these details can be found in the "Log.final.out" file. For RNA-Seq on Human samples, for good quality data, about 70 - 90% of the reads should match somewhere on the genome. If the data in question are of good quality but < 60 % of the reads are mapping to the genome, it is worth evaluating the parameters, and testing unmapped reads for presence of potential contaminants.


##### **Protocol 2: Estimated counts using graph-based or similar approaches**

A more recent alternative to alignment-based methods both dramatically increases the speed of the analysis and resolves some of issues that make analysis of transcripts or gene families problematic using alignment-based approaches.  These methods use an alternative approach that performs essentially a very lightweight 'alignment' that speeds up analysis, sometimes by orders of magnitude.  These approaches also generate (as output) estimated counts that can be imported into standard R-based workflows, thus combining the initial two steps in the alignment-based approach above.

The speedups are based on the tool being used and are accomplished in slightly different ways, such as mapping kmers from the reads to a transcriptome-based de bruijn graph (exemplified by [kallisto](https://pachterlab.github.io/kallisto/)) or 'quasi-mapping' of reads to transcript positions in a simple transcriptome-based index (e.g. [Salmon](https://combine-lab.github.io/salmon/)).  These are normally followed by an expectation maximization (EM) step to resolve ambiguous assignments, re-proportioning reads based on evidence from the overall analysis.  Estimates of read counts to the transcripts can then be generated and used in downstream analyses.  For more background, a good independent summary of the 'pseudo-alignment' approach is found [here](http://tinyheero.github.io/2015/09/02/pseudoalignments-kallisto.html).  The EM step in particular has proven useful in finding additional genes or transcripts in data that were missed using alignment-based approaches, which normally skip ambiguously mapped sequences.

A key difference in these procedures from the alignment approach is the tools require a transcriptome data set (not a reference genome).  This may be a problem if your reference genome annotation isn't of reasonably high quality, for instance if the transcripts described aren't well-annotated or incomplete.  However, these tools are of great use for well-characterized genomes such as human and mouse, and can also be used with transcriptome assemblies.  

Note, as these analyses generate estimated read counts as part of their output, you can skip Step 2.2


#### _Step 2.2: Count generation_ {#step-2-2-count-generation}


##### **Protocol 1**

Once it is determined that the alignment step was successful, the next step is to enumerate the number of reads that are associated with the genes. There are multiple tools to perform this step (e.g. HTSeq's htseq-count, Subread's featureCounts), in addition, there are statistical analysis tools that do not require this step (e.g. Cufflinks' Cuffdiff and Ballgown). Both scenarios will be tackled in Phase III.

We recommend using featureCounts to collect the raw gene count information; this tool will require the alignment file (**BAM**), and the associated gene annotation file (**GTF**). Irrespective of the counting tool used, it needs the following information:



*   Different sources have slightly different formats, so it is essential to specify how the counting needs to be performed, irrespective of the counting tool employed. Gene counts should be collected for each gene ("-g" set as gene_id, for featureCounts), and at the level of the exon ("-t" set as exon, for featureCounts).
*   Are the data stranded, if so, was the standard dUTP method employed ("-s" set as 2 for reverse, using featureCounts)?

The final output of any of these programs is a tab-delimited file with gene names in column one and counts in the second column. featureCounts will return an additional file that ends in ".summary" that specifies the number of reads that did not map only to one gene, split into various categories. It is normal for the total sum of all the rows in this file to be higher than the number of aligned reads for a sample, because if one read maps to two locations, featureCounts will count it twice in the "Unassigned_MultiMapping" category.


##### **Protocol 2**

Counts are already generated and can be skipped


#### _Step 2.3: Collecting and tabulating alignment stats_ {#step-2-3-collecting-and-tabulating-alignment-stats}


##### **Protocol 1**

For a given RNA-Seq run it is valuable to collect several stats related to the alignment and counting steps; this is an important step in the evaluation process. There are many tools that gather this information from the **SAM** or **BAM **output file, e.g samtools' flagstat, Picard's CollectAlignmentSummaryMetrics. However there are quirks with each reporting tool, hence it is recommended to collect this information as described in the table below.


<table>
  <tr>
   <td>Gather numbers for the following categories
   </td>
   <td>Calculating or gathering the information
   </td>
   <td>Tool generating the information
   </td>
  </tr>
  <tr>
   <td><p style="text-align: right">
<strong>Total reads</strong></p>

   </td>
   <td>"Total Sequences" in the "Basic Statistics" section
   </td>
   <td>FastQC
   </td>
  </tr>
  <tr>
   <td><p style="text-align: right">
<strong>Total reads after trimming</strong></p>

   </td>
   <td>"Total Sequences" in the "Basic Statistics" section (FastQC) <strong>OR</strong> "Number of input reads" in "Log.final.out" file (STAR)
   </td>
   <td>FastQC <strong>OR</strong> STAR (<em>prefix</em>_Log.final.out)
   </td>
  </tr>
  <tr>
   <td><p style="text-align: right">
<strong>Unmapped reads</strong></p>

   </td>
   <td>Take the "Number of input reads" and subtract "Uniquely mapped reads number", "Number of reads mapped to multiple loci", and "Number of reads mapped to too many loci"
   </td>
   <td>STAR (<em>prefix</em>_Log.final.out)
   </td>
  </tr>
  <tr>
   <td><p style="text-align: right">
<strong>Reads mapped to genome</strong></p>

   </td>
   <td>Add "Uniquely mapped reads number", "Number of reads mapped to multiple loci", and "Number of reads mapped to too many loci"
   </td>
   <td>STAR (<em>prefix</em>_Log.final.out)
   </td>
  </tr>
  <tr>
   <td><p style="text-align: right">
<strong>Multiply mapped reads</strong></p>

   </td>
   <td>Add "Number of reads mapped to multiple loci", and "Number of reads mapped to too many loci"
   </td>
   <td>STAR (<em>prefix</em>_Log.final.out)
   </td>
  </tr>
  <tr>
   <td><p style="text-align: right">
<strong>Reads mapped to genes</strong></p>

   </td>
   <td>This number is listed as  "Assigned" in ".summary" file
   </td>
   <td>featureCounts (<em>prefix</em>.txt.summary)
   </td>
  </tr>
  <tr>
   <td><p style="text-align: right">
<strong>Uniquely mapped reads without an associated gene</strong></p>

   </td>
   <td>This number is listed as "Unassigned_NoFeatures" in ".summary" file
   </td>
   <td>featureCounts (<em>prefix</em>.txt.summary)
   </td>
  </tr>
  <tr>
   <td><p style="text-align: right">
<strong>Uniquely mapped reads with an ambiguous gene assignment</strong></p>

   </td>
   <td>This number is listed as "Unassigned_Ambiguity" in ".summary" file
   </td>
   <td>featureCounts (<em>prefix</em>.txt.summary)
   </td>
  </tr>
</table>


You may also find it easier to collate all of the summary data from FastQC, STAR, featureCounts, and Picard's metrics by running MultiQC on all of the reports generated by these three programs. MultiQC will generate an html file that visually summarizes this data as well as tab-delimited files containing all the stats produced by these programs.


##### **Protocol 2**

Apart from FASTQC, other standard QC metrics that rely on an alignment are not available, such as Picard's tools, or a more complete assessment of read fates.  Salmon and kallisto both provide output that give basic overall statistics such and the number of reads mapped, and MultiQC can also summarize this information for all samples.




### _<span style="text-decoration:underline;">Phase 3: Statistical analysis</span>_ {#phase-3-statistical-analysis}


#### _Step 3.1: QC and outlier/batch detection_ {#step-3-1-qc-and-outlier-batch-detection}

Alignment-based counts can normally be imported directly into R/Bioconductor, but estimated counts generated using Protocol 2 above will require importing using the _[tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html)_ Bioconductor library.  This can also be used to summarize transcript-level information to gene counts.

As mentioned above, estimated transcript-level counts can be used in standard differential gene expression analyses steps, but these should be imported using tools like _[tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html)_ (for R/Bioconductor).   We also recommend using these tools to sum transcript-level information to gene-level counts for initial analyses.

Once the counts have been generated, it is good practice to do some QC checks in addition to the ones listed above and to cluster the samples to see if there are any outliers or batch effects. Batch effects are usually caused by obtaining or processing the samples in batches and can obscure detection of expression differences if not adjust for statistically.

Below are 2 examples of variation between samples plotted using the plotDensity function from the "affy" package in R (see note about R below). Whatever the shape of the distribution, ideally it will be about the same for all samples. One or two samples that are very different could be outliers, but if there are two or more distinct groups, see if they correspond to treatment groups or known or unknown batch effects.

There are other QC steps that are recommended, like simple hierarchical clustering (hclust + plot functions in base R) or Principal Component Analysis clustering (plotPCA function in the affy package). This will help validate the presence of outliers.

_Note that while StringTie has facilities for count generation, normalization and statistical analysis, it does not have any internal methods for sample QC or clustering. Instead, the R statistical software and add-on packages from Bioconductor are an excellent way to handle all aspects of statistical analysis of RNA-Seq data. They are free and available for any computer platform, although the command line interface can have a steep learning curve. Learning how to use R is well worth the time investment as it is a general tool for any sort of data manipulation, statistical analyses and graphing needs._


#### _Step 3.2: Remove low count genes and normalize_ {#step-3-2-remove-low-count-genes-and-normalize}

The QC investigations in Step 3.1 should be done on the output from htseq-count, which typically is the entire transcriptome (all genes). However, for most eukaryotic species, only ~40-60% of genes are expressed in any given cell, tissue or age and so a large proportion of the genes may contain 0 or only a few reads. It is recommended to remove these genes because they do not contain enough information to be statistically valid and it reduces the amount of multiple hypothesis testing (Step 3.3). A common practice is to require a minimum number of reads in a minimum number of samples. Often 1 CPM (count per million, relative to the total number of reads in the sample) is used as the minimum threshold but this may need to be adjusted up or down depending on the total number of reads per sample to get a threshold between 5 and 20 reads per gene.

The CPM adjustment is a minimal normalization that is necessary because of the normal large variation in total reads per sample. However, this assumes that the total number of reads _should_ have been the same for all samples. This assumption can often be violated by extremely high expression in a few genes in one treatment or if the DE genes are predominantly in one direction. The traditional FPKM normalization used in older tools such as Cufflinks' Cuffdiff cannot adjust for this and was one of the reasons FPKM was found to be inferior to other normalization methods, including Transcripts Per Million (or TPM; Dillies et al. 2012).

Which extra normalization, DESeq or TMM, to use in R depends on which package, DESeq2 or edgeR, you prefer to use in R for statistical analysis. Both use extra normalization methods that are comparable and adjust for moderate biases in the number and direction of gene expression changes. Both are based on negative binomial statistical modeling and were found to compare quite evenly. Both have functions for outputting normalized expression values for use in QC and downstream visualizations. These normalized expression values should be put through the same QC metrics in Step 3.1 to see if the extra normalization rectified any problems indicated. Severe outliers should be removed and the data re-normalized. Remaining batch effects should be adjusted for in the statistical model.


#### _Step 3.3: Statistics for differential expression_ {#step-3-3-statistics-for-differential-expression}

The main goal of most RNA-Seq is detection of differential expression between two or more groups. This is done for thousands of genes with often only a few replicates per group, so a statistical method must be used. This field is still under development and various statistical methods have been proposed. Both DESeq2 and edgeR can handle paired samples, other batch effects and complex designs. The choice between the two can simply be made based on which one is easier to understand and implement for the user.

In any statistical test of differential expression between groups, the amount of change from group to group must be evaluated by how much variation there is among the replicates of a group and how many replicates were used. These three factors, the amount of change, the amount of variation and the number of replicates are combined by the statistical test into one "p-value" which assesses the amount of evidence for differential expression. The traditional p-value cutoff of 0.05 for significance means that if you were to randomly sample two sets of replicates from the sample population 100 times, only 5 times would the value of the test statistic be larger than what you did see. This is a reasonable threshold when you only test a single gene, but in RNA-Seq you are measuring and testing thousands of genes at the same time. Therefore, some sort of adjustment of the p-values needs to be done. The most popular method is the False Discovery Rate method (Benjamimi & Hochberg, 1991), which adjusts the p-values so that the entire set of genes with values less than 0.05 is expected to contain 5% false-positives. Depending on the goals of your experiment, the FDR p-value threshold can be raised to gain more true-positives at the expense of a higher false-positive rate. However, keep in mind that the FDR correction depends on how many genes have low p-values compared with the number expected to to have low p-values and the lack of any genes with reasonably low p-values does not mean that no genes were truly changing.




#### _Working with Galaxy_ {#galaxy}

_A note about Galaxy_

If it is desirable to perform all processing in Galaxy, it should not be a problem for smaller experiments with a 1:1 comparisons between samples. For experiments with a large number of samples, and also for complex comparisons (e.g. 2x2 factorial design), Galaxy may not work as well; we instead recommend learning and using the command line tools and R/Bioconductor. However, Galaxy can be used to test parameter settings on a subset of the data, prior to switching over to command line for the whole analysis. Another issue of note, though Galaxy is quite good with data provenance,the latest versions of tools may not be available in Galaxy's Tool Shed (in particular those available in R/Bioconductor); In most cases Galaxy's excellent data provenance tracking will keep track of the versions used, so please make a note of the version numbers.  

The [Galaxy Tool Shed](https://toolshed.g2.bx.psu.edu) lists all tools currently available through Galaxy.  Tools available in Galaxy mentioned above include:



*   FASTQC
*   STAR
*   featureCounts
*   Trimmomatic, Trim_galore
*   Picard
*   edgeR, DESeq2
*   Kallisto
*   Salmon


## References {#references}

[https://genome.ucsc.edu/ENCODE/protocols/dataStandards/RNA_standards_v1_2011_May.pdf](https://genome.ucsc.edu/ENCODE/protocols/dataStandards/RNA_standards_v1_2011_May.pdf)
