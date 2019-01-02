---
title: GWAS data processing 
keywords: GWAS, population structure, association testing 
tags: [genomics_analysis]
last_updated: December 23, 2018 
summary: "This SOP is intended to layout best practices for GWAS data processing, especially for groups undertaking H3ABioNet accreditation exercises"
sidebar: gwas_sidebar
permalink: GWAS.html  
author_profile: true
authors:
 - Scott_Hazelhurst
 - Segun_Fatumo
---

## Introduction {#introduction}

GWAS is a key workflow of the H3A. Our recommended approach to SOP is shown below. Note, this is not an algorithm -- there are multiple ways of doing QC which are good. Also, what is shown below is not necessarily linear. For example, some of the analyses at the plate level might only be possible after doing an initial QC and then coming back to look at plate and batch effects.

While this SOP is designed for bioinformaticists doing the GWAS analysis once the data is ready, we recommend that it be studied thoroughly to help plan the entire process. This GWAS is written as a general guide for bioinformaticists, and in particular to assist groups undertaking H3ABioNet accreditation, both to prepare for and do the accreditation exercisde.

*Disclaimer: * Although we hope this SP is educational and will help groups learn to do GWAS, it is not meant as a tutorial or a complete checklist. *Groups undertaking accreditation should realise that accreditation is undertaking by an independent international evaluation committee, and while they will have regard to the SOP, they make their decision at their own discretion.*



### Tool and file format {#tools-file-format}
This SOP assumes that is that much data is in the binary format used by the PLINK software suite. However, this is not mean to say that PLINK is the only or even best tool to be used.  The PLINK binary format (hereafter referred to as bped) encodes a dataset as a set of three files, with the following suffixes to their names:

* .bed: this is the binary genotype data, stored as 2 bits per genotype per sample
* .bim: a text file containing marker information, one line per marker, in genomic order
* .fam: a text file containing sample pedigree information


### Glossary of associated terms and jargon  {#glossary-of-associated-terms-and-jargon}


## Genotype calling {#calling}

The raw output of the genotyping process are image files -- for Illumina products, these are IDAT files. Each image gives the data for one SNP, and each individual will be have one dot on the image. For haploid SNPs, typically the dots will visually cluster into three groups: those homoazygous for the reference allele, those heterozygous, those homozygous for the alternate allele. Note that the image files just describe the positions of each individual for each SNP: usually the clusters are visually obviously, and the process of *calling* is to rigorously put each call in one of three clusters.

Typically most SNPs will have well behaved clusters and most individuals will clearly associate with one of the clusters, but there are always errors, anomalies and some borderline cases. Different clustering algorithms will produce different results. Typically, genotyping centres will provide a cluster report which describes for each SNP and each individual what the calls are. While, well reputed centres will produce high quality calls, it is always a  good idea to try alternate approaches. Well known calling tools are Illumina's GenomeStudion and crlmm

Our SOP does not yet take this process forward.

As is noted later, part of the QC will require coming back to the raw data.

## Sources of error {#errors}

Error can occur at multiple places, and these errors can easily introduce false signals or mask real signals

* poor record keeping in the field, in transport, in storage, DNA extraction, plating and shipping
* sample contamination
* DNA extraction
* Problems with plates

Effort should be made to detect and eliminate any such problems.

Problems may occur randomly but there may also be due to non-random effects. For examples, a problem with a particular consignment of reagents may make two or three days' worth of DNA extraction unreliable. The better record keeping there is of samples, the easier it is to work out problems. Particular considerations are.

### Plate Level errors {#plate-errors}

Currently (2019), samples are processed in batches of 96 wells (8x12), with one sample per well.  Comparative QC should be done across the plates according to the QC parameters described below. If an entire plate is found to be found to behave poorly, it should be removed. 

### Batch level {#batch-errors}

A batch is a set of samples that are processed together. Ideally for consistency, all samples would be processed at the same   time (collected, DNA extraction, genotyping) but for any real study this is not practical. As an example, in one H3A project, the samples were divided into four roughly equal batches, each batch corresponding to one campaign of DNA extractions. Batches 1-3 were shipped together to the genotyping centre and genotyped, and about 6 months later Batch 4 was shipped and genotyped. In this case, the group had to compare the four batches with each other and compare batch 1-3 with batch 4.

There may be different ways of dividing samples into batches (e.g. by recruitment stratum, when DNA extracted etc).

### Site level {#site-errors}

If participants are recruited in different places, differences may occur due to the way in which samples are handled at these places.

### Case/control or phenotype  differences {#case-control-errors}

Of course, the whole purpose of a GWAS is to discover genetic differences between cases and controls. However, in a well designed study it very highly unlikely that the condition or disease of interest will cause observable differences at the genomic level. If there are overall genomic differences then this is most probably due to the way in which participants were sampled or the different processes in which DNA was extracted or genotyped. For example, a study might recruit 3000 cases and use results from a matched cohort for controls. As different people/labs might be responsible for collection and processing, even if indentical arrays are used and genotyping is done at the same centre, there may be systemic errors.

One useful phenotype to use for QC is sex. For many studies, it would not be expected that in the autosomal SNPs that there are differences at the genomic ldevel. 

### Population structure {#population-structure}

As discussed later, population structure is a potential confounding factor. This may also make QC more difficult. In a highly homogeneous sample, one would not expect there to be significant genomic level difference between plates (e.g., when doing a PC analysis). However, in a large multi-country project with diverse ancestries, this would not be unexpected. This is especially because unless a deliberate effort is made to randomise (which may be not be practical) it is likely that bacthes of samples from the same site are likely to be processed close to each other.

A particular concern in a study with population structure is if the number of cases and controls is not balanced across different sub-groups and even more so when differences in both environmental and genomic background might lead to different outcomes



### Random errors {#random-errors}

Many errors will be due to random effects at either the SNP or individual level. The systematic errors above are often detected by comparing groups of things with each other, while random effects are usually detected using some experiment-wide cut-off


## Calling quality {#quality}

For Illumina data, the key quality measure is called _GenCall_. Once clustering is done (recall, clustering is done per SNP, using the (x,y) coordinates of all the individuals for that SNP), a GenCall score is given for each sample/person of that SNP which a confidence score of the correctness of the call. This confidence score is  a function of the the clustering algorithm used. From this is derived
* The _call rate_ or _call frequency_.  The person doing the calling picks a lower value for an acceptable _GenCall_value, and any calls lower than this are rejected.
* The 10%_GenCall value: The GenCall value of the SNP at the 10% percentile of GenCall value (e.g, the SNP/Indvidual for which 90% of other individuals have a _better_ GenCall score). For a well genotyped SNP or individual this should be a very high quality absolute value even though it is relatively poor.


These can be measured per sample or per SNP.  Good graphs to look at are _call rate_ versus sample number, and _call rate_ versus 10%_GenCall value, both of which may be used to identify poor performing individuals.

If biological replicates have been included in your study (a good idea), the concordance between these replicates should be examined as these will give good empirical evidence of thquality of genotyping. Similarly, if there are big differences in genoyping for biological replicates, these may indicate problems in sample handling.

## Quality control {#qc}

The following should be checked in the QC process:

* Missingness 
  * Missingness at the SNP-level. SNPs with high missingness should be removed. Usually this is set at the 1-2% level. (PLINK `--geno` option)
  * Missingness at the individual level. Individuals with high missingness should be removed. Usually this is set at the 1-2% level. (PLINK `--mind` option)

   Of course, the two are inter-related so a small benefit may be gained by iteratively removing badly behaved SNPs poor individuals. Note that if a significant number of poorly genotyped samples are found, it may be desirable to remove the individuals and re-call the genotypes since poorly performing samples may significantly affect how clustering is done.

* Sex concordance. Sex concordance is typically used as a proxy that there hasn't been swapping of samples. We check that the individual's biological sex as recorded in our meta-data about the individual matches what the individual's genotype says, typically as measured by the computing the in-breeding co-efficient of the X-chromosome. PLINK computes this with the F-statistic. By default, individuals with an F statistic less than 0.2 are regarded as female, and individuals with F-statistic greater than 0.8 are male.  Of course the the scores are 0.2 and and 0.8 are arbitrary so care must be taken. In large samples, there may be people with unusual F-statistics for a variety of reasons. However, if a project has an individual recorded as male and their F-statistic is 0, then the most likely explanation is that there has been sample mishandling at some stage and that results we have are not for the individual we think it is. Any individuals who fail the F-statistic like this should be removed. This test is a proxy for testing correct labellings and record keeping: it cannot detect swaps where M/M or F/F swaps have taken place. If  _p_ sex concordance errors are found, it is likely that there are _p_ other swap erros that have not been detected. Pay attention to the distribution of errors -- are they randomly distributed or bunched in batches or groups. In some extreme cases  it may be necessary to remove an entire group

* Population structure concordance. Using population structure approaches (discussed below), individuals should cluster appropriately.  Some insight and judicious judgement is required depending on the sample chosen and what is known about the individuals. The more meta-data, the more accurate our testing can be done. But for example in a multi-centre trial, if an individual has told us that all four of her grand-parents were Zulu-speaking and identified as Zulu, but genetically the individual clusters with other individuals from a West African site then (although there are all sorts of interesting life stories) there is a strong possiblity that this is a sampling error. As discussed above, population structure differences between batches or plates may also indicate errors, but in a heterogeneous study or multi-centre study this may be hard to test for. Individuals whose PC position falls far from the centroid of the group are likely a result of mislabelling or even genotyping problems.

* Minor allele frequency. The ability to reliably detect and quantify variation is limited by the sample size and error rate of genotyping. A minimum minor allele frequency must be imposed since with very low MAF even a low error rate in genotyping can cause cases to be appear different to controls.

* Difference in MAFs between sub-groups of the study (either tested directly using chi^2 or similar test or PC analysis) may also indicate batch effects of different sorts. As mentioned above this is complicated by population structure.

* Deviation from expected heterozygosity. The theory of Hardy-Weinberg equilibrium indicates that if the MAF of a SNP is _f_, the  distribution of [homozygous minor allele, heterozygous, homozygous major alleles] should be [f^2, (1-f)f, (1-f)^2]. Deviation from this may imply problems with the array. Batch effects of different sorts should be considered.  Deviations from heterozygosity can be detected at the SNP and sample level. However, there are complexities.
  * in a recently admixed population, deviations from HWE are likely;
  * SNPs implicated in a disease may well be out of deviation from HWE, especially in cases;
  * Typically we test for deviation from HWE by performing  a chi^2 or similar test for each SNP. Since we may be testing hundreds of thousands or millions of SNPs we may need to adjust for multiple testing.

  Thus there are reasons where deviation from HWE may not be a sign of error; but deviation from HWE may also indicate a problem with the assay. This makes setting the correct cut-off value for HWE testing a challenge.

  PLINK's `--hardy` and `--hwe` options will be very useful.

* Relatedness. An investigation of relatedness between the samples is important for two reasons. There are several ways of measuring relatedness and you should consider the most appopriate for your study. One way of measuring relatedness is the PI-hat score computed by PLINK. There are two reasons for testing for relatedness (the discussion below uses PI-hat  but similar considerations apply to other approaches, _mutatis _mutandis).
  * If there are pairs of individuals with a PI-hat score close to 1 then this indicates (a) a known biological replicate; or  (b) a sample handing error; (c) the same individual being recruited more than once into the study. In case (a) one of the pair should be removed; in case (b)  and (c) it is likely that both pairs will need to be removed. PI-hat values significantly above 0.5 but less than 1 are highly unlikely biologically and if you find any such pairs these are probably caused by contamination of some sort
   * The level of relatedness in the study, which will determine analytic approach. For family based studies a high level of relatedness is expected and it is important to check that the genotype-computed relatedness matches what is expected in the study. For most GWASes samples are  ideally randomly recruited and some methods rely on low levels of relatedness. For samples recruited in large urban areas, random sampling should result in low levels of relatedness and if there are a significant number of related pairs, there may be a problem in experimental design and execution. In smaller regions, depending on region histry and cultural practices, a relatively high number of related pairs is possible.

    If methods such as PLINK's chi^2 or linear or logistic regression are used, a cut-off value of relatedness should be used and one member of each pair with a PI-hat value higher than this cut-off should be removed. It is hard to say what a correct cut-off value is -- the literature has studies where PI-hat values ranging from very low (e.g. 0.01) to 0.2 and the choice may be an expedient one. Recall that PI-hat of 0.125 corresponds to first cousin relatedness and 0.25 to grandparent-grandchild relatedness. There are several studies that use 0.18 as a cut-off. Higher than this may be problematic

  Alternatively, a method that can handle relatedness should be used. For example, mixed models approaches such as GEMMA and BOLTLMM can handle related samples if used appropriately. 

  Pairwise relatedness can be computed using PLINK’s `­­genome`xs option

* Samples with chromosomal abnormalities
  Chromosomal abnormalities such as aneuploidy or long stretches of homozygosity should be tested for and the causes identified. These cases should be examined with a geneticist to determined possible causes and decide whether the sample should be removed.


### Strand errors {#strand-errors}

As a byproduct of the genotyping technology, and the annotation data that accompanies the chips, SNP data from both Illumina and Affymetrix platforms may be reported as the allele on either the “forward” or “reverse” strand. Although the information to orient these calls properly is contained within the annotation data, conversion from the native format to PLINK format can be tricky, and errors will not be evident in “ambiguous” A­T/G­C SNPs. Therefore it is prudent to check that the alleles reported in the bim file match the known alleles in dbSNP or on the Ensembl genome browser. Another easy check is to identify the control samples that are often included in the dataset (eg samples from HapMap) and compare these to their known genotypes.

This is a particular problem when you will be merging your data with other sets.


### Build errors {#build-errors}


### Merging data from different genotype experiments {#merging-data}



## Association Study {#assoc}


The following document describes some recommended steps for processing genotype data
from a genome­wide association study, as will be generated by many of the H3Africa projects. Although it is written with data from SNP genotyping arrays in mind, many of the steps also apply to full genome sequence and exome sequencing data. Many of the tests can be performed with the PLINK software suite, but where other software is required this will be indicated.

Some steps below may be repeat QC steps mentioned above. If QC has been done thoroughly above then they may be omitted now.






### Population stratification {#strafication}
Population stratification refers to structure within the sample group that is the result of systematic genetic differences between individuals that correlate with the phenotypic data. This could result in allele frequency differences that are due to ancestral proportion differences between cases and controls being mistaken for an association with the phenotype.

Population structure be can used
* to detect errors (as describe above)
* remove individual outliers -- these are not _errors_ but individuals whose genetic background are so different from others in the group that they may confound result
* adjust for genetic differences between and within groups

Allele frequency based methods, such as STRUCTURE or Admixture can be used.. Add known reference populations such as 1000 Genomes population. Determine most likely _K_ and compute the ancestral components of each individual


PCA based methods, such as Eigensoft. Eigensoft’s smartpca, and examine the top _n_ eigenvectors. If any of these eigenvectors display segmentation based on a potential confounder, these factors should be examine further, and possibly corrected for. Individuals that differ significantly from the rest in terms of cluster membership should be investigated for possible removal. A common threshold in Eigensoft is to remove individuals that are outliers by more than 6 standard deviations.

Both these methods should be used on the data in an exploratory fashion, to try identify possible confounding factors *and* to undertstand how to handle population structure in the GWAS.

Different types and levels of structure can be observed. At one extreme, sampling from a very homogeneous group, relatively isolated from others will lead to a very tighly cluster observable in a PCA. In such a case, population structure may be ignored. In an admixed group formed 10 generations or so ago in which random mating has taken place, a single cluster may be observed that is fuzzy and loose  may be observed. In a GWAS recruited from multiple sites, multiple distinct clusters may be observed with few intermediate individuals, but more complex patterns can be seen, particularly in cosmpolitan, African regions. Managing population by adjusting
* GC control
* PCA analysis: including the right number of eigenvectors in the analysis
* Using mixed models approach
* Meta-analysis
* Come combination of the above




### Autosomal, X, Y and MT SNPs {#snps}

Many GWASes only consider the haploid chromosomes and so only autosomal SNP are considered.

 
### Association testing for single locus  {#single-locus}

The most common tests are single locus testing. Most GWAS tools provide multiple tests and insight into your biological question is needed. For example, if it is known that the disease being studied is recessive or dominant the correct variation on the association test can be used

Some of the approaches are 
* chi^2, logistic regression (case control)
* linear regression (quantitative study)
* mixed models approaches (e.g., GEMMA, BOLTLMM) which 


### Multiple testing {#multiple-testing}

Multiple testing should be dealt with, for example using Bonferroni correction or permuation testing.

### Covariates {#covariates}

Covariates include population structure, sex, smoking , age, socio-economic status. How to handle covariates is complex. Crudely there are two types of covariates:
* confounding covariates where not including the co-variate may cause false signals. Population structure is a good example of  this. These are often covariates which mediate the genomic background (like population structure)
* covariates where not taking them into account will hide signals

### Imputation {#imputation}

Imputation can be done in two modes -- first, one can use imputation just to fill in missing genotypes. Second, imputation can be done to impute values for as many postions as pssoble

### Multi-locus testing {#multi-locus-testing}

### Replication {#replication}
The relatively low power of most GWAS designs means that a substantial number of false positives can be expected to be generated, so validation via replication in an independent population is an important part of most studies.
A common strategy is to genotype a subset of the samples on a high coverage array, and then follow up with targeted genotyping of markers that appear interesting, on a larger number of samples

### Meta analyses {#meta-analyses}


## References {#references}


* Highland et al, 2018. Quality Control Analysis of Add Health GWAS Data. https://www.cpc.unc.edu/projects/addhealth/documentation/guides/AH_GWAS_QC.pdf    (Example of a good QC)


* Illumina Inc, 2014. Infinium Genotyping Data Analysis. https://www.illumina.com/Documents/products/technotes//technote_infinium_genotyping_data_analysis.pdf



* Turner et al. [Quality Control Procedures for Genome Wide Association Studies](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3066182/). _Current Protoclols in Human Genetics_, 2011,

* Zeggini and Morris 2010. Analysis of Complex Disease Association Studies: A Practical Guide. Academic Press.
