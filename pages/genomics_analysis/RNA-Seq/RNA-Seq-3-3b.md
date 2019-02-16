---
title: RNA-Seq data processing and gene expression analysis
keywords: rna-seq, ngs
tags: [genomics_analysis]
last_updated: July 4, 2018
toc: false
sidebar: rna_sidebar
permalink: RNA-Seq-3-3b.html  
folder: genomics_analysis/RNA-Seq

complex_map: true
map_name: map_RNA-Seq_phase2a
box_number: 1

author_profile: true
authors:
 - Jenny_Zadeh
 - Radhika_Khetani
 - Jessica_Holmes
 - Chris_Fields
 - Meng-Chun Tseng
---
### _Step 2.3: Collecting and tabulating alignment stats_ {#step-2-3-collecting-and-tabulating-alignment-stats}


#### **Protocol 2**

Apart from `FASTQC`, other standard QC metrics that rely on an alignment are not available, such as Picard's tools, or a more complete assessment of read fates.  Salmon and kallisto both provide output that give basic overall statistics such and the number of reads mapped, and `MultiQC` can also summarize this information for all samples.

## Bibliography {#bibliography}

[^1]: Robinson, Mark D., and Alicia Oshlack. [A scaling normalization method for differential expression analysis of RNA-seq data.](https://doi.org/10.1186/gb-2010-11-3-r25) Genome biology 11.3 (2010): R25.

[^2]: Trapnell, C., Roberts, A., Goff, L., Pertea, G., Kim, D., Kelley, D. R., ... & Pachter, L. (2012). [Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks.](https://doi.org/10.1038/nprot.2012.016) Nature protocols, 7(3), 562.

[^3]: Love MI, Anders S, Kim V and Huber W. [RNA-Seq workflow: gene-level exploratory analysis and differential expression.](https://doi.org/10.12688/f1000research.7035.2) F1000Research 2016, 4:1070.

[^4]: Law CW, Alhamdoosh M, Su S et al. [RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR.](https://doi.org/10.12688/f1000research.9005.3) F1000Research 2018, 5:1408.

[^5]: Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L. (2016). [Near-optimal probabilistic RNA-seq quantification.](https://doi.org/10.1038/nbt.3519) Nature biotechnology, 34(5), 525.

[^6]: Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). [Salmon provides fast and bias-aware quantification of transcript expression.](https://doi.org/10.1038/nmeth.4197) Nature methods, 14(4), 417.

[^7]: Andrews, S. (2010). [FastQC: a quality control tool for high throughput sequence data.]( https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

[^8]: Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). [MultiQC: summarize analysis results for multiple tools and samples in a single report.](https://doi.org/10.1093/bioinformatics/btw354) Bioinformatics, 32(19), 3047-3048.

[^9]: Bolger, A. M., Lohse, M., & Usadel, B. (2014). [Trimmomatic: a flexible trimmer for Illumina sequence data.](https://doi.org/10.1093/bioinformatics/btu170) Bioinformatics, 30(15), 2114-2120.

[^10]: [Trim Galore.](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)

[^11]: [BBMap]( https://sourceforge.net/projects/bbmap/)

[^12]: Dodt, M., Roehr, J. T., Ahmed, R., & Dieterich, C. (2012). [FLEXBAR—flexible barcode and adapter processing for next-generation sequencing platforms.](https://doi.org/10.3390/biology1030895 Biology, 1(3), 895-905. )

[^13]: Gordon, A., & Hannon, G. (2010). [Fastx-toolkit. FASTQ/A short-reads pre-processing tools.]( http://hannonlab.cshl.edu/fastx_toolkit.) Unpublished

[^14]: (PRINSEQ) Schmieder, R., & Edwards, R. (2011). [Quality control and preprocessing of metagenomic datasets.](https://doi.org/10.1093/bioinformatics/btr026) Bioinformatics, 27(6), 863-864.

[^15]: Cox, M. P., Peterson, D. A., & Biggs, P. J. (2010). [SolexaQA: At-a-glance quality assessment of Illumina second-generation sequencing data.](https://doi.org/10.1186/1471-2105-11-485) BMC bioinformatics, 11(1), 485.

[^16]: Patro, R., Mount, S. M., & Kingsford, C. (2014). [Sailfish enables alignment-free isoform quantification from RNA-seq reads using lightweight algorithms.](https://doi.org/10.1038/nbt.2862) Nature biotechnology, 32(5), 462.

[^17]: Bray, N. L., Pimentel, H., Melsted, P., & Pachter, L. (2016). [Near-optimal probabilistic RNA-seq quantification.](https://doi.org/10.1038/nbt.3519) Nature biotechnology, 34(5), 525.

[^18]: Patro, R., Duggal, G., Love, M. I., Irizarry, R. A., & Kingsford, C. (2017). [Salmon provides fast and bias-aware quantification of transcript expression.](https://doi.org/10.1038/nmeth.4197) Nature methods, 14(4), 417.

[^19]: Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., ... & Gingeras, T. R. (2013). [STAR: ultrafast universal RNA-seq aligner.](https://doi.org/10.1093/bioinformatics/bts635) Bioinformatics, 29(1), 15-21.

[^20]: Kim, D., Langmead, B., & Salzberg, S. L. (2015). [HISAT: a fast spliced aligner with low memory requirements.](https://doi.org/10.1038/nmeth.3317) Nature methods, 12(4), 357.

[^21]: Wu, T. D., & Nacu, S. (2010). [Fast and SNP-tolerant detection of complex variants and splicing in short reads.](https://doi.org/10.1093/bioinformatics/btq057) Bioinformatics, 26(7), 873-881.

[^22]: Huang, S., Zhang, J., Li, R., Zhang, W., He, Z., Lam, T. W., ... & Yiu, S. M. (2011). [SOAPsplice: genome-wide ab initio detection of splice junctions from RNA-Seq data.](https://doi.org/10.3389/fgene.2011.00046) Frontiers in genetics, 2, 46.

[^23]: Anders, S., Pyl, P. T., & Huber, W. (2015). [HTSeq—a Python framework to work with high-throughput sequencing data.](https://doi.org/10.1093/bioinformatics/btu638) Bioinformatics, 31(2), 166-169.

[^24]: Liao, Y., Smyth, G. K., & Shi, W. (2013). [featureCounts: an efficient general purpose program for assigning sequence reads to genomic features.](https://doi.org/10.1093/bioinformatics/btt656) Bioinformatics, 30(7), 923-930.

[^25]: (Cufflinks) Trapnell, C., Williams, B. A., Pertea, G., Mortazavi, A., Kwan, G., Van Baren, M. J., ... & Pachter, L. (2010). [Transcript assembly and quantification by RNA-Seq reveals unannotated transcripts and isoform switching during cell differentiation.](https://doi.org/10.1038/nbt.1621) Nature biotechnology, 28(5), 511.

[^26]: Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., ... & Durbin, R. (2009). [The sequence alignment/map format and SAMtools.](https://doi.org/10.1093/bioinformatics/btp352) Bioinformatics, 25(16), 2078-2079.

[^27]: [Picard](http://broadinstitute.github.io/picard/)

[^28]: Gentleman, R. C., Carey, V. J., Bates, D. M., Bolstad, B., Dettling, M., Dudoit, S., ... & Hornik, K. (2004). [Bioconductor: open software development for computational biology and bioinformatics.](https://doi.org/10.1186/gb-2004-5-10-r80) Genome biology, 5(10), R80.

[^29]: Soneson, C., Love, M. I., & Robinson, M. D. (2015). [Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences.](https://doi.org/10.12688/f1000research.7563.2) F1000Research, 4.

[^30]: Gautier L, Cope L, Bolstad BM, Irizarry RA (2004). [affy—analysis of Affymetrix GeneChip data at the probe level.](https://doi.org/10.1093/bioinformatics/btg405) Bioinformatics, 20(3), 307–315. ISSN 1367-4803,

[^31]: Pertea, M., Pertea, G. M., Antonescu, C. M., Chang, T. C., Mendell, J. T., & Salzberg, S. L. (2015). [StringTie enables improved reconstruction of a transcriptome from RNA-seq reads.](https://doi.org/10.1038/nbt.3122) Nature biotechnology, 33(3), 290.

[^32]: Dillies, M. A., Rau, A., Aubert, J., Hennequet-Antier, C., Jeanmougin, M., Servant, N., ... & Guernec, G. (2013). [A comprehensive evaluation of normalization methods for Illumina high-throughput RNA sequencing data analysis.](https://doi.org/10.1093/bib/bbs046) Briefings in bioinformatics, 14(6), 671-683.

[^33]: Anders S, Huber W (2010). [Differential expression analysis for sequence count data.](https://doi.org/ 10.1186/gb-2010-11-10-r106) Genome Biology, 11, R106.

[^34]: Love M.I., Huber W., Anders S. (2014). [Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.](https://doi.org/10.1186/s13059-014-0550-8) Genome Biology, 15, 550.

[^35]: Robinson MD, McCarthy DJ, Smyth GK (2010). [edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.](https://doi.org/10.1093/bioinformatics/btp616) Bioinformatics, 26(1), 139-140.

[^36]: McCarthy, J. D, Chen, Yunshun, Smyth, K. G (2012). [Differential expression analysis of multifactor RNA-Seq experiments with respect to biological variation.](https://doi.org/10.1093/nar/gks042) Nucleic Acids Research, 40(10), 4288-4297.

[^37]: Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the royal statistical society. Series B (Methodological), 289-300.

[^38]: Afgan, E., Baker, D., Batut, B., Van Den Beek, M., Bouvier, D., Čech, M., ... & Guerler, A. (2018). [The Galaxy platform for accessible, reproducible and collaborative biomedical analyses: 2018 update.](https://doi.org/10.1093/nar/gky379) Nucleic acids research, 46(W1), W537-W544.

[^39]: Leek, J.T., & Peng, R. D. (2015). [Reproducible research can still be wrong.](https://doi.org/10.1073/pnas.1421412111) Proceedings of the National Academy of Sciences 112 (6) 1645-1646

[^40]: [Reproducibility in Science: A Guide to enhancing reproducibility in scientific results and writing.](http://ropensci.github.io/reproducibility-guide/)

[^41]: Peng, R. D. (2011). [Reproducible Research in Computational Science.](https://doi.org/10.1126/science.1213847) Science Vol. 334, Issue 6060, pp. 1226-1227.

[^42]: Ashburner, M., Ball, C. A., Blake, J. A., Botstein, D., Butler, H., Cherry, J. M., ... & Harris, M. A. (2000). [Gene Ontology: tool for the unification of biology.](https://doi.org/10.1038/75556) Nature genetics, 25(1), 25.

[^43]: Kanehisa, M., & Goto, S. (2000). KEGG: kyoto encyclopedia of genes and genomes. Nucleic acids research, 28(1), 27-30.

[^limma]: Ritchie, ME, Phipson, B, Wu, D, Hu, Y, Law, CW, Shi, W, and Smyth, GK (2015). [limma powers differential expression analyses for RNA-sequencing and microarray studies](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4402510/). Nucleic Acids Research 43(7), e47.

[^conessa]: Conesa, A., Madrigal, P., Tarazona S., ... & Mortazavi, A. 2016. [A survey of best practices for RNA-Seq data analysis.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8) Genome Biology 17:13. https://doi.org/10.1186/s13059-016-0881-8


[//]: <> (These are common abbreviations in the page.)
*[H3ABioNet]: The Bioinformatics Network within the H3Africa Consortium
*[FASTQ]: Standard format of raw sequence data. Quality scores assigned in the FASTQ files represent the probability that a certain base was called incorrectly. These scores are encoded in various ways and it is important to know the type of encoding for a given FASTQ file.
*[TPM]:Transcripts Per Million
*[RPKM]:Reads Per Kilobase per Million mapped reads (for SE data)
*[FPKM]:Fragments Per Kilobase per Million mapped reads (for PE data)
*[CPM]:Counts Per Million, relative to the total number of reads in the sample
*[SE]: single-end
*[PE]: paired-end
*[FASTQ]: The standard format of raw sequence data
*[TMM]: A type of normalization and is an acronym for "Trimmed Mean of Ms
*[Adaptor]: Artificial pieces of DNA introduced prior to sequencing to ensure that the DNA fragment being sequenced attaches to the sequencing flow cell.
*[adaptor]: Artificial pieces of DNA introduced prior to sequencing to ensure that the DNA fragment being sequenced attaches to the sequencing flow cell.
*[EM]: Expectation maximization
*[FDR]: False Discovery Rate
*[H3ABioNet]: The Bioinformatics Network within the H3Africa Consortium
*[NGS]: Next Generation sequencing
*[BAM]: Binary Alignment Format
*[SAM]: Sequence Alignment Format
*[VCF]: Variant Calling Format
