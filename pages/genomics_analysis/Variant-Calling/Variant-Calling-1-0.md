---
title: Variant calling in human whole genome/exome sequencing data
keywords: wgs, wes, ngs
tags: [genomics_analysis]
last_updated: Fall, 2018

sidebar: varcall_sidebar
toc: false
permalink: Variant-Calling-1-0.html
folder: genomics_analysis/Variant-Calling
author_profile: true
authors:
 - Azza_Ahmed
 - Matthew_Weber
 - Faisal_Fadlelmola
 - Luidmila_Mainzer
---


## Introduction {#introduction}

This document briefly outlines the essential steps in the process of making genetic variant calls, and recommends tools that have gained community acceptance for this purpose. It is assumed that the purpose of the study is to detect short germline or somatic variants in a single sample. Recommended coverage for acceptable quality of calls in a research setting is around 30-50x for whole genome and 70-100x for exome sequencing, but lower coverage is discussed as well.

The procedures outlined below are recommendations to the H3ABioNet groups planning to do variant calling on human genome data, and are not meant to be prescriptive. Our goal is to help the groups set up their procedures and workflows, and to provide an overview of the main steps involved and the tools that can be used to implement them. For optimizing a workflow or an individual analysis step, the reader is referred to [^1]<sup>,</sup>[^2]<sup>,</sup>[^3].

## Bibliography {#bibliography}


[^1]: Heldenbrand, J. R. et al. [Performance benchmarking of GATK3.8 and GATK4](https://www.biorxiv.org/content/10.1101/348565v1). BioRxiv (2018). doi:10.1101/348565

[^2]: Laurie, S. et al. [From Wet-Lab to Variations: Concordance and Speed of Bioinformatics Pipelines for Whole Genome and Whole Exome Sequencing](https://onlinelibrary.wiley.com/doi/full/10.1002/humu.23114). Hum. Mutat. 37, 1263–1271 (2016).

[^3]: Hwang, S., Kim, E., Lee, I. & Marcotte, E. M. [Systematic comparison of variant calling pipelines using gold standard personal exome variants](https://www.nature.com/articles/srep17875). Sci. Rep. 5, 17875 (2015).

[^4]: GATK Dictionary. at <https://software.broadinstitute.org/gatk/documentation/topic?name=dictionary>

[^5]: Myllykangas, S., Buenrostro, J. & Ji, H. P. Overview of sequencing technology platforms. InBioinformatics for High Throughput Sequencing 11–25 (Springer New York, 2012). doi:10.1007/978-1-4614-0782-9_2

[^6]: Schiemer, J. Illumina TruSeq DNA Adapters De-Mystifie available from <http://tucf-genomics.tufts.edu/documents/protocols/TUCF_Understanding_Illumina_TruSeq_Adapters.pdf> .

[^7]: Goodwin, S., McPherson, J. D. & McCombie, W. R. [Coming of age: ten years of next-generation sequencing technologies](https://www.nature.com/articles/nrg.2016.49). Nat. Rev. Genet. 17, 333–351 (2016).

[^8]: Head, S. R. et al. [Library construction for next-generation sequencing: overviews and challenges](https://www.future-science.com/doi/10.2144/000114133). BioTechniques 56, 61–4, 66, 68, passim (2014).

[^9]: Regier, A. A. et al. [Functional equivalence of genome sequencing analysis pipelines enables harmonized variant calling across human genetics projects](https://www.nature.com/articles/s41467-018-06159-4). BioRxiv (2018). doi:10.1101/269316

[^10]: broadinstitute/gatk: Official code repository for GATK versions 4 and up. at <https://github.com/broadinstitute/gatk>

[^11]: Pabinger, S. et al. [A survey of tools for variant analysis of next-generation genome sequencing data](https://academic.oup.com/bib/article/15/2/256/210976). Brief. Bioinformatics 15, 256–278 (2014).

[^12]: Nielsen, R., Paul, J. S., Albrechtsen, A. & Song, Y. S. [Genotype and SNP calling from next-generation sequencing data](https://www.nature.com/articles/nrg2986). Nat. Rev. Genet. 12, 443–451 (2011).

[^13]: Cock, P. J. A., Fields, C. J., Goto, N., Heuer, M. L. & Rice, P. M. [The Sanger FASTQ file format for sequences with quality scores, and the Solexa/Illumina FASTQ variants](https://academic.oup.com/nar/article/38/6/1767/3112533). Nucleic Acids Res. 38, 1767–1771 (2010).

[^14]: Bolger, A. M., Lohse, M. & Usadel, B. [Trimmomatic: a flexible trimmer for Illumina sequence data](https://academic.oup.com/bioinformatics/article/30/15/2114/2390096). Bioinformatics 30, 2114–2120 (2014).

[^15]: Dodt, M., Roehr, J. T., Ahmed, R. & Dieterich, C. [FLEXBAR-Flexible Barcode and Adapter Processing for Next-Generation Sequencing Platforms](https://www.mdpi.com/2079-7737/1/3/895). Biology (Basel) 1, 895–905 (2012).


[^16]: Tools to remove adapter sequences from next-generation sequencing data, Genomics Gateway. at <http://bioscholar.com/genomics/tools-remove-adapter-sequences-next-generation-sequencing-data/>

[^17]: Adapter trimming software tools, WGS analysis - OMICtools. at <https://omictools.com/adapter-trimming-category>

[^18]: Del Fabbro, C., Scalabrin, S., Morgante, M. & Giorgi, F. M. [An extensive evaluation of read trimming effects on Illumina NGS data analysis](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0085024). PLoS ONE 8, e85024 (2013).

[^19]: Babraham Bioinformatics - FastQC A Quality Control tool for High Throughput Sequence Data. at <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>

[^20]: Schmieder, R. & Edwards, R. [Quality control and preprocessing of metagenomic datasets](https://academic.oup.com/bioinformatics/article/27/6/863/236283). Bioinformatics 27, 863–864 (2011).

[^21]: Cox, M. P., Peterson, D. A. & Biggs, P. J. [SolexaQA: At-a-glance quality assessment of Illumina second-generation sequencing data](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-485). BMC Bioinformatics 11, 485 (2010).

[^22]: GATK Doc Article #1213. at <https://software.broadinstitute.org/gatk/documentation/article.php?id=1213>

[^23]: Li, H. [Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM](https://arxiv.org/abs/1303.3997). (2013).

[^24]: Langmead, B. & Salzberg, S. L. [Fast gapped-read alignment with Bowtie 2](https://www.nature.com/articles/nmeth.1923). Nat. Methods 9, 357–359 (2012).

[^25]: Li, H. & Homer, N. [A survey of sequence alignment algorithms for next-generation sequencing](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r99). Brief. Bioinformatics 11, 473–483 (2010).

[^26]: Thankaswamy-Kosalai, S., Sen, P. & Nookaew, I. [Evaluation and assessment of read-mapping by multiple next-generation sequencing aligners based on genome-wide characteristics](https://www.ncbi.nlm.nih.gov/pubmed/28286147). Genomics 109, 186–191 (2017).

[^27]: Li, H. et al. [The Sequence Alignment/Map format and SAMtools](https://academic.oup.com/bioinformatics/article/25/16/2078/204688). Bioinformatics 25, 2078–2079 (2009).

[^28]: How PCR duplicates arise in next-generation sequencing. at <http://www.cureffi.org/2012/12/11/how-pcr-duplicates-arise-in-next-generation-sequencing/>

[^29]: Faust, G. G. & Hall, I. M. [SAMBLASTER: fast duplicate marking and structural variant read extraction](https://academic.oup.com/bioinformatics/article/30/17/2503/2748175). Bioinformatics 30, 2503–2505 (2014).

[^30]: Tarasov, A., Vilella, A. J., Cuppen, E., Nijman, I. J. & Prins, P. [Sambamba: fast processing of NGS alignment formats](https://academic.oup.com/bioinformatics/article/31/12/2032/214758). Bioinformatics 31, 2032–2034 (2015).

[^31]: NOVOCRAFT TECHNOLOGIES SDN BHD. Novocraft. NOVOCRAFT TECHNOLOGIES SDN BHD at <http://www.novocraft.com/>

[^32]: Picard Tools - By Broad Institute. at <https://broadinstitute.github.io/picard/>

[^33]: DePristo, M. A. et al. [A framework for variation discovery and genotyping using next-generation DNA sequencing data](https://www.nature.com/articles/ng.806). Nat. Genet. 43, 491–498 (2011).

[^34]: Albers, C. A. et al. [Dindel: accurate indel calls from short-read data](https://genome.cshlp.org/content/21/6/961.long). Genome Res. 21, 961–973 (2011).

[^35]: Homer, N. & Nelson, S. F. [Improved variant discovery through local re-alignment of short-read next-generation sequencing data using SRMA](https://academic.oup.com/bioinformatics/article/31/12/2032/214758). Genome Biol. 11, R99 (2010).

[^36]: Van der Auwera, G. A. et al. [From FastQ data to high confidence variant calls: the Genome Analysis Toolkit best practices pipeline](https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/0471250953.bi1110s43). Curr. Protoc. Bioinformatics 11, 11.10.1-11.10.33 (2013).

[^37]: Rimmer, A. et al. [Integrating mapping-, assembly- and haplotype-based approaches for calling variants in clinical sequencing applications](https://www.nature.com/articles/ng.3036). Nat. Genet. 46, 912–918 (2014).

[^38]: Garrison, E. & Marth, G. [Haplotype-based variant detection from short-read sequencing](https://arxiv.org/abs/1207.3907). arXiv (2012).

[^39]: Poplin, R. et al. [Scaling accurate genetic variant discovery to tens of thousands of samples](https://www.biorxiv.org/content/10.1101/201178v3). BioRxiv (2017). 

[^40]: Ebbert, M. T. W. et al. [Evaluating the necessity of PCR duplicate removal from next-generation sequencing data and a comparison of approaches](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1097-3). BMC Bioinformatics 17 Suppl 7, 239 (2016).

[^41]: Olson, N. D. et al. [Best practices for evaluating single nucleotide variant calling methods for microbial genomics](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4493402/). Front. Genet. 6, 235 (2015).

[^42]: Cabanski, C. R. et al. [ReQON: a Bioconductor package for recalibrating quality scores from next-generation sequencing data](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-221). BMC Bioinformatics 13, 221 (2012).

[^43]: Bansal, V. [A statistical method for the detection of variants from next-generation resequencing of DNA pools](https://academic.oup.com/bioinformatics/article/26/12/i318/285976). Bioinformatics 26, i318-24 (2010).

[^44]: Cibulskis, K. et al. [Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples](https://www.nature.com/articles/nbt.2514). Nat. Biotechnol. 31, 213–219 (2013).

[^45]: Danecek, P. et al. [The variant call format and VCFtools](https://academic.oup.com/bioinformatics/article/27/15/2156/402296). Bioinformatics 27, 2156–2158 (2011).

[^46]: Li, H. [Toward better understanding of artifacts in variant calling from high-coverage samples](https://academic.oup.com/bioinformatics/article/30/20/2843/2422145). Bioinformatics 30, 2843–2851 (2014).

[^47]: hail-is/hail: Scalable genomic data analysis. at <https://github.com/hail-is/hail>

[^48]: Paila, U., Chapman, B. A., Kirchner, R. & Quinlan, A. R. [GEMINI: integrative exploration of genetic variation and genome annotations](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003153). PLoS Comput. Biol. 9, e1003153 (2013).

[^49]: Peng, G. et al. [Rare variant detection using family-based sequencing analysis](https://www.pnas.org/content/110/10/3985). Proc Natl Acad Sci USA 110, 3985–3990 (2013).

[^50]: Boyle, A. P. et al. [Annotation of functional variation in personal genomes using RegulomeDB](https://genome.cshlp.org/content/22/9/1790.long). Genome Res. 22, 1790–1797 (2012).

[^51]: What's in the resource bundle and how can I get it? — GATK-Forum. at <https://gatkforums.broadinstitute.org/gatk/discussion/1213/whats-in-the-resource-bundle-and-how-can-i-get-it>

[^52]: Afgan, E. et al. [The Galaxy platform for accessible, reproducible and collaborative biomedical analyses: 2018 update](https://academic.oup.com/nar/article/46/W1/W537/5001157). Nucleic Acids Res. 46, W537–W544 (2018).

[^53]: Stephens, Z. D. et al. [Simulating Next-Generation Sequencing Datasets from Empirical Mutation and Sequencing Models](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0167047). PLoS ONE 11, e0167047 (2016).



[//]: <> (These are common abbreviations in the page.)
*[H3ABioNet]: The Bioinformatics Network within the H3Africa Consortium
*[FASTQ]: Standard format of raw sequence data. Quality scores assigned in the FASTQ files represent the probability that a certain base was called incorrectly. These scores are encoded in various ways and it is important to know the type of encoding for a given FASTQ file.
*[Adapter]: Short nucleotide sequences added on to the ends of the DNA fragments that are to be sequenced.
*[Lane]: The basic machine unit for sequencing. The lane reflects the basic independent run of an NGS machine. For Illumina machines, this is the physical sequencing lane.
*[Library]: A unit of DNA preparation that at some point is physically pooled together.  Multiple lanes can be run from aliquots from the same library. The DNA library is the natural unit that is being sequenced. For example, if the library has limited complexity, then many sequences are duplicated and will result in a high duplication rate across lanes
*[NGS]: Next Generation sequencing
*[WGS]: Whole Genome Sequencing
*[sample]: A single individual, such as human CEPH NA12878. Multiple libraries with different properties can be constructed from the original sample DNA source. Here we treat samples as independent individuals whose genome sequence we are attempting to determine. From this perspective, tumor/normal samples are different despite coming from the same individual.
*[SNV]: Single nucleotide variant, whether in a non-coding region or in a coding regoin such that it is synonymous or nonsynonymous (which is then classified as missense or nonsense variant).
*[Functional equivalence]:  Specifications intended to eliminate batch effects and promote data interoperability by standardizing pipeline implementations: used tools, versions of these tools, and versions of reference genomic files. Large genomic databases, like gnomAD and TOPmed are being processed by pipelines adhering to these specifications.
*[Functional Equivalence]:  Specifications intended to eliminate batch effects and promote data interoperability by standardizing pipeline implementations: used tools, versions of these tools, and versions of reference genomic files. Large genomic databases, like gnomAD and TOPmed are being processed by pipelines adhering to these specifications.
*[GATK]: The Genome Analysis Toolkit, a commonly used framework and toolbox for Variant Calling
*[BAM]: Binary Alignment Format
*[VCF]: Variant Calling Format
