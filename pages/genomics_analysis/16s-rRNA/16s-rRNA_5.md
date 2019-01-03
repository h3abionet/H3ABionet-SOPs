---
title: 16s data processing and microbiome analysis
keywords: 16s rRNA, microbiome
tags: [genomics_analysis]
last_updated: September 17, 2014
summary: "This document briefly outlines the processing of 16s rRNA data, and provides guidance to key concepts and terms used"
sidebar: 16Smicro_sidebar
permalink: 16s-rRNA_5.html
folder: genomics_analysis/16s-rRNA/
author_profile: true
authors:
 - Gerrit_Botha
---



## Appendices {#appendices}

### **_Tools referred to in SOP_** {#tools}

*   FASTQC - [http://www.bioinformatics.babraham.ac.uk/projects/fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc)
*   PRINSEQ - [http://edwards.sdsu.edu/cgi-bin/prinseq/prinseq.cgi](http://edwards.sdsu.edu/cgi-bin/prinseq/prinseq.cgi)
*   SolexaQA - [http://www.biomedcentral.com/1471-2105/11/485](http://www.biomedcentral.com/1471-2105/11/485)
*   PEAR - [http://bioinformatics.oxfordjournals.org/content/early/2013/10/18/bioinformatics.btt593.full.pdf](http://bioinformatics.oxfordjournals.org/content/early/2013/10/18/bioinformatics.btt593.full.pdf)
*   PANDASeq - [http://www.biomedcentral.com/1471-2105/13/31](http://www.biomedcentral.com/1471-2105/13/31)
*   FLASH - [http://bioinformatics.oxfordjournals.org/content/early/2011/09/07/bioinformatics.btr507.full.pdf](http://bioinformatics.oxfordjournals.org/content/early/2011/09/07/bioinformatics.btr507.full.pdf)
*   UCHIME - [http://drive5.com/usearch/manual/uchime_algo.html](http://drive5.com/usearch/manual/uchime_algo.html)
*   ChimeraSlayer - [http://nebc.nox.ac.uk/bioinformatics/docs/chimeraslayer.html](http://nebc.nox.ac.uk/bioinformatics/docs/chimeraslayer.html)
*   Perseus - [http://www.biomedcentral.com/1471-2105/12/38/](http://www.biomedcentral.com/1471-2105/12/38/)
*   UPARSE - [http://www.drive5.com/uparse/](http://www.drive5.com/uparse/)
*   UCLUST - [http://www.drive5.com/uclust/downloads1_2_22q.html](http://www.drive5.com/uclust/downloads1_2_22q.html)
*   RDP classifier - [http://sourceforge.net/projects/rdp-classifier/files/rdp-classifier/](http://sourceforge.net/projects/rdp-classifier/files/rdp-classifier/)
*   RTAX - [https://github.com/davidsoergel/rtax](https://github.com/davidsoergel/rtax)
*   PyNAST - [http://www.ncbi.nlm.nih.gov/pubmed/19914921](http://www.ncbi.nlm.nih.gov/pubmed/19914921)
*   INFERNAL - [http://infernal.janelia.org/](http://infernal.janelia.org/)
*   FastTree - [http://www.microbesonline.org/fasttree/](http://www.microbesonline.org/fasttree/)
*   IM-TORNADO - [http://sourceforge.net/projects/imtornado/](http://sourceforge.net/projects/imtornado/)
*   Mothur - [http://www.mothur.org/](http://www.mothur.org/)
*   QIIME - [http://qiime.org/](http://qiime.org/)
*   UNIFRAC - [http://bmf.colorado.edu/unifrac/](http://bmf.colorado.edu/unifrac/)
*   R packages
    *   phyloseq - [http://www.bioconductor.org/packages/release/bioc/html/phyloseq.html](http://www.bioconductor.org/packages/release/bioc/html/phyloseq.html)
    *   ade4 - [http://cran.r-project.org/web/packages/ade4/index.html](http://cran.r-project.org/web/packages/ade4/index.html)

### **_Databases referred to in SOP_** {#databases}

*   SILVA - [http://www.arb-silva.de/](http://www.arb-silva.de/)
*   Greengenes - [http://greengenes.lbl.gov/](http://greengenes.lbl.gov/)
*   RDP classifier - [http://rdp.cme.msu.edu/](http://rdp.cme.msu.edu/)





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
