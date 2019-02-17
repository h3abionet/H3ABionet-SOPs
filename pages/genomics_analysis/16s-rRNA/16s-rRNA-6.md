---
title: 16s data processing and microbiome analysis
keywords: 16s rRNA, microbiome
tags: [genomics_analysis]
last_updated: September 17, 2014
toc: false 
sidebar: 16Smicro_sidebar
permalink: 16s-rRNA-6.html
folder: genomics_analysis/16s-rRNA/
author_profile: true
authors:
 - Gerrit_Botha
---



### **_Tools referred to in SOP_** {#tools}

*   FASTQC - [http://www.bioinformatics.babraham.ac.uk/projects/fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc)
*   MultiQC - [https://multiqc.info/](https://multiqc.info/)
*   PRINSEQ - [http://edwards.sdsu.edu/cgi-bin/prinseq/prinseq.cgi](http://edwards.sdsu.edu/cgi-bin/prinseq/prinseq.cgi)
*   SolexaQA - [http://www.biomedcentral.com/1471-2105/11/485](http://www.biomedcentral.com/1471-2105/11/485)
*   PEAR - [http://bioinformatics.oxfordjournals.org/content/early/2013/10/18/bioinformatics.btt593.full.pdf](http://bioinformatics.oxfordjournals.org/content/early/2013/10/18/bioinformatics.btt593.full.pdf)
*   PANDASeq - [http://www.biomedcentral.com/1471-2105/13/31](http://www.biomedcentral.com/1471-2105/13/31)
*   FLASH - [http://bioinformatics.oxfordjournals.org/content/early/2011/09/07/bioinformatics.btr507.full.pdf](http://bioinformatics.oxfordjournals.org/content/early/2011/09/07/bioinformatics.btr507.full.pdf)
*   UCHIME - [http://drive5.com/usearch/manual/uchime_algo.html](http://drive5.com/usearch/manual/uchime_algo.html)
*   ChimeraSlayer - [http://nebc.nox.ac.uk/bioinformatics/docs/chimeraslayer.html](http://nebc.nox.ac.uk/bioinformatics/docs/chimeraslayer.html)
*   Perseus - [http://www.biomedcentral.com/1471-2105/12/38/](http://www.biomedcentral.com/1471-2105/12/38/)
*   UPARSE - [http://www.drive5.com/uparse/](http://www.drive5.com/uparse/)
*   vsearch = [https://github.com/torognes/vsearch](https://github.com/torognes/vsearch)
*   UCLUST - [http://www.drive5.com/uclust/downloads1_2_22q.html](http://www.drive5.com/uclust/downloads1_2_22q.html)
*   RDP classifier - [http://sourceforge.net/projects/rdp-classifier/files/rdp-classifier/](http://sourceforge.net/projects/rdp-classifier/files/rdp-classifier/)
*   RTAX - [https://github.com/davidsoergel/rtax](https://github.com/davidsoergel/rtax)
*   PyNAST - [http://www.ncbi.nlm.nih.gov/pubmed/19914921](http://www.ncbi.nlm.nih.gov/pubmed/19914921)
*   INFERNAL - [http://infernal.janelia.org/](http://infernal.janelia.org/)
*   FastTree - [http://www.microbesonline.org/fasttree/](http://www.microbesonline.org/fasttree/)
*   IM-TORNADO - [http://sourceforge.net/projects/imtornado/](http://sourceforge.net/projects/imtornado/)
*   Mothur - [http://www.mothur.org/](http://www.mothur.org/)
*   QIIME - [http://qiime.org/](http://qiime.org/)
*   QIIME2 - [https://qiime2.org/](https://qiime2.org/)
*   UNIFRAC - [http://bmf.colorado.edu/unifrac/](http://bmf.colorado.edu/unifrac/)
*   R packages
    *   phyloseq - [http://www.bioconductor.org/packages/release/bioc/html/phyloseq.html](http://www.bioconductor.org/packages/release/bioc/html/phyloseq.html)
    *   ade4 - [http://cran.r-project.org/web/packages/ade4/index.html](http://cran.r-project.org/web/packages/ade4/index.html)
    *   DADA2 - [https://benjjneb.github.io/dada2](https://benjjneb.github.io/dada2)

### **_Databases referred to in SOP_** {#databases}

*   SILVA - [http://www.arb-silva.de/](http://www.arb-silva.de/)
*   Greengenes - [http://greengenes.lbl.gov/](http://greengenes.lbl.gov/)
*   RDP classifier - [http://rdp.cme.msu.edu/](http://rdp.cme.msu.edu/)

## Bibliography {#bibliography}

[^1]: Bokulich, Nicholas A., et al. ["Quality-filtering vastly improves diversity estimates from Illumina amplicon sequencing."](https://www.nature.com/articles/nmeth.2276) Nature methods 10.1 (2013): 57.

[^2]:  Edgar, Robert C., et al. ["UCHIME improves sensitivity and speed of chimera detection."](https://academic.oup.com/bioinformatics/article/27/16/2194/255262) Bioinformatics 27.16 (2011): 2194-2200.

[^3]: Quince, Christopher, et al. ["Accurate determination of microbial diversity from 454 pyrosequencing data."](https://www.nature.com/articles/nmeth.1361) Nature methods 6.9 (2009): 639.

[^4]: Callahan, B. J., et al. ["DADA2: High-resolution sample inference from Illumina amplicon data."](https://www.nature.com/articles/nmeth.3869) Nature methods, 13.7 (2016), 581-3.

[^5]: Callahan, B. J., et al. ["Exact sequence variants should replace operational taxonomic units in marker-gene data analysis."](https://www.nature.com/articles/ismej2017119) The ISME journal, 11(12), 2639-2643.


[//]: <> (Below are the common abbreviations in the page.)
*[SOPs]: Standard Operating Procedures
*[16S rRNA]: 16S ribosomal RNA
*[16S rRNA gene]: The gene that is responsible for the coding of the 16S ribosomal RNA. The gene is used in constructing phylogenies.
*[Barcodes]: Short nucleotide sequences added onto the ends of the DNA fragments that are to be sequenced. It allows for indexing of samples, so multiple DNA libraries can be mixed together into one sequencing lane.
*[Variable region]: 16S rRNA gene sequences contain hypervariable regions that can provide clade-specific signature sequences useful for bacterial identification.
*[Demultiplex]: This is a process of binning reads based on barcodes, primarily used to split them amongst samples.
*[OTU]: An operational taxonomic unit is an operational definition of a species or group of species often used when only DNA sequence data is available.
*[Rarefaction analysis]: Rarefaction is a process used to estimate the true diversity of a sample by extracting random subsets of sequences. The analysis estimates diversity from subsets of different sizes and extrapolates the resulting rarefaction curve to an infinite number of sequences. Rarefaction is also used to determine whether the sequencing depth achieved (number of reads per sample) is sufficient to capture the diversity within a sample. A plot is generated showing increase in number of species (or other diversity metrics) as the number of sequence reads increase. A curve that is reaching asymptote indicates that no further diversity would be expected if sequencing depth was increased.
*[Phred quality score]: Measure of sequencing accuracy. Logarithmically related to the probability that a base is called incorrectly by a sequencer. Example: a Phred score of 30 (Q30) means that the probability of an incorrect call for that base is 1 in 1000 and the base call accuracy is 99.9%. The base call accuracy for a base with Q10 is 90%, Q20 is 99%, Q30 is 99.9%, Q40 is 99.99%, and Q50 = 99.999%.
*[Q score]: Measure of sequencing accuracy. Logarithmically related to the probability that a base is called incorrectly by a sequencer. Example: a Phred score of 30 (Q30) means that the probability of an incorrect call for that base is 1 in 1000 and the base call accuracy is 99.9%. The base call accuracy for a base with Q10 is 90%, Q20 is 99%, Q30 is 99.9%, Q40 is 99.99%, and Q50 = 99.999%.
*[Adapter]: Platform specific nucleotide sequence added to the ends of DNA molecule to facilitate sequencing e.g. in Illumina the adapter facilitates binding to the complementary target sequences mobilised on the flow cell.
*[Chimera]: PCR artefact. Chimeras are potentially formed during PCR when incompletely extended DNA fragments from different templates anneal due to closely related sequences generating recombinants between starting templates. Chimeras can greatly impact estimates of diversity (generally overestimate).
*[Alpha diversity]: Diversity within a single sample. Diversity can be characterised using the number of different species (richness), the abundance of each species (evenness), with indices that combine richness and evenness, and with divergence-based methods (phylogenetic diversity).
*[Beta diversity]: Comparison of diversity between samples.
*[Unifrac]: Beta diversity distance metric based on the phylogenetic distance between the members of communities/samples. Unifrac captures the total amount of evolution that is unique to each sample.
*[ASV]: Amplicon sequence variance

