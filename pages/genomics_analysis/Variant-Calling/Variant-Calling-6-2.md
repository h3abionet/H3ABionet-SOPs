---
title: Variant calling in human whole genome/exome sequencing data
keywords: wgs, wes, ngs
tags: [genomics_analysis]
last_updated: Fall, 2018

sidebar: varcall_sidebar
toc: false
permalink: Variant-Calling-5-2.html
folder: genomics_analysis/Variant-Calling
author_profile: true
authors:
 - Azza_Ahmed
 - Matthew_Weber
 - Faisal_Fadlelmola
 - Luidmila_Mainzer
---
### Alphabetized list of recommended tools {#alphabetized-list-of-recommended-tools}

 
<div markdown="span" class="alert alert-info" role="alert"><i class="fa fa-info-circle"></i> <b>Note on Galaxy:</b> 
If it is desirable to perform all processing in Galaxy [^52], then it is possible to construct a complete workflow by including the needed tools from its toolshed. The majority of the tools below can be found in the toolshed and incorporated readily. A complete coverage is beyond the scope of this SOP, and hence, the interested reader is referred to the galaxy project site for more details (https://galaxyproject.org/ )
</div>



*   
**ANNOVAR**--------------------------------------------------http://annovar.openbioinformatics.org/en/latest/
Utilizes update-to-date information to generate gene-based, region-based and filter-based functional annotations of genetic variants detected from diverse genomes.

_Usage</_



1.  Annotate the variants in ex1.human file and classify them as intergenic, intronic, non-synonymous SNP, frameshift deletion, large-scale duplication, etc. The ex1.human file is in text format, one variant per line. The annotation procedure takes seconds on a typical modern computer.

    ```
perl annotate_variation.pl -geneanno example/ex1.human humandb/
```


1.  Download cytogenetic band annotation databases from the UCSC Genome Browser and save it to the humandb/ directory as hg18_cytoBand.txt file, then annotate variants in ex1.human file and identify the cytogenetic band for these variants.

    ```
perl annotate_variation.pl -downdb cytoBand humandb/
perl annotate_variation.pl -regionanno -dbtype cytoBand example/ex1.human humandb/
```


1.  Download 1000 Genome Projects allele frequency annotations, then identify a subset of variants in ex1.human that are not observed in 1000G CEU populations and those that are observed with allele frequencies.

    ```
perl annotate_variation.pl -downdb 1000g humandb/
perl annotate_variation.pl -filter -dbtype 1000g_ceu example/ex1.human humandb/
```

<div markdown="span" class="alert alert-info" role="alert"><i class="fa fa-info-circle"></i> <b>Note:</b>  By default, all the above commands work on variants files in hg18 (human genome NCBI build 36) coordinate. If your file is in hg19 coordinate, add "-buildver hg19" in every command. </div>
 



*   
**Bed tools**-------------------------------------------------------------http://bedtools.readthedocs.io/en/latest/
A flexible suite of utilities for comparing genomic features.

_Filtering SNPs</_

To remove known SNPs, use intersectBed (bedtools). Known SNPs can be downloaded from Ensembl or NCBI/UCSC (eg. dbSNP)


```
intersectBed -a accepted_hits.raw.filtered.vcf \ 
             -b Mus_musculus.NCBIM37.60.bed \ 
             -wo \ 
             > filteredSNPs.vcf
```




*   
**Bioconductor**---------------------------------------------------------------------------[http://bioconductor.org/](http://bioconductor.org/)
Provides tools for the analysis and comprehension of high-throughput genomic data. Bioconductor uses the R statistical programming language, and is open source and open development.

_Sequence trimming</_

Use trimLRPatterns for adaptor trimming. From [http://manuals.bioinformatics.ucr.edu/home/ht-seq](http://manuals.bioinformatics.ucr.edu/home/ht-seq):


```
# Create sample sequences.
myseq <-DNAStringSet(c("CCCATGCAGACATAGTG", "CCCATGAACATAGATCC", "CCCGTACAGATCACGTG"));
names(myseq) <- letters[1:3];

# Remove the specified R/L-patterns. The number of maximum mismatches can be specified
# for each pattern individually. Indel matches can be specified with the arguments: 
# with.Lindels and with.Rindels.
trimLRPatterns(Lpattern ="CCC", Rpattern="AGTG", subject=myseq, max.Lmismatch = 0.33,
               max.Rmismatch = 1) 

# To remove partial adaptors, the number of mismatches for all possible partial matches
# can be specified by providing a numeric vector of length nchar(mypattern). 
# The numbers specifiy the number of mismatches for each partial match. 
# Negative numbers are used to prevent trimming of a minimum fragment length, 
# e.g. most terminal nucleotides.
trimLRPatterns(Lpattern = "CCC", Rpattern="AGTG", subject=myseq,
               max.Lmismatch=c(0,0,0), max.Rmismatch=c(1,0,0)) 

# With the setting 'ranges=TRUE' one can retrieve the corresponding 
# trimming coordinates.
trimLRPatterns(Lpattern = "CCC", Rpattern="AGTG", subject=myseq,
               max.Lmismatch=c(0,0,0), max.Rmismatch=c(1,0,0), ranges=T)
```


_Base quality score recalibration</_

Use ReQON to recalibrate quality of nucleotides for aligned sequencing data in BAM format: [http://bioconductor.org/packages/2.12/bioc/html/ReQON.html](http://bioconductor.org/packages/2.12/bioc/html/ReQON.html) See ReQON tutorial for detailed examples of usage: [http://bioconductor.org/packages/devel/bioc/vignettes/ReQON/inst//doc/ReQON.pdf](http://bioconductor.org/packages/devel/bioc/vignettes/ReQON/inst//doc/ReQON.pdf) 

 



*   
**blat**----------------------------------------------------------------[http://genome.ucsc.edu/FAQ/FAQblat.html](http://genome.ucsc.edu/FAQ/FAQblat.html)
       ------------------------------------------------[http://genome.ucsc.edu/goldenPath/help/blatSpec.html](http://genome.ucsc.edu/goldenPath/help/blatSpec.html)

Blat is an alignment tool like BLAST, but it is structured differently. On DNA, Blat works by keeping an index of an entire genome in memory. From a practical standpoint, Blat has several advantages over BLAST:



1.  speed (no queues, response in seconds) at the price of lesser homology depth
1.  the ability to submit a long list of simultaneous queries in fasta format
1.  five convenient output sort options
1.  a direct link into the UCSC browser
1.  alignment block details in natural genomic order
1.  an option to launch the alignment later as part of a custom track

_Usage</_


```
blat <database> <query> [-ooc=11.ooc] output.psl
```


Database and query are each either a `.fa` , `.nib` or `.2bit` file, or a list these files.

The option `-ooc=11.ooc` tells the program to load over-occurring 11-mers from an external file.  This will increase the speed by a factor of 40 in many cases, but is not required.

 



*   
**Bowtie2**---------------------------------------------[http://bowtie-bio.sourceforge.net/bowtie2/index.shtml](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
An ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. It is particularly good at aligning reads of about 50 up to 100s or 1,000s of characters to relatively long (e.g. mammalian) genomes.

_Usage</_


```
bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r>} -S [<sam>]

-x <bt2-idx>  --  The basename of the index for the reference genome. The basename is the name of any of the index files up to but not including the final .1.bt2 / .rev.1.bt2 / etc. bowtie2 looks for the specified index first in the current directory, then in the directory specified in the BOWTIE2_INDEXES environment variable.

-1 <m1>  --  Comma-separated list of files containing mate 1s (filename usually includes _1), e.g. -1 flyA_1.fq,flyB_1.fq. Sequences specified with this option must correspond file-for-file and read-for-read with those specified in <m2>. Reads may be a mix of different lengths. If - is specified, bowtie2 will read the mate 1s from the "standard in" or "stdin" filehandle.

-2 <m2>  --  Comma-separated list of files containing mate 2s (filename usually includes _2), e.g. -2 flyA_2.fq,flyB_2.fq. Sequences specified with this option must correspond file-for-file and read-for-read with those specified in <m1>. Reads may be a mix of different lengths. If - is specified, bowtie2 will read the mate 2s from the "standard in" or "stdin" filehandle.

-U <r>   --  Comma-separated list of files containing unpaired reads to be aligned, e.g. lane1.fq,lane2.fq,lane3.fq,lane4.fq. Reads may be a mix of different lengths. If - is specified, bowtie2 gets the reads from the "standard in" or "stdin" filehandle.

-S <sam> --  File to write SAM alignments to. By default, alignments are written to the "standard out" or "stdout" filehandle (i.e. the console).
```


 



*   
**BWA**------------------------------------------------------------------------------[http://bio-bwa.sourceforge.net/](http://bio-bwa.sourceforge.net/)
Burrows-Wheeler Alignment Tool. A software package for mapping low-divergent sequences against a large reference genome, such as the human genome. It consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence reads up to 100bp, while the rest two for longer sequences ranged from 70bp to 1Mbp. BWA-MEM and BWA-SW share similar features such as long-read support and split alignment, but BWA-MEM, which is the latest, is generally recommended for high-quality queries as it is faster and more accurate. BWA-MEM also has better performance than BWA-backtrack for 70-100bp Illumina reads.

_Usage</_


```
Index database sequences:
bwa index ref.fa
 
Align 70bp-1Mbp query sequences with the BWA-MEM algorithm:
bwa mem ref.fa reads.fq > aln-se.sam #Single end reads
bwa mem ref.fa read1.fq read2.fq > aln-se.sam #For functional equivalence 9, add the options: -K 100000000 -Y
 
Find the SA coordinates of the input reads.
bwa aln ref.fa short_read.fq > aln_sa.sai
 
Generate alignments in the SAM format given single-end reads shorter than ~ 70bp

 bwa aln ref.fa short_reads.fq > reads.sai 
 bwa samse ref.fa reads.sai short_reads.fq > aln-se.sam

 
Generate alignments in the SAM format given paired-end reads shorter than ~70bp.
bwa aln ref.fa read1.fq > read1.sai; bwa aln ref.fa read2.fq > read2.sai
bwa sampe ref.fa read1.sai read2.sai read1.fq read2.fq > aln-pe.sam
```


 



*   
**FASTX-Toolkit**-------------------------------------------------------[http://hannonlab.cshl.edu/fastx_toolkit/](http://hannonlab.cshl.edu/fastx_toolkit/)
A collection of command line tools for Short-Reads FASTA/FASTQ files preprocessing.

_fastq_quality_filter</_

Filters sequences based on quality.  [http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastq_quality_filter_usage](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastq_quality_filter_usage) 


```
fastq_quality_filter [-q N] [-p N] [-z] [-i INFILE] [-o OUTFILE]

[-q N]  =  Minimum quality score to keep.
[-p N]  =  Minimum percent of bases that must have [-q] quality.
[-z]    =  Compress output with GZIP.
[-i INFILE]  =  FASTA/Q input file. default is STDIN.
[-o OUTFILE]  =  FASTA/Q output file. default is STDOUT.
```


_fastx-collapser</_

Collapses identical sequences in a FASTQ/A file into a single sequence). [http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastx_collapser_usage](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastx_collapser_usage) 


```
fastx_collapser [-i INFILE] [-o OUTFILE]
```




*   
**Flexbar**----------------------------------------------------------[ https://github.com/seqan/flexbar](http://sourceforge.net/p/flexbar/wiki/Manual/)
Flexbar — flexible barcode and adapter removal. Flexbar is a software to preprocess high-throughput sequencing data efficiently. It demultiplexes barcoded runs and removes adapter sequences. Trimming and filtering features are provided. Flexbar increases mapping rates and improves genome and transcriptome assemblies. It supports next-generation sequencing data from Illumina, Roche 454, and the SOLiD platform. Recognition is based on exact overlap sequence alignment.

_Usage</_


```
flexbar -t target -r reads [-b barcodes] [-a adapters] [options]

Mandatory parameters:
input file with reads
input format of reads
target name for output prefix
```


 



*   
**GATK**------------------------------------------------------------------------[http://www.broadinstitute.org/gatk/](http://www.broadinstitute.org/gatk/)
Genome Analysis Toolkit - a software package developed at the Broad Institute to analyse next-generation resequencing data. The toolkit offers a wide variety of tools, with a primary focus on variant discovery and genotyping as well as strong emphasis on data quality assurance. Its robust architecture, powerful processing engine and high-performance computing features make it capable of taking on projects of any size.

There is a slight change when calling tools from either GATK version:


```
GATK3 invocation: 
java -jar <jvm args like -Xmx4G go here> GenomeAnalysisTK.jar -T <ToolName>

GATK4 invocation:
gatk [--java-options <jvm args like -Xmx4G go here>] <ToolName> [GATK args go here]
	
               
```


_BQSR</_             [http://gatkforums.broadinstitute.org/discussion/44/base-quality-score-recalibration-bqsr](http://gatkforums.broadinstitute.org/discussion/44/base-quality-score-recalibration-bqsr)

The tools in this package recalibrate base quality scores of sequencing-by-synthesis reads in an aligned BAM file. After recalibration, the quality scores in the QUAL field in each read in the output BAM are more accurate in that the reported quality score is closer to its actual probability of mismatching the reference genome. Moreover, the recalibration tool attempts to correct for variation in quality with machine cycle and sequence context, and by doing so provides not only more accurate quality scores but also more widely dispersed ones. The system works on BAM files coming from many sequencing platforms: Illumina, SOLiD, 454, Complete Genomics, Pacific Biosciences, etc. Invocation based on recommended functional equivalent parameters [^9] is below: 


```
java -Xmx4g -jar GenomeAnalysisTK.jar \
               -T BaseRecalibrator
               -R ${ref_fasta} \
               -I ${input_bam} \
               -O ${recalibration_report_filename} \
               -knownSites "Homo_sapiens_assembly38.dbsnp138.vcf" \
               -knownSites "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
               -knownSites "Homo_sapiens_assembly38.known_indels.vcf.gz"


```


To create a recalibrated BAM you can use GATK's PrintReads (GATK3), or ApplyBQSR (GATK4) with the engine on-the-fly recalibration capability. Here is the recommended command line to do so while complying with the functional equivalence specification [^9]:  


```
java -jar GenomeAnalysisTK.jar \	
               -T PrintReads \
               -R reference.fasta \
               -I input.bam \    
               -BQSR recalibration_report.grp \    
               -o output.bam
               -SQQ 10 -SQQ 20 -SQQ 30 \
               --disable_indel_quals

```


_IndelRealigner</_

Performs local realignment of reads to correct misalignments due to the presence of indels. Realginemet is no longer needed with a haplotype-aware variant caller like the HaplotypeCaller  [^39] , but this is provided for completion [https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_indels_IndelRealigner.php) 


```
java -Xmx4g -jar GenomeAnalysisTK.jar \
               -T IndelRealigner \
               -R ref.fasta \
               -I input.bam \
               -targetIntervals intervalListFromRTC.intervals \
               -o realignedBam.bam \
              [-known /path/to/indels.vcf]
```


_HaplotypeCaller</_

Call SNPs and indels simultaneously via local de-novo assembly of haplotypes in an active region. [https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) 


```
java -jar GenomeAnalysisTK.jar
 	        -T HaplotypeCaller \
               -R reference/human_g1k_v37.fasta \
               -I sample1.bam [-I sample2.bam ...] \
               --dbsnp dbSNP.vcf \
               -stand_call_conf [50.0] \
               -stand_emit_conf 10.0 \
              [-L targets.interval_list] \
               -o output.raw.snps.indels.vcf
```


_MuTect2</_

Call SNPs and indels simultaneously via local de-novo assembly of haplotypes, combining the MuTect genotyping engine and the assembly-based machinery of HaplotypeCaller. For cancer variant discovery. Tumor/Normal or Normal-only calls.

[https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_cancer_m2_MuTect2.php)


```
java -jar GenomeAnalysisTK.jar
               -T MuTect2 \
               -R reference.fa \
               -I:tumor tumor.bam \
               -I:normal normal.bam \
               [--dbsnp sbSNP.vcf] \
               [--cosmic COSMIC.vcf] \
               [-L targets.interval_list] \
               -o output.vcf
```


_UnifiedGenotyper</_

A variant caller which unifies the approaches of several disparate callers -- Works for single-sample and multi-sample data. In most applications, the recommendation is to use the newer HaplotypeCaller [https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.php](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.php) 


```
java -jar GenomeAnalysisTK.jar
               -R resources/Homo_sapiens_assembly18.fasta \
               -T UnifiedGenotyper \
               -I sample1.bam [-I sample2.bam ...] \
               --dbsnp dbSNP.vcf \
               -o snps.raw.vcf \
               -stand_call_conf [50.0] \
               -stand_emit_conf 10.0 \
               -dcov [50 for 4x, 200 for >30x WGS or Whole exome] \
              [-L targets.interval_list]
```


_GATK Variant Quality Score Recalibrator</_

Creates a Gaussian mixture model by looking at the annotations values over a high quality subset of the input call set and then evaluate all input variants. [http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html](http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html)


```
 java -Xmx4g -jar GenomeAnalysisTK.jar
   -T VariantRecalibrator \
   -R reference/human_g1k_v37.fasta \
   -input NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.b37.vcf \
   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 \ 
                                                     hapmap_3.3.b37.sites.vcf \
   -resource:omni,known=false,training=true,truth=false,prior=12.0 \
                                                     1000G_omni2.5.b37.sites.vcf \
   -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 dbsnp_135.b37.vcf \
   -an QD -an HaplotypeScore -an MQRankSum \
   -an ReadPosRankSum -an FS -an MQ -an InbreedingCoeff \
   -mode SNP \
   -recalFile path/to/output.recal \
   -tranchesFile path/to/output.tranches \
   -rscriptFile path/to/output.plots.R
```




*   
**HuVariome**--------------------------------------------------------------------[http://huvariome.erasmusmc.nl/](http://huvariome.erasmusmc.nl/)
The HuVariome project aims to determine rare and common genetic variation in a Northern European population (Benelux) based on whole genome sequencing results. Variations, their population frequencies and the functional impact are stored in the HuVariome Database. Users who provide samples for inclusion within the database are able to access the full content of the database, whilst guests can access the public set of genomes published by Complete Genomics.

Has a straightforward web interface.

 



*   
**KGGSeq**---------------------------------------------------[http://statgenpro.psychiatry.hku.hk/limx/kggseq/](http://statgenpro.psychiatry.hku.hk/limx/kggseq/)
A biological Knowledge-based mining platform for Genomic and Genetic studies using Sequence data.

_Usage</_

To prioritize variants based on the hg18 assembly for a rare Mendelian disease, need input files:



1.  a VCF file (rare.disease.hg18.vcf); and
1.  a linkage pedigree file (rare.disease.ped.txt).

    ```
java  -jar   -Xms256m   -Xmx1300m  kggseq.jar   examples/param.rare.disease.hg18.txt
```





*   
**NEAT-genReads**---------------------------------------[https://github.com/zstephens/neat-genreads.git](https://github.com/zstephens/neat-genreads.git)
NEAT-genReads is a fine-grained read simulator published in PLoS One  [^53]. GenReads simulates real-looking data using models learned from specific datasets. Simulated reads can be whole genome, whole exome, specific targeted regions, or tumor/normal, with optional vcf and bam file outputs. Additionally, there are several supporting utilities for generating models used for simulation. Requires Python 2.7 and Numpy 1.9.1+. 

_Usage:</_ 


```
Simulate whole genome dataset with random variants inserted according to the default model. Generates paired-end reads, bam, and vcf output.

python genReads.py                  \
        -r hg19.fa                  \
        -R 101                      \
        -o /home/me/simulated_reads \
        --bam                       \
        --vcf                       \
        --pe 300 30
```



```
Simulate a targeted region of a genome, e.g. exome, with -t. Generates paired-end reads, bam, and vcf output.

python genReads.py                  \
        -r hg19.fa                  \
        -R 101                      \
        -o /home/me/simulated_reads \
        --bam                       \
        --vcf                       \
        --pe 300 30                 \
        -t hg19_exome.bed
```




*   
**Picard MarkDuplicates**----- [https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.4.0/picard_sam_markduplicates_MarkDuplicates.php](https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.4.0/picard_sam_markduplicates_MarkDuplicates.php) 
Examines aligned records in the supplied SAM or BAM file to locate duplicate molecules. All records are then written to the output file with the duplicate records flagged.

Usage</

Starting with GATK4.0, all Picard tools are directly available from the GATK suit itself. This means that for older GATK releases, one would call a tool as below:


```
java -jar picard.jar MarkDuplicates [options] INPUT=File OUTPUT=File
```


Whereas for GATK4 releases, the same tool would be called as below:


```
gatk MarkDuplicates [options] --INPUT File --OUTPUT File
```


For a comprehensive list of options, see the utility's documentation page listed above.

 



*   
**Multiple databases can be used to annotate variants**--------------[http://regulome.stanford.edu/](http://regulome.stanford.edu/)
A database that annotates SNPs with known and predicted regulatory elements in the intergenic regions of the H. sapiens genome. Known and predicted regulatory DNA elements include regions of DNAase hypersensitivity, binding sites of transcription factors, and promoter regions that have been biochemically characterized to regulation transcription. Source of these data include public datasets from GEO, Chromatin States from the Roadmap Epigenome Consortium, the ENCODE project, and published literature.

Has a straightforward web interface.



*   
**SAMtools**----------------------------------------------------------------------[http://samtools.sourceforge.net/](http://samtools.sourceforge.net/)
                ----------------------------------------------------------------------[http://www.htslib.org](http://www.htslib.org) 

SAM Tools provide various utilities for manipulating alignments in the SAM format, including sorting, merging, indexing and generating alignments in a per-position format.

_SNP calling with mpileup</_                                               [http://samtools.sourceforge.net/mpileup.shtml](http://samtools.sourceforge.net/mpileup.shtml)


```
# Use Q to set base quality threshold:
samtools mpileup -Q 25 -ugf <reference.fasta> <file.bam> | \
         bcftools view -bvcg - > accepted_hits.raw.bcf

# Convert bcf to vcf (vairant call format) 
# and filter using varFilter (from vcfutils), if needed; 
# use Q to set mapping quality; use d for minimum read depth:
bcftools view accepted_hits.raw.bcf | vcfutils.pl varFilter -d 5 -Q 20 > \
         accepted_hits.raw.vcf
```




*   
**SolexaQA++**----------------------------------------------------------------------[http://solexaqa.sourceforge.net/](http://solexaqa.sourceforge.net/)
SolexaQA is a software package to calculate sequence quality statistics and create visual representations of data quality for Illumina's second-generation sequencing technology (historically known as "Solexa").

_Usage</_

Running directly on Illumina FASTQ files:


```
./SolexaQA++ analysis FASTQ_input_files \
     [-p|probcutoff 0.05] \
     [-h|phredcutoff] \
     [-v|variance] \
     [-m|minmax] \
     [-s|sample 10000] \
     [-b|bwa] \
     [-d|directory path] \
     [-sanger -illumina -solexa]
```


 



*   
**snpEff**------------------------------------------------------------------------------[http://snpeff.sourceforge.net/](http://snpeff.sourceforge.net/)
Genetic variant annotation and effect prediction toolbox. Version 2.0.5 is is included in the GATK suit, but not newer versions. For human data, the recommendation is to use the **GRCh37.64 database**

_Usage</_


```
# If you don't already have the database installed:
java -jar snpEff.jar download -v GRCh37.66

# Annotate the file:
# Use appropriate i) organism  and ii) annotation (eg. UCSC, RefSeq, Ensembl)
# E.g. mouse: mm37 (UCSC/RefSeq), mm37.61 (Ensembl)
# Human: hg37 (UCSC/RefSeq), hg37.61 (Ensembl)
# Output is created in several files: 
# an html summary file and text files with detailed information.
java -jar snpEff.jar -vcf4 GRCh37.75  accepted_hits.raw.filtered.vcf >accepted_hits.raw.filtered.snpEff
```




*   
**SnpSift**-------------------------------------------------------------[http://snpeff.sourceforge.net/SnpSift.html](http://snpeff.sourceforge.net/SnpSift.html)
Helps filter and manipulate annotated files. Included in GATK.

Once your genomic variants have been annotated, you need to filter them out in order to find the "interesting / relevant variants". Given the large data files, this is not a trivial task (e.g. you cannot load all the variants into XLS spreasheet). SnpSift helps to perform this VCF file manipulation and filtering required at this stage in data processing pipelines.

_Usage</_

To filter out samples with quality less than 30:


```
cat variants.vcf   |   java -jar SnpSift.jar filter " ( QUAL >= 30 )"   >  filtered.vcf
```


To do the same but keep InDels that have quality 20 or more:


```
cat variants.vcf   |   java -jar SnpSift.jar filter \
                       "(( exists INDEL ) & (QUAL >= 20)) | (QUAL >= 30 )" \
                       >   filtered.vcf
```




*   

The Variant Annotation, Analysis and Search Tool) is a probabilistic search tool for identifying damaged genes and their disease-causing variants in personal genome sequences. VAAST builds upon existing amino acid substitution (AAS) and aggregative approaches to variant prioritization, combining elements of both into a single unified likelihood-framework that allows users to identify damaged genes and deleterious variants with greater accuracy, and in an easy-to-use fashion. VAAST can score both coding and non-coding variants, evaluating the cumulative impact of both types of variants simultaneously. VAAST can identify rare variants causing rare genetic diseases, and it can also use both rare and common variants to identify genes responsible for common diseases. VAAST thus has a much greater scope of use than any other existing methodology. 

It has outstanding quickstart quide:  [http://www.yandell-lab.org/software/VAAST/VAAST_Quick-Start-Guide.pdf](http://www.yandell-lab.org/software/VAAST/VAAST_Quick-Start-Guide.pdf) 

_Usage</_

The set of genomes being analyzed for disease causing features (cases) are referred to as the target genomes. The set of healthy genomes (controls) that the target genomes are being compared to are referred to as the background genomes. The basic inputs to VAAST consist of:



1.  A set of target (case) variant files in either VCF or GVF format.
1.  A set of background (control) variant files in either VCF or GVF format.
1.  A set features to be scored, usually genes; this file should be in gff3 format.
1.  A multi-fasta file of the reference genome

    ```
 ../bin/VAT -f data/easy-hg18-chr16.gff3 \
            -a data/hg18-chr16.fasta data/miller-1.gvf > data/miller-1.vat.gvf
```



VCF can be easily converted to GVF using the vaast_converter script, included with the distribution: VAAST/bin/vaast_tools/vaast_converter.



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
