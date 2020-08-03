---
title: Node Accreditation Reviewer Guide
keywords: reviewer, report
tags: [misc]
last_updated: August 2, 2020
toc: true
folder: misc/
permalink: reviewers.html
author_profile: true
authors:
 - Chris_Fields
---

<!-- 
sidebar: 16Smicro_sidebar
 -->

## Introduction to node accreditation

Node accreditation exercises have been developed by the H3Africa Bioinformatics Network (H3ABioNet) to assess the ability of participating Nodes to perform common bioinformatics tasks. This is to ensure that these Nodes can undertake processing and analysis of data produced by the H3Africa consortium. Work towards this purpose was initiated during the first phase of H3ABioNet, with an initial assessment published in 2017 [(Jongeneel et al. 2017)](https://paperpile.com/c/zs0EoN/n7Xl).  This included the generation of:


1. standard operating procedures (SOPs) for four commonly performed tasks: variant calling, genome-wide association studies (GWAS), gene expression analysis using RNA-Seq, and 16S microbiome analysis;
2. data sets for each exercise, including a small set that can be used for initial test purposes and more realistic data used for the actual accreditation task;
3. a procedure outlining the process for undertaking the accreditation exercise, including the writeup of the report that summarizes the results and the review of the exercises.

Below we focus briefly on the SOPs, the procedure for taking the accreditation exercise, the report, and review expectations.


## Standard operating procedures

The SOPs are now hosted on Github:

[https://h3abionet.github.io/H3ABionet-SOPs/index.html](https://h3abionet.github.io/H3ABionet-SOPs/index.html)

The SOPs are not meant to be prescriptive but outline the overall procedure for analysis of each data type with some specific examples of useful tools where appropriate.  The SOPs are not completely focused on just procedure but also include thinking about the resources that may be needed for data processing and analysis, an important consideration as nodes may have limited access to high performance resources and thus should be capable of assessing these needs. It should be noted that H3ABioNet funded the purchase of capable servers for many Nodes, and that these should be shared among Nodes. 

As the documentation is now hosted on Github and are thus version controlled, the SOPs are regularly updated based on input from working groups in the Node Accreditation Task Force (NATF).  Each exercise also has some guiding questions that nodes who wish to undertake the exercise may want to consider when writing their report. 


## The accreditation exercise

Each exercise comprises three phases.



1. The candidate Node downloads the small sample datasets and the SOP for the method it wants to validate. Typically, a small team is formed that will take on the exercise together. During this first phase, the team puts in place and tests the workflows they are planning to use. There are predefined workflows available from the H3ABioNet consortium, or the team can build their own. Often, the team will get additional training at one of the H3ABioNet technical workshops. There is no time limit on this phase.
2. The candidate Node announces that it is ready for the actual exercise. It is provided with a unique link to download the relevant dataset for the exercise. It has six weeks from the day of download to perform the analysis,prepare a report and submit through a unique uplink path or an attachment to a provided email address. No extensions are granted except under truly exceptional circumstances.
3. Once the report is received, the HPCBio Node at the U. of Illinois gives it a first reading to ensure that it is complete and does not contain any gross errors; if the report does not follow these standards the exercise is considered as ‘not passing’.  If the exercise passes the initial review, the report is passed on to external reviewers for additional evaluation and comment. Our goal is to ensure that the external review process does not take more than four weeks, with the external reviewer’s report going a long way to determine the results of the exercise.


## The Node accreditation report

The report is a document that summarizes the overall (step-by-step) process and outcomes for the exercise.  This of course includes the results, but also should convey that the node understands the overall procedure, demonstrating that they are knowledgeable regarding the processing of the data and interpreting results. In the end they should be able to communicate how the data were processed (including code where appropriate), why they performed certain steps such as filtering and quality assessment of the sample data, demonstrate how they perform any downstream analyses, and then proceed to basic interpretation of results.  They also may include some benchmarking assessments on the tools used, as well as whether they utilized a particular workflow (and how it was implemented). 

However, it is important to keep in mind the nodes are largely left to decide for themselves how to prepare for and undertake the accreditation exercise.  As noted in the 2017 paper, they can choose different strategies and different methodologies.  Some nodes may have better overall infrastructure and more personnel and can thus dedicate more time and resources to training prior to the exercise, while others may have only one or two students with access to a shared server. We strongly recommend the candidate Nodes that more than one individual be involved, as this decreases the chances of mistakes due to gaps in individual training, and mitigates the loss of know-how upon departure / graduation of any one individual.


## Reviewer expectations

Our external reviewers are key to the success of the Node Accreditation program, as they alone can vouch for the fact that candidate Nodes have in fact attained international standards of excellence and can be trusted to analyze complex genomic datasets in a fully professional manner. Developing a pool of African scientists that can be trusted to produce such high-quality analyses is a key step in achieving scientific self-sufficiency and transitioning away from dependence on First World scientists and their “helicopter science”.

The reviewers are generally recruited via recommendations and word of mouth, with more active recruitment starting in spring 2019.  Minimally three reviewers are needed for each exercise, with at least one external review in addition to internal review performed by HPCBio.  We prefer to have two external reviewers when possible, though this largely depends on availability of reviewers during the period needed. 

Reviewers are requested to carefully read the report and address the following points:



1. Does the Node appear to understand the overall process for analysis of this dataset, as outlined in the SOPs?
2. Have they performed the analysis satisfactorily?  Have any critical errors occurred?
3. Have they carefully and thoroughly documented the steps they have taken, including any specific ad hoc code, software versions, inputs and outputs, and workflows?
4. Have they justified their choices regarding which tools to use, based on objective criteria?
5. Could the Node be considered experienced enough to analyze additional H3Africa datasets for this data type and to publish the results in a reputable Journal?
6. How could their overall analyses be improved?


## Reviewer summary report

**Time is of the essence! Please return your evaluation within four weeks of receiving the report, or let us know if you cannot meet that deadline.**

We do not provide a fixed template for the review, but request that it flags major and minor issues, makes suggestions on how to improve the analysis if needed, and provides an assessment of the readiness of the Node to participate in the analysis of real research data. The report should be evaluated on the same criteria as if it were produced by a group in your home country.

We rely on the opinions of two external reviewers to make a decision. This can be a pass (sometimes with congratulations for a job particularly well done), a conditional pass with an allowance of two weeks to revise the report based on reviewer comments, and a fail.


    [Jongeneel, C. Victor, Ovokeraye Achinike-Oduaran, Ezekiel Adebiyi, Marion Adebiyi, Seun Adeyemi, Bola Akanle, Shaun Aron, et al. 2017. “Assessing Computational Genomics Skills: Our Experience in the H3ABioNet African Bioinformatics Network.” PLoS Computational Biology 13 (6): e1005419.](http://paperpile.com/b/zs0EoN/n7Xl)

