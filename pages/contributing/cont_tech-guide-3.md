---
title: Technical Guide 
keywords: contribution, guide, tips 
tags: [documentation, formatting]
last_updated: 1 2, 2019
toc: false
sidebar: cont_sidebar
permalink: cont_tech-guide-3.html  
folder: contributing
author_profile: true
authors:
 - Azza_Ahmed 
---

## Project tree {#tree} 

Below is a tree of key files and folders within this site that would need editing for adding new content. This is a simplistic tree to suggest where to place content, and is not meant to help with functional changes to the theme itself.

``` bash
.
├── _config.yml          #Specific site configurations
├── index.md             #Main landing page of the site
├── authors.yml          #Contact details for authors- inactive now
├── README.md            #Github readme file
├── template_frontmatter.md   #Sample page frontmatter file
├── pages                     #Folder of site pages, arranged per the top navigation bar categories
│   ├── genomics_analysis
│   │   ├── 16s-rRNA
│   │   │   └── 16s-rRNA.md
│   │   └── Variant-Calling
│   │       └── Variant-Calling.md
│   └── tags                 #Folder of all tags in all pages
│       ├── tag_documentation.md
│       └── tag_genomics_analysis.md
├── _data
│   ├── sidebars            #Folder of all sidebars in the site
│   │   ├── sops_sidebar.yml
│   │   └── varcall_sidebar.yml
│   ├── tags.yml            #File of all tags present in the site
│   └── topnav.yml          #File of top navigation contents
├── assets
│   └── images              #Folder of images used within the site
│       └── VarCall.png
├── scripts                 #Various scripts to help with automation of tasks 
│   └── split_by_h2.sh     
├── _includes
│   └── custom             
│       ├── usermap.html   #Template file for simple workflow maps
│       └── usermapcomplex.html  #Template file for complex maps
├── pdf
│   └── mydoc.pdf
├── pdfconfigs
│   ├── config_mydoc_pdf.yml
│   ├── prince-list.txt
│   ├── titlepage.html
│   └── tocpage.html
└── pdf-mydoc.sh 

```

