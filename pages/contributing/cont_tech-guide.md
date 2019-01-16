---
title: Technical Guide 
keywords: contribution, guide, tips 
tags: [documentation, formatting]
last_updated: 1 2, 2019
summary: "This page describes a few technical tips for authoring content on this site" 
sidebar: cont_sidebar
permalink: cont_tech-guide.html  
folder: contributing
author_profile: true
authors:
 - Azza_Ahmed 
---

## How to contribute {#overview}

Thank you for your interest in contributing to this effort! These SOPs are published openly with the aim of helping the greater bioinformatics community, so we love to recieve contributions. There are many ways to do so, from reporting typos, issues, suggesting content, ... etc.

 To make it easier, we created this guide- let's dive in!

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
│   └── images
│       ├── 16s-overview.png
│       ├── author1_image.png
│       ├── author2_image.jpeg
│       └── VarCall.png
├── pdf
│   └── mydoc.pdf
├── pdfconfigs
│   ├── config_mydoc_pdf.yml
│   ├── prince-list.txt
│   ├── titlepage.html
│   └── tocpage.html
└── pdf-mydoc.sh 
```
## Road map {#road-map}

- Split lenghty pages into smaller chunks
- Create diagrams for SOP pages
- Add authors' contact details next to the pages they contribute to
- Generate pdf files of the SOPs

## Tips for contributing new content: {#handy-tips}


<div markdown="span" class="alert alert-success" role="alert"><i class="fa fa-check-square-o"></i> <b>Tip: </b> This site is built using [Jekyll](https://jekyllrb.com/) and hosted on [github pages](https://pages.github.com/). A complete and detailed coverage of the Jekyll theme used for building this site, the Documentation theme, is available [here](https://idratherbewriting.com/documentation-theme-jekyll/). </div>

- Format your content in markdown (the theme uses the `kramdown` flavour), and place your markdown file within the `pages` folder. If you have been working with Google docs, then the [GD2md-html Google Docs add-on](https://github.com/evbacher/gd2md-html/wiki) can do the conversion efficiently. If your content is lenghthy, then you may wish to split it up for more managable navigation.

- Append appropriate frontmatter to the top of your markdown. You can simply copy-paste from the file [`template_sop.md`](https://github.com/h3abionet/H3ABionet-SOPs/blob/master/template_sop.md) as a starter.

- For the `tags` to work, you still need to define them in `_data/tags.yml` file. Additionally, you need to create corresponding markdown pages in `pages/tags` with the tag name and 

- For the sidebar, it should refer to a file of the intended table of contents found in the directory `_data/sidebars`.  An example file is provided [sops_sidebar.yml](https://github.com/h3abionet/H3ABionet-SOPs/blob/master/_data/sidebars/sops_sidebar.yml). Simply, put the name of this file in your frontmatter above for it to work. For more information about its format, refer to [this guide](https://idratherbewriting.com/documentation-theme-jekyll/#configure-the-sidebar) and [this one](https://idratherbewriting.com/documentation-theme-jekyll/#sidebar-syntax).


- For top navigation, simply edit the file: `_data/topnav.yml`

- It is desirable to have a `hover-over` tooltip helpers next to key abbreviations or terms in each SOP. The convention here is to put such content at the end of a given SOP file, and start the section with `[\\]:`. Next, put each term in tits own line, `*[term]: explanation`. Again, the [template page](https://github.com/h3abionet/H3ABionet-SOPs/blob/master/template_sop.md) has some examples.

- `References` are included as footnotes in this site. In markdown, this means using a `[^label]` for inline citations, and `[^label]: complete citation details in Vancuver style?` for biblography. Note that the `label` does not need to be numeric, as this will be done automatically by the engine

<div markdown="span" class="alert alert-success" role="alert"><i class="fa fa-check-square-o"></i> <b>Tip: </b> I find it easier to numerically number referneces when migrating a lenghty page (say from a Google word document) into the site. However, once all is settled, and future edits are made, labelling referneces (in a `latex` manner for example) is easier to do (and maintain)! </div>

<div markdown="span" class="alert alert-danger" role="alert"><i class="fa fa-exclamation-circle"></i> <b>Warning: </b> Footnotes do not appear within tables. </div>

- To insert special boxes to your page: warning, note, tip, info, ... etc, you may use the following `div`s in your markdown:

|**Type**         | **class** |
| --------------- | -----------
| tip             | `<div markdown="span" class="alert alert-success" role="alert"><i class="fa fa-check-square-o"></i> <b>Tip: </b> CONTENT </div>`
| note            | `<div markdown="span" class="alert alert-info" role="alert"><i class="fa fa-info-circle"></i> <b>Note: </b> CONTENT </div>`
| important       | `<div markdown="span" class="alert alert-warning" role="alert"><i class="fa fa-warning"></i> <b>Important: </b> CONTENT </div>`
| warning         | `<div markdown="span" class="alert alert-danger" role="alert"><i class="fa fa-exclamation-circle"></i> <b>Warning: </b> CONTENT </div>`
| callout_danger  | `<div markdown="span" class="bs-callout bs-callout-danger"> CONTENT </div>`
| callout_default | `<div markdown="span" class="bs-callout bs-callout-default"> CONTENT </div>`
| callout_primary | `<div markdown="span" class="bs-callout bs-callout-primary"> CONTENT </div>`
| callout_success | `<div markdown="span" class="bs-callout bs-callout-success"> CONTENT </div>`
| callout_info    | `<div markdown="span" class="bs-callout bs-callout-info"> CONTENT </div>`
| callout_warning | `<div markdown="span" class="bs-callout bs-callout-warning">`
| hr_faded        | `<hr class="faded"/>`
| hr_shaded       | `<hr class="shaded"/>`


## Changelog {#changelog}
