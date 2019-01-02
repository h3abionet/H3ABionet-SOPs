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

- Add authors' contact details next to the pages they contribute to
- Generate pdf files of the SOPs

## Tips for contributing new content: {#handy-tips}


<div markdown="span" class="alert alert-success" role="alert"><i class="fa fa-check-square-o"></i> <b>Tip: </b> This site is built using [Jekyll](https://jekyllrb.com/) and hosted on [github pages](https://pages.github.com/). A complete and detailed coverage of the Jekyll theme used for building this site, the Documentation theme, is available [here](https://idratherbewriting.com/documentation-theme-jekyll/). </div>

- Format your content in markdown (the theme uses the `kramdown` flavour), and place your markdown file within the `pages` folder. If you have been working with Google docs, then the [GD2md-html Google Docs add-on](https://github.com/evbacher/gd2md-html/wiki) can do the conversion efficiently. If your content is lenghthy, then you may wish to split it up for more managable navigation.

- Append appropriate frontmatter to the top of your markdown. Below is a generic example, that you may copy-paste from the file [`template_frontmatter.md`](https://github.com/h3abionet/H3ABionet-SOPs/blob/master/template_frontmatter.md) as well.


```

---
title: <title of your page- will appear at the browser window>
tags: [comma seperated list of tags] 
keywords: <relevant keywords to your content, for search engine> 
last_updated: <date of last update- appears at the bottom of the pate>
summary: <summary of the page if desired> 
sidebar: <The *.yml navigation file in _data/sidebars (see below for details)> 
hide_sidebar: true #Use this if you don't wish to include a sidebar. Otherwise, the default sops_sidebar.yml will appear
permalink: <the same name as your file, with .html extension>
folder: <the name of the folder within the pages directory that contains your markdown content>
author_profile: true # if you wish to add authors details
authors:
 - <list of authors, whose details are in the authors.yml file>
---

<Your intended page content, formatted in kramdown. >

```

- For the `tags`, you still need to define them in `_data/tags.yml` file. Additionally, you need to create corresponding markdown pages in `pages/tags` with the tag name and 

- For the sidebar, it should refer to a file of the intended table of contents to be found in _data/sidebars  An example is sops_sidebar.yml, and it can be added as below. For more information about its format, refer to [this guide](https://idratherbewriting.com/documentation-theme-jekyll/#configure-the-sidebar) and [this one](https://idratherbewriting.com/documentation-theme-jekyll/#sidebar-syntax).


```
sidebar: sops_sidebar
```

- For top navigation, simply edit the file: `_data/topnav.yml`

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
