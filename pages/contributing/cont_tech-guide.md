---
title: Technical Guide 
keywords: contribution, guide, tips 
tags: [documentation, formatting]
last_updated: 1 2, 2019
summary: "This page describes a few technical tips for authoring content on this site" 
sidebar: cont_sidebar
permalink: cont_tech-guide.html  
author_profile: true
authors:
 - Azza_Ahmed 
---

## How to contribute {#overview}

Thank you for your interest in contributing to this effort! These SOPs are published openly with the aim of helping the greater bioinformatics community, so we love to recieve contributions. There are many ways to do so, from reporting typos, issues, suggesting content, ... etc.

 To make it easier, we created this guide- let's dive in!

## Project tree {#tree} 


## Road map {#road-map}
[] Add authors' contact details next to the pages they contribute to


## Tips for contributing new content: {#handy-tips}


<div markdown="span" class="alert alert-success" role="alert"><i class="fa fa-check-square-o"></i> <b>Tip: </b> This site is built using [Jekyll](https://jekyllrb.com/) and hosted on [github pages](https://pages.github.com/). A complete and detailed coverage of the Jekyll theme used for building this site, the Documentation theme, is available [here](https://idratherbewriting.com/documentation-theme-jekyll/). </div>

- To add ad-hoc pages, just place the markdown file in the root directory. Add a frontmatter similar to this at top (or copy-paste and modify the contents of the file `template_frontmatter.md`). Your content should follow.

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
author_profile: true # if you wish to add authors details
authors:
 - <list of authors, whose names are in the authors.yml file>
---

<Your intended page content, formatted in markdown. >

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
