---
title: Technical Guide 
keywords: contribution, guide, tips 
tags: [documentation, formatting]
last_updated: 1 2, 2019
toc: false
sidebar: cont_sidebar
permalink: cont_tech-guide-2.html  
folder: contributing
author_profile: true
authors:
 - Azza_Ahmed 
---

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

- For a better user experience, it is desirable to limit the lenght of a given SOP page. The helper scripts: [scripts/split_by_h2.sh](https://github.com/h3abionet/H3ABionet-SOPs/blob/master/scripts/split_by_h2.sh) and [scripts/split_by_h3.sh](https://github.com/h3abionet/H3ABionet-SOPs/blob/master/scripts/split_by_h3.sh) split up a properly formatted markdown page by level 2 and level 3 respectively. <a href="Variant-Calling.html">This</a> lengthy page was split up into its level2 and 3 components <a href="Variant-Calling-1-0.html">in those pages</a>

- Splitting lenghty pages as per the previous tip makes it easier to create _workflow maps_ in your page. Now, you only need 2 things:
  2. Modify the boxes that are to appear in your page by copying and then modifying the contents of [this map file](https://github.com/h3abionet/H3ABionet-SOPs/blob/master/_includes/custom/usermap.html). You are only limited to 5 boxes in this simple design. For more complex maps (or more details, please refer to [the original documentation theme guide](https://idratherbewriting.com/documentation-theme-jekyll/mydoc_workflow_maps.html).

  1. In your intended page, add the following lines to your front matter:
```
simple_map: true
map_name: usermap  // this is the name of your map file 
 box_number: 1      // this is the order of the box in your map
```

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
