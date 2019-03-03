---
title: Technical Guide
keywords: contribution, guide, tips
tags: [documentation, formatting]
last_updated: 1 2, 2019
sidebar: cont_sidebar
permalink: cont_tech-guide-1.html  
folder: contributing
author_profile: true
authors:
 - Azza_Ahmed
---

## How to contribute {#overview}

Thank you for your interest in contributing to this effort! These SOPs are published openly with the aim of helping the greater bioinformatics community, so we love to recieve contributions. There are many ways to do so, from reporting typos, issues, suggesting content, ... etc.

 To make it easier, we created this guide- let's dive in!

#### Credit and Contributions

We find it important that all contributors get credit for their work. Along with including your name in the `yml` header of your contributed effort, please make sure you also update the [authors.yml](https://github.com/h3abionet/H3ABionet-SOPs/blob/master/authors.yml) to include your contact details. 

## Inner working of the site
This website is based on the [Documentation Jekyll theme](http://idratherbewriting.com/documentation-theme-jekyll/). [Jekyll](https://jekyllrb.com/), written in Ruby, is a "simple, blog-aware, static site generator".

Before proposing any changes, test your changes first by building a local copy in your computer and viewing from your own browser as below:

1. Clone down the repository (`git clone https://github.com/h3abionet/H3ABionet-SOPs.git`)
2. `cd` into the cloned directory
3. Assuming Ruby is installed[^1], at the prompt, simply run `bundle`[^2].
4. Run `bundle exec jekyll serve` to start the preview server
5. Visit [`localhost:4000`](http://localhost:4000) in your browser to preview the site


### Road map {#road-map}

- ~~Split lenghty pages into smaller chunks~~
- ~~Create diagrams for SOP pages~~
- Add authors' contact details next to the pages they contribute to
- ~~Generate pdf files of the SOPs~~: Requires some commercial license, so not going to happen
- ~~Integrate with Disqus for easy commenting~~
- ~~Add license information to pages~~
- ~~Add proper H3ABioNet branding~~

#### Notes

[^1]: Jekyll is written in Ruby, and you need it installed for it to function (at least Ruby 2.1.0 is needed). On Ubuntu systems, you may run `ruby -v` at the prompt to check Ruby and its installed version.

[^2]: The [Bundler](https://bundler.io/) may be installed using: `gem install bundler`. On some Ubuntu 16.04 systems, installing the bundler as above (`gem install bundler`), may give errors of the form: `ERROR:  While executing gem ... (Gem::FilePermissionError)`. Some useful pointers for handling this case can be found in the stack overflow entries [here](https://stackoverflow.com/questions/37720892/you-dont-have-write-permissions-for-the-var-lib-gems-2-3-0-directory?answertab=votes#tab-top)
