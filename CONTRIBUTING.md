# How to contribute

Thank you for your interest in contributing to this effort! This SOP is published openly with the aim of helping the greater bioinformatics community, so we love to recieve contributions. There are many ways to do so, from reporting typos, issues, suggesting content, ... etc.

 To make it easier, we created this guide- let's dive in!

# Coding/Naming conventions and style

# Documentation

# Road map (milestones we like to eventually achieve)

# PRs

# Changelog <probably this should go elsewhere?>

# Where to start/ Project organization:

(Specific for the documentation theme. This is a quick guide adopted from https://idratherbewriting.com/documentation-theme-jekyll/)

- To add ad-hoc pages, just place the markdown file in the root directory. Add a frontmatter similar to this at top (or copy-paste and modify the contents of the file `frontmatter.md`). Your content should follow.

```

---
title: <title of your page- will appear at the browser window>
tags: [formatting]
keywords: <relevant keywords to your content> 
last_updated: <date of last update- appears at the bottom of the pate>
summary: <summary of the page if desired> 
sidebar: <For sidebar navigation, create a *.yml navigation file in _data/sidebars and put its name here. Only 2 level navigation allowed>
hide_sidebar: true #Use this if you don't wish to include sidebar, or else, the default sops_sidebar.yml will appear
permalink: <>
---

<Your intended page content, formatted in markdown. >

```

- For the sidebar, it should refer to a file of the intended table of contents to be found in _data/sidebars  An example is sops_sidebar.yml. Simply, add this line to the frontmatter of your page: 

```
sidebar: mydoc_sidebar
```
