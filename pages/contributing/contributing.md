---
title: Technical Guide
keywords: contribution, guide, tips
tags: [documentation, formatting]
last_updated: May 23, 2020
toc: false
sidebar: cont_sidebar
permalink: contributing.html
folder: contributing
author_profile: true
authors:
 - Edward Lukyamuzi
---

## Contributing to SOPS

H3ABioNet welcomes contributions from everyone. Here are a few guidelines and instructions if you're thinking of helping with the curation of SOPS.

## Getting started

### Install jekyll
* Make sure you have jekyll 2.4.0 or greater installed. You can find information on how to install jekyll [here](https://jekyllrb.com/docs/)

### Fork project repo
* After following the [jekyll installation guide](https://jekyllrb.com/docs/installation/), Go to the H3ABionet-SOPs repository on Github [here](https://github.com/h3abionet/H3ABionet-SOPs) and click on Fork. Fork will create a copy of the repository in your Github account so that you can make changes to the project

### Clone repo
* Go to your own Github repository and you will see a repository named H3ABionet-SOPs. Make a local copy of it on your computer by hitting the clone button
* Copy the Https URL, Open your terminal or git bash window and Change directory to where you want to create a copy of this project 
* Run `git clone` appending the Https URL, looks like `git clone https://github.com/userid/H3ABionet-SOPs.git` where `userid` is your github username. This will automatically clone this repository. You now have a local copy on your computer.The working directory on your computer will have the same name as the remote repository on Github, H3ABionet-SOPs
* Add your fork to git's remotes:
  * If you use SSH authentication: `git remote add <your username> git@github.com:<your username>/H3ABionet-SOPs.git`.
  * Otherwise: `git remote add <your username>https://github.com/<your username>/H3ABionet-SOPs.git`.

### Make edits to the SOPs

* Switch into S0Ps directory, `cd H3ABionet-SOPs`.
* Using your favorite text editor make changes to the code and files as you deem necessary

### Test your changes locally

* Build a local copy of site in your computer and view from your own browser as below:
  * At the prompt, simply run `bundle`
  * Run `bundle exec jekyll serve` to start the preview server
  * Visit `localhost:4000` in your browser to preview the site

### Commit

* Run `git init` command
* Run `git status` to see the current status of your local repo, showing that there are changes not staged to be committed with the modified files in red
* Run `git add .` to add these changes
* Re-run `git status`, your changes have now been staged and are ready to be committed
* Now run git commit -m to commit these changes. Remember to describe the specific changes you made in the git commit message
* Run git status again to check that everything is up to date

### Push changes to your fork

* We can now push our code to the forked repository. Before we do this, run `git remote -v` to crosscheck that the remote repository on Github our local repository is connected to is indeed the forked repo
* Run `git push origin master` to push our changes to this repository
* Check that your remote repository has updated on Github

### Submit Pull Requests

* You’re now all ready to submit the improvements/changes you’ve made to the project’s maintainers for approval
* **Create a pull request against `master`**, and a contributor will come by and review your submission. They may ask for some changes, and hopefully your contribution will be merged to the `master` branch!

## Communication

If you need any help contributing to SOPs, several [contributors](https://github.com/h3abionet/H3ABionet-SOPs/graphs/contributors) are available
