This website is based on the [Documentation Jekyll theme](http://idratherbewriting.com/documentation-theme-jekyll/). [Jekyll](https://jekyllrb.com/), written in Ruby, is a "simple, blog-aware, static site generator".


The live version of this website is accessible via: [https://h3abionet.github.io/H3ABionet-SOPs/](https://h3abionet.github.io/H3ABionet-SOPs/). However, you may host a local copy of this repo and view from from your browser- which is helpful in proposing changes to the site (See following sections).

## Locally building

1. Clone down the repository (`git clone https://github.com/h3abionet/H3ABionet-SOPs.git`) 
2. `cd` into the cloned directory
3. Assuming Ruby is installed<sup>[1](#footnote1)</sup>, at the prompt, simply run `bundle`<sup>[2](#footnote2)</sup>.
4. Run `bundle exec jekyll serve` to start the preview server
5. Visit [`localhost:4000`](http://localhost:4000) in your browser to preview the site

## Contributing

Feedback is most welcome, and so are contributions of any kind (e.g. reporting typos, issues, suggesting content, ... etc). Please see our [contribution guide](https://h3abionet.github.io/H3ABionet-SOPs/cont_tech-guide-1) and the [comprehensive contribution guide](https://h3abionet.github.io/H3ABionet-SOPs/contributing) intended for those not very comfortable with git/jekyll 

## Notes

<a name="footnote1">1</a>: Jekyll is written in Ruby, and you need it installed for it to function (at least Ruby 2.1.0 is needed). On Ubuntu systems, you may run `ruby -v` at the prompt to check Ruby and its installed version.

<a name="footnote2">2</a>: The [Bundler](https://bundler.io/) may be installed using: `gem install bundler`. On some Ubuntu 16.04 systems, installing the bundler as above (`gem install bundler`), may give errors of the form: `ERROR:  While executing gem ... (Gem::FilePermissionError)`. Some useful pointers for handling this case can be found in the stack overflow entries [here](https://stackoverflow.com/questions/37720892/you-dont-have-write-permissions-for-the-var-lib-gems-2-3-0-directory?answertab=votes#tab-top)
