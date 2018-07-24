This website is based on the [Minimal Mistakes Jekyll theme](https://mmistakes.github.io/minimal-mistakes/). [Jekyll](https://jekyllrb.com/), written in Ruby, is a "simple, blog-aware, static site generator".

Besides viewing online, you may host a local copy of this repo and view from from your browser- which is helpful in proposing changes to the site.

## Locally building

1. Clone down the repository (`git clone https://github.com/HPCBio/H3ABionet-SOPs.git`) 
2. `cd` into the cloned directory
3. Assuming Ruby is installed[^1], at the prompt, simply run `bundle`
4. Run `bundle exec jekyll serve` to start the preview server
5. Visit [`localhost:4000`](http://localhost:4000) in your browser to preview the site 

## Notes

1. Jekyll is written in Ruby, and you need it installed beforehand. You need at least Ruby 2.1.0. The [Bundler](https://bundler.io/) may be installed using: `gem install bundler`

2. On some Ubuntu 16.04 systems, installing the bundler as above (`gem install bundler`), may give errors of the form: `ERROR:  While executing gem ... (Gem::FilePermissionError)`. Some useful pointers for handling this case can be found in the stack overflow entries [here](https://stackoverflow.com/questions/37720892/you-dont-have-write-permissions-for-the-var-lib-gems-2-3-0-directory?answertab=votes#tab-top)
 
