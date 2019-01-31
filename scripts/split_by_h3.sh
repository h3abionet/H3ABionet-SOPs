#!/bin/bash

## Usage: ./split_by_h3.sh input_file.md
##
## input_file.md a properly formatted markdown file of a working page 
### - It starts with proper frontmatter, and 
### - is divided into h2 and h3 headers. 
### - The last header section contains References and abbreviations of the page  (The References section is named "Bibliography")

## The result is a group of files split by level 3 headers
### - all files are consecutively named after the basename of input_file
### - frontmatter of all files is numbered and named accordingly
### - References and abbreviations section appended to each

## You may wish to delete input_file.md when done!

filename=$1

cp ${filename} bak.${filename}

sed -n '/^---$/,/^---$/p' ${filename} > yml.header
sed -i '/^---$/,/^---$/d' ${filename} 

sed -n '/## Bib.*/,//p' ${filename}> bib_abbrevs
sed -i '/## Bib.*/,//d' ${filename}

csplit --prefix ${filename%.*}- -n 1 ${filename}  '/^### .*$/' '{*}'

subfiles=`expr $( ls -1 ${filename%.*}-*   | wc -l) - 0`

for i in `seq 0 $subfiles` ; do
    cat yml.header ${filename%.*}-$i > ${filename%.*}-$i.md  #Copy header
    cat bib_abbrevs >> ${filename%.*}-$i.md #Copy abbreviations
    rm  ${filename%.*}-$i
    sed -i "s/${filename%.*}.html/${filename%.*}-$i.html/g" ${filename%.*}-$i.md #change permalink in yml frontmatter
    sed -i "s/toc: true/toc: false/g" ${filename%.*}-$i.md #remove table of content from frontmatter
    sed -i "s/summary: .*//g" ${filename%.*}-$i.md #remove summary from frontmatter 

done

rm bak.${filename} ${filename%.*}-$subfiles.md yml.header bib_abbrevs  ${filename}
