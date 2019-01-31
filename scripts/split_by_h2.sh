#!/bin/bash

## Usage: ./split_by_h2.sh input_file.md
##
## input_file.md a properly formatted markdown file of a working page 
### - It starts with proper frontmatter, and 
### - is divided into h2 and h3 headers. 
### - The last header section contains References and abbreviations of the page 

## The result is a group of files split by level 2 headers
### - all files are consecutively named after the basename of input_file
### - frontmatter of all files is numbered and named accordingly
### - References and abbreviations section appended to each

## You may wish to delete input_file.md when done!

filename=$1

cp ${filename} bak.${filename}

csplit --prefix ${filename%.*}- -n 1 ${filename}  '/^## .*$/' '{*}'

subfiles=$( ls -1 ${filename%.*}-*   | wc -l)

for i in `seq 1 $subfiles` ; do
    cat ${filename%.*}-0 ${filename%.*}-$i > ${filename%.*}-$i.md  #Copy header
    cat ${filename%.*}-`expr $subfiles - 1` >> ${filename%.*}-$i.md #Copy abbreviations
    rm  ${filename%.*}-$i
    sed -i "s/${filename%.*}.html/${filename%.*}-$i.html/g" ${filename%.*}-$i.md
done



rm bak.${filename} ${filename%.*}-0 ${filename%.*}-`expr $subfiles - 1`.md ${filename%.*}-`expr $subfiles`.md 
