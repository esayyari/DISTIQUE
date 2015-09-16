#!/bin/bash


shopt -s expand_aliases
source ~/.bash_profile
DIR=$( dirname "${BASH_SOURCE[0]}" )
source $DIR/setup.sh
name=$( basename "${1}")
echo $name
path=~/data/results
cat $1* > $path/$name-res.txt


