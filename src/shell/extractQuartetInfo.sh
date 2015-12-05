#!/bin/bash
test $# == 2 || { echo USAGE: $0 file_with_quartet_info output_file; exit 1; }

DIR=$( cd "$(dirname ${BASH_SOURCE[0]})" && pwd)
cat $1 | grep -e '^{[0-9]' | grep -oe '\[.*' | sed -e 's/[{} ]//g' | sed -e 's/:.*//g' > "`echo $2`"-bipartition

cat $1 | grep -oe '^{[0-9].*' | sed -e 's/\[.*//g' | sed -e 's/ //g'  >  "`echo $2`"-quadpartition

cat $1 | grep -e '^{[0-9]' | grep -oe ':.*' | sed -e 's/ //g' | sed -e 's/://g' >  "`echo $2`"-posterior
