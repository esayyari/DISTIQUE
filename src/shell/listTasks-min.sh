#!/bin/bash

DIR_DISTIQUE=$WS_HOME/DISTIQUE/src/shell
path_11=/oasis/projects/nsf/uot136/esayyari/data/avian
method=min
summethod=fm
out=/oasis/projects/nsf/uot136/esayyari/outputs/
path_tmp=$path_11
for file in `find $path_tmp -type f -name "*.gt"`; do
        tmp_out=$( dirname "${file}" )
        tmp_out=$( echo $tmp_out | sed -e "s/^.*data\///" )
        f_tmp=$( basename "${file}" )
        if [ -d $out/$tmp_out/distique-1/$method ]; then
                printf "the folder $out/$tmp_out/distique-1/$method already exists\n"
        else
                printf "the folder $out/$tmp_out/distique-1/$method created\n"
                mkdir -p $out/$tmp_out/distique-1/$method
        fi
        if [ ! -s $out/$tmp_out/distique-cons/$method/distance.d_distique-1_tree.nwk ]; then
                printf "mkdir -p ./$out/$tmp_out/distique-1/$method/; source /etc/profile.d/modules.sh; module load python; module load scipy; VIRTUAL_ENV_DISABLE_PROMPT=1; source $WS_HOME/python27-gordon/bin/activate; /usr/bin/time -po ./$out/$tmp_out/distique-1/$mehtod/$f_tmp-log.info $DIR_DISTIQUE/summarize.tree.sh -s fm -d $method -w 1  -r ./$out/$tmp_out/distique-1/$method $file; cp -r ./$out/$tmp_out/distique-1/$method $out/$tmp_out/distique-1/$method\n">>tasks.massive-$method.3
        fi
done
