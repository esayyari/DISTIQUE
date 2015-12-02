#!/bin/bash
DIR_Astral=$WS_HOME/ASTRAL/Astral
path_11=/oasis/projects/nsf/uot136/esayyari/data/mammalian
out=/oasis/projects/nsf/uot136/esayyari/outputs/
if [ -s $DIR/tasks.massive ]; then
rm $DIR/tasks.massive
fi
path_tmp=$path_11
for file in `find $path_tmp -type f -name "*.gt"`; do
	tmp_out=$( dirname "${file}" )
	tmp_out=$( echo $tmp_out | sed -e "s/^.*data\///" )
	
	f_tmp=$( basename "${file}" )
	if [ -d $out/$tmp_out/astral ]; then
		printf "the folder $out/$tmp_out/astral already exists\n"
	else
		printf "the folder $out/$tmp_out/astral created\n"
		mkdir -p $out/$tmp_out/astral
	fi
	if [  -s $out/$tmp_out/astral/distance.d_astral_tree.nwk ]; then 
	printf " mkdir -p ./$out/$tmp_out/astral;  /usr/bin/time -po ./$out/$tmp_out/astral/$f_tmp-log.info java -Xmx2000M -jar $DIR_Astral/astral.4.8.3.jar -t 2 -i $file -o ./$out/$tmp_out/astral/astral_tree.nwk -q $out/$tmp_out/astral/distance.d_astral_tree.nwk;cp -r ./$out/$tmp_out/astral/* $out/$tmp_out/astral/ \n">>tasks.masive
	fi
done
