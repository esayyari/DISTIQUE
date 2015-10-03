#!/bin/bash
DIR=/home/esayyari/repository/DISTIQUE/src/shell
DIR=/home/esayyari/repository/ASTRID/
DIR=/home/esayyari/repository/ASTRAL/Astral
path=/oasis/projects/nsf/uot136/esayyari/data/mammalian/
method=prod
summethod=fm
out=/oasis/projects/nsf/uot136/esayyari/outputs/
if [ -s $DIR/tasks.massive ]; then
rm $DIR/tasks.massive
fi
for file in `find $path -type f -name "genetrees.gt"`; do
tmp_out=$( dirname "${file}" )
tmp_out=$( echo $tmp_out | sed -e "s/^.*data//g" )
f_tmp=$( basename "${file}" )
if [ -d $out/$tmp_out/astral ]; then
	printf "the folder $out/$tmp_out/astral already exists\n"
else
	printf "the folder $out/$tmp_out/astral created\n"
	mkdir -p $out/$tmp_out/astral
fi
if [ "$summethod" == "fm" -a -s $out/$tmp_out/astral/distance.d ]; then
echo the results for $file exists
else
printf "/usr/bin/time -po $out/$tmp_out/astral/$f_tmp-log.info java -Xmx25000M -jar $DIR/astral.4.7.8.jar -i $file -o $out/$tmp_out/astral/distance.d_astral_tree.nwk\n" >> tasks.massive
#printf "/usr/bin/time -po $out/$tmp_out/$f_tmp-log.info $DIR/summarize.tree.sh -s $summethod -d $method -w 1 -r $out/$tmp_out $file\n">>tasks.massive
#printf "for i in {1..100}; do sleep 10; echo $i; done\n" >> tasks.massive
fi
done
