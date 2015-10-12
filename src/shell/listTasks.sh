#!/bin/bash
DIR_DISTIQUE=$WS_HOME/DISTIQUE/src/shell
DIR_NJst=$WS_HOME/ASTRID/src
DIR_Astral=$WS_HOME/ASTRAL/Astral
path_11=/oasis/projects/nsf/uot136/esayyari/data/11-taxon
DIR_UTILS=$WS_HOME/DISTIQUE/src/utils
method=prod
summethod=fm
out=/oasis/projects/nsf/uot136/esayyari/outputs/

if [ -s $DIR/tasks.massive ]; then
rm $DIR/tasks.massive
fi
path_tmp=$path_11
for file in `find $path_tmp -type f -name "*.gt"`; do
	tmp_out=$( dirname "${file}" )
	tmp_out=$( echo $tmp_out | sed -e "s/^.*data\///" )
	f_tmp=$( basename "${file}" )
	if [ -d $out/$tmp_out/distique-cons/$method ]; then
		printf "the folder $out/$tmp_out/distique-cons/$method already exists\n"
	else
		printf "the folder $out/$tmp_out/distique-cons/$method created\n"
		mkdir -p $out/$tmp_out/distique-cons/$method
	fi
	if [ -d $out/$tmp_out/astral ]; then
		printf "the folder $out/$tmp_out/astral already exists\n"
	else
		printf "the folder $out/$tmp_out/astral created\n"
		mkdir -p $out/$tmp_out/astral
	fi
	if [ -d $out/$tmp_out/njst ]; then
		printf "the folder $out/$tmp_out/njst already exists\n"
	else
		printf "the folder $out/$tmp_out/njst created\n"
		mkdir -p $out/$tmp_out/njst
	fi
#	if [ ! -s $out/$tmp_out/astral/distance.d_astral_tree.nwk ]; then 
#	printf "source /etc/profile.d/modules.sh; module load python; module load scipy; VIRTUAL_ENV_DISABLE_PROMPT=1; source $WS_HOME/python27-gordon/bin/activate; /usr/bin/time -po $out/$tmp_out/astral/$f_tmp-log.info java -Xmx2000M -jar $DIR_Astral/astral.4.7.8.jar -i $file -o $out/$tmp_out/astral/distance.d_astral_tree.nwk\n" >> tasks.massive7
#	fi
	if [ ! -s $out/$tmp_out/distique-cons/$method/distance.d_distique_tree.nwk ]; then
		printf "source /etc/profile.d/modules.sh; module load python; module load scipy; VIRTUAL_ENV_DISABLE_PROMPT=1; source $WS_HOME/python27-gordon/bin/activate; /usr/bin/time -po $out/$tmp_out/$f_tmp-log.info $DIR_UTILS/distique.py -g $file -m $method  -o $out/$tmp_out/distique-cons/$method \n">>tasks.massive-$method.2
	fi
#	if [ ! -s $out/$tmp_out/distance.d_njst_tree.nwk ]; then
#	printf "source /etc/profile.d/modules.sh; module load python; module load scipy; VIRTUAL_ENV_DISABLE_PROMPT=1; source $WS_HOME/python27-gordon/bin/activate; /usr/bin/time -po $out/$tmp_out/njst/$f_tmp-log.info python $DIR_NJst/ASTRID.py -i $file -m fastme2 -o $out/$tmp_out/njst/distance.d_njst_tree.nwk -c $out/$tmp_out/njst/CACHE.csv\n" >>tasks.massive9
#	fi
done
