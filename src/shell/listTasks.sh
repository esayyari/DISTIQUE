#!/bin/bash
DIR_DISTIQUE=$WS_HOME/DISTIQUE/src/shell
DIR_NJst=$WS_HOME/ASTRID/src
DIR_Astral=$WS_HOME/ASTRAL/Astral
path_11=/oasis/projects/nsf/uot136/esayyari/data/11-taxon
DIR_UTILS=$WS_HOME/DISTIQUE/src/utils
method=prod
summethod=fm
out=/oasis/projects/nsf/uot136/esayyari/outputs/
out=/oasis/scratch/esayyari/temp_project/oasis/projects/nsf/uot136/esayyari/outputs/
if [ -s $DIR/tasks.massive ]; then
rm $DIR/tasks.massive
fi
path_tmp=$path_11
for file in `find $path_tmp -type f -name "*.gt"`; do
	tmp_out=$( dirname "${file}" )
	tmp_out=$( echo $tmp_out | sed -e "s/^.*data\///" )
	
	find $out/$tmp_out/distique-2/ -name "distance.d_*_tree.nwk" | wc -l
	f_tmp=$( basename "${file}" )
	if [ -d $out/$tmp_out/distique-2/gmean/$method ]; then
		printf "the folder $out/$tmp_out/distique-2/gmean/$method already exists\n"
	else
		printf "the folder $out/$tmp_out/distique-2/gmean/$method created\n"
		mkdir -p $out/$tmp_out/distique-2/gmean/$method
	fi
	if [ -d $out/$tmp_out/distique-2/mean/$method ]; then
		printf "the folder $out/$tmp_out/distique-2/mean/$method already exists\n"
	else
		printf "the folder $out/$tmp_out/distique-2/mean/$method created\n"
		mkdir -p $out/$tmp_out/distique-2/mean/$method
	fi
	if [ -d $out/$tmp_out/distique-2/rmsquare/$method ]; then
		printf "the folder $out/$tmp_out/distique-2/rmsquare/$method already exists\n"
	else
		printf "the folder $out/$tmp_out/distique-2/rmsquare/$method created\n"
		mkdir -p $out/$tmp_out/distique-2/rmsquare/$method
	fi
	if [ -d $out/$tmp_out/distique-2/mean-low-thr/$method ]; then
		printf "the folder $out/$tmp_out/distique-2/mean-low-thr/$method already exists\n"
	else
		printf "the folder $out/$tmp_out/distique-2/mean-low-thr/$method created\n"
		mkdir -p $out/$tmp_out/distique-2/mean-low-thr/$method
	fi
	if [ -d $out/$tmp_out/distique-2/mean-vlow-thr/$method ]; then
		printf "the folder $out/$tmp_out/distique-2/mean-vlow-thr/$method already exists\n"
	else
		printf "the folder $out/$tmp_out/distique-2/mean-vlow-thr/$method created\n"
		mkdir -p $out/$tmp_out/distique-2/mean-vlow-thr/$method
	fi
#	if [ -d $out/$tmp_out/astral ]; then
#		printf "the folder $out/$tmp_out/astral already exists\n"
#	else
#		printf "the folder $out/$tmp_out/astral created\n"
#		mkdir -p $out/$tmp_out/astral
#	fi
#	if [ -d $out/$tmp_out/njst ]; then
#		printf "the folder $out/$tmp_out/njst already exists\n"
#	else
#		printf "the folder $out/$tmp_out/njst created\n"
#		mkdir -p $out/$tmp_out/njst
#	fi
	q=$(echo $file | sed -e 's/\/genetrees.gt//g')
	if [ -s $q/quartets.q ]; then
		printf "found\n"
		qfile=$q/quartets.q
	elif [ -s $q/quartet_tmp.q ]; then
		printf "found\n"
		qfile=$q/quartet_tmp.q
	else
		exit 1
	fi
	if [ ! -s $out/$tmp_out/distique-2/gmean/$method/distance.d_distique_tree.nwk ]; then
		printf "mkdir -p /oasis/scratch/esayyari/temp_project/$out; mkdir -p ./$out/$tmp_out/distique-2/gmean/$method/; source /etc/profile.d/modules.sh; module load python; module load scipy; VIRTUAL_ENV_DISABLE_PROMPT=1; source $WS_HOME/python27-gordon/bin/activate; /usr/bin/time -po ./$out/$tmp_out/distique-2/gmean/$mehtod/$f_tmp-log.info $DIR_UTILS/distique-2.py -a gmean -f $qfile -g $file -m $method  -o ./$out/$tmp_out/distique-2/gmean/$method; cp -r ./$out/* /$out/ \n">>tasks.massive
	fi
	if [ ! -s $out/$tmp_out/distique-2/mean/$method/distance.d_distique_tree.nwk ]; then
               printf "mkdir -p /oasis/scratch/esayyari/temp_project/$out; mkdir -p ./$out/$tmp_out/distique-2/mean/$method/; source /etc/profile.d/modules.sh; module load python; module load scipy; VIRTUAL_ENV_DISABLE_PROMPT=1; source $WS_HOME/python27-gordon/bin/activate; /usr/bin/time -po ./$out/$tmp_out/distique-2/mean/$mehtod/$f_tmp-log.info $DIR_UTILS/distique-2.py -a mean -f $qfile -g $file -m $method  -o ./$out/$tmp_out/distique-2/mean/$method; cp -r ./$out/* /$out \n">>tasks.massive
        fi
	if [ ! -s $out/$tmp_out/distique-2/rmsquare/$method/distance.d_distique_tree.nwk ]; then
                printf "mkdir -p /oasis/scratch/esayyari/temp_project/$out; mkdir -p ./$out/$tmp_out/distique-2/rmsquare/$method/; source /etc/profile.d/modules.sh; module load python; module load scipy; VIRTUAL_ENV_DISABLE_PROMPT=1; source $WS_HOME/python27-gordon/bin/activate; /usr/bin/time -po ./$out/$tmp_out/distique-2/rmsquare/$mehtod/$f_tmp-log.info $DIR_UTILS/distique-2.py -f $qfile -a rmsquare -g $file -m $method  -o ./$out/$tmp_out/distique-2/rmsquare/$method; cp -r ./$out/* /$out \n">>tasks.massive
        fi
	if [ ! -s $out/$tmp_out/distique-2/mean-low-thr/$method/distance.d_distique_tree.nwk ]; then
                printf "mkdir -p /oasis/scratch/esayyari/temp_project/$out; mkdir -p ./$out/$tmp_out/distique-2/mean-low-thr/$method/; source /etc/profile.d/modules.sh; module load python; module load scipy; VIRTUAL_ENV_DISABLE_PROMPT=1; source $WS_HOME/python27-gordon/bin/activate; /usr/bin/time -po ./$out/$tmp_out/distique-2/mean-low-thr/$mehtod/$f_tmp-log.info $DIR_UTILS/distique-2.py -f $qfile -a mean -t 0.4 -g $file -m $method  -o ./$out/$tmp_out/distique-2/mean-low-thr/$method; cp -r ./$out/* /$out \n">>tasks.massive
        fi
	if [ ! -s $out/$tmp_out/distique-2/mean-vlow-thr/$method/distance.d_distique_tree.nwk ]; then
                printf "mkdir -p /oasis/scratch/esayyari/temp_project/$out; mkdir -p ./$out/$tmp_out/distique-2/mean-vlow-thr/$method/; source /etc/profile.d/modules.sh; module load python; module load scipy; VIRTUAL_ENV_DISABLE_PROMPT=1; source $WS_HOME/python27-gordon/bin/activate; /usr/bin/time -po ./$out/$tmp_out/distique-2/mean-vlow-thr/$mehtod/$f_tmp-log.info $DIR_UTILS/distique-2.py -f $qfile -a mean -t 0.1 -g $file -m $method  -o ./$out/$tmp_out/distique-2/mean-vlow-thr/$method; cp -r ./$out/* /$out \n">>tasks.massive
        fi
	#if [ ! -s $out/$tmp_out/distance.d_njst_tree.nwk ]; then
	#printf "mkdir -p /oasis/scratch/esayyari/temp_project/$out; mkdir -p ./$out/$tmp_out/njst; source /etc/profile.d/modules.sh; module load python; module load scipy; VIRTUAL_ENV_DISABLE_PROMPT=1; source $WS_HOME/python27-gordon/bin/activate; /usr/bin/time -po $out/$tmp_out/njst/$f_tmp-log.info python $DIR_NJst/ASTRID.py -i $file -m fastme2 -o ./$out/$tmp_out/njst/distance.d_njst_tree.nwk -c ./$out/$tmp_out/njst/CACHE.csv; cp -r ./$out/* /oasis/scratch/esayyari/temp_project/cp -r ./$out/* /oasis/scratch/esayyari/temp_project/$out/ \n" >>tasks.massive9
	#fi
	#if [ ! -s $out/$tmp_out/astral/distance.d_astral_tree.nwk ]; then 
	#printf "mkdir -p /oasis/scratch/esayyari/temp_project/$out; mkdir -p ./$out/$tmp_out/astral;  /usr/bin/time -po ./$out/$tmp_out/astral/$f_tmp-log.info java -Xmx2000M -jar $DIR_Astral/astral.4.7.8.jar -i $file -o ./$out/$tmp_out/astral/distance.d_astral_tree.nwk;cp -r ./$out/* /oasis/scratch/esayyari/temp_project/$out \n">>tasks.massive-astral.7-$method
	#fi
done
