# FIND MAX
max="0"
for i in `ls $@`;
do
	#echo $i
	temp=`LC_ALL=en_US awk 'BEGIN{max=0;} {if($4=="setrgbcolor"){ if(max<$1)max=$1;}} END{print max}' < $i`
	#echo $temp
	#echo "$max < $temp" | bc
	if [ `echo "$max < $temp" | bc` -eq 1 ] 
	then
		max=$temp
	fi
	echo $max
done
echo $max

# RESCALE
index=0
for i in `ls $@`;
do
	if   [ $x -lt 10   ] ; then filename=markerMRI0000${index}.rescaled
	elif [ $x -lt 100  ] ; then filename=markerMRI000${index}.rescaled
	elif [ $x -lt 1000 ] ; then filename=markerMRI00${index}.rescaled
	elif [ $x -lt 10000 ]; then filename=markerMRI0${index}.rescaled
	else filename=markerMRI${index}.rescaled
	fi 

	LC_ALL=en_US awk 'BEGIN{max=4000} {if($4=="setrgbcolor"){ if(max<$1)max=$1; printf("%lf %lf %lf %s\n",$1/max,$2/max,$3/max,$4)}else{print }}' max=1 $i > ${filename}.eps

	convert ${filename}.eps ${filename}.png
done
