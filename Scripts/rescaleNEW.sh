rm *.rescaled.eps

# FIND MAX
echo "FIND MAX"
max="0"
for i in `ls $@`;
do
	temp=`LC_ALL=en_US awk 'BEGIN{max=0;} {if($4=="setrgbcolor"){ if(max<$1)max=$1;}} END{print max}' < $i`
	if [ `echo "$max < $temp" | bc` -eq 1 ] 
	then
		max=$temp
	fi
done
echo "Max. Intensity is $max"


# RESCALE
echo "RESCALE & RENAME (leading zeros)"
index=0
for i in `ls $@`;
do
	# extract numeration from filename
	index=$(echo "$i" | sed "s/[^0-9]//g")

	# create new filename with new numeration (leading zeros)
	if   [ ${index} -lt 10   ] ; then filename=markerMRI0000${index}.rescaled
	elif [ ${index} -lt 100  ] ; then filename=markerMRI000${index}.rescaled
	elif [ ${index} -lt 1000 ] ; then filename=markerMRI00${index}.rescaled
	elif [ ${index} -lt 10000 ]; then filename=markerMRI0${index}.rescaled
	else filename=markerMRI${index}.rescaled
	fi 
	
	# rescale color intensities (.eps)
	echo "$i > ${filename}.eps"
	LC_ALL=en_US awk '{if($4=="setrgbcolor"){ if(max<$1)max=$1; printf("%f %f %f %s\n",$1/max,$2/max,$3/max,$4)}else{print }}' "max=$max" $i > ${filename}.eps

	# convert to .png
<<<<<<< .mine
        convert ${filename}.eps -define png:bit-depth=16 -define png:color-type=Grayscale ${filename}.png
=======
	convert ${filename}.eps -colorspace Gray -define png:bit-depth=16 ${filename}.png
>>>>>>> .r30
	index=`expr ${index} + 1`
done

# CONVERT TO MOVIE
#convert -delay 10 -size 500x500 markerMRI?????.rescaled.eps markerMRI_int.mpeg
echo "CONVERT TO MOVIE (interpolated)"
convert -delay 10 -resize 500x500 markerMRI*.rescaled.eps markerMRI_int.mpeg
echo "CONVERT TO MOVIE (pixeled)"
convert -delay 10 -resize 500x500 -density 500 markerMRI*.rescaled.eps markerMRI.mpeg


# RESCALE
#echo "RESCALE"
#for i in `ls $@`;
#do
#	LC_ALL=en_US awk '{if($4=="setrgbcolor"){ if(max<$1)max=$1; printf("%lf %lf %lf %s\n",$1/max,$2/max,$3/max,$4)}else{print }}' "max=$max" $i > ${i}.rescaled.eps
#done
#convert -delay 10 -size 500x500 markerMRI?.rescaled.eps markerMRI??.rescaled.eps markerMRI???.rescaled.eps markerMRI_int.mpeg
