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
