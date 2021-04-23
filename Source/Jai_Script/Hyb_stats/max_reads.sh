#!/bin/bash
#
# copyright 2017 Jai J. Tree 
# GNU Lesser General Public License
#
# shell script to retrieve max read depth @ hybrid coordinates
# inputs are: HYB PWS
# max_reads.sh HYB PWS


HYB=$1 
PWS=$2
pws_len=$(tail -n 1 $PWS | cut -f 2)
num_cpu=$(nproc)

echo "inputs are: " $HYB " " $PWS

#convert to GTF format to access each hybrid half
hyb2gtf.py $HYB > $HYB.gtf

echo $HYB "converted to gtf"

#sort hybrids and pws to allow file splitting
sort -k4n,4 $HYB.gtf > $HYB.gtf.sort
sort -k2n,2 $PWS > $PWS.sort

echo "files sorted" $PWS".sort " $HYB".gtf.sort"

split_gtf.py -s 100000 -f $HYB.gtf.sort
split_gtf.py -s 100000 -t pws -f $PWS.sort

echo "files split"

#create a counter and run through each chromosome fragment reporting max read depth at the hybrid coordinate
#fragment the file to reduces search space and speed up processing

make_counter.py 100000 $pws_len > counter.txt

# divide the job across available CPUs

split -a 4 -d -l $num_cpu counter.txt counter_

echo "counter made and running while loop on HYB.gtf and PWS"
 
for x in counter_* ; 
	do 
	while read i ; 
		do 	
		echo "processing " $HYB.gtf.sort$i.gtf $PWS.sort$i.pws 
		hyb_gtf2max_read_pws.py $HYB.gtf.sort$i.gtf $PWS.sort$i.pws &
		done < $x 
	wait
	done 

wait

echo "while loop finished"

#put the file fragments back together and remove header line
cat $HYB.gtf.sort*.gtf_correlate.txt | awk ' !/#/ {print} ' | sort -k5,5 -k7n,7 -k9,9 > $HYB.correlate.txt

echo "files reassembled"

#reassemble hybrid halfs based on ID
awk ' fnd[$1] {print fnd[$1] "\t"  $0 } !fnd[$1] {fnd[$1] = $0 } ' $HYB.correlate.txt > $HYB.correlate.reassembled.txt

echo "hybrid halfs reassembled"

# make a header

echo -e "ID\tmax_reads\thyb_count\trna_class\tchromo\tname\tstart\tend\tstrand\tID\tmax_reads\thyb_count\trna_class\tchromo\tname\tstart\tend\tstrand\tmin\tmax" > $HYB.max_reads.txt

#add min and max fields
awk '{ $2 > $11 ? max = $2 : max = $11 ; $2 > $11 ? min = $11 : min = $2 ; print $0 "\t" min "\t" max }' $HYB.correlate.reassembled.txt >> $HYB.max_reads.txt

echo "min and max variables added"
echo "time for a tidy up"

#tidy up step to remove intermediate files
rm counter*
rm $HYB.gtf.sort*.gtf
rm $PWS.sort*.pws
rm $HYB.gtf.sort*.gtf_correlate.txt
rm $HYB.correlate.txt
rm $HYB.correlate.reassembled.txt
