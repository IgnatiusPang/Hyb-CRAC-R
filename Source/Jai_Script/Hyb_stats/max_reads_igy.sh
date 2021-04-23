HYB=$1
PWS=$2
MAX_READS=$3
STATS_SCRIPT_LOC=$4

pws_len=$(tail -n 1 $PWS | cut -f 2)
num_cpu=$(nproc)

echo "inputs are: " $HYB " " $PWS

LOC_TEMP_DIR=$(dirname $HYB) #$TMPDIR

HYB_NAME=${HYB##*/}
PWS_NAME=${PWS##*/}

#convert to GTF format to access each hybrid half
${STATS_SCRIPT_LOC}/hyb2gtf.py $HYB > $LOC_TEMP_DIR/$HYB_NAME.gtf
echo $HYB_NAME "converted to gtf"

#sort hybrids and pws to allow file splitting
sort -k4n,4 $LOC_TEMP_DIR/$HYB_NAME.gtf > $LOC_TEMP_DIR/$HYB_NAME.gtf.sort
sort -k2n,2 $PWS > $LOC_TEMP_DIR/$PWS_NAME.sort

echo "files sorted" $PWS_NAME".sort " $HYB_NAME".gtf.sort"

${STATS_SCRIPT_LOC}/split_gtf.py -s 100000 -f $LOC_TEMP_DIR/$HYB_NAME.gtf.sort
${STATS_SCRIPT_LOC}/split_gtf.py -s 100000 -t pws -f $LOC_TEMP_DIR/$PWS_NAME.sort

echo "files split"

#create a counter and run through each chromosome fragment reporting max read depth at the hybrid coordinate
#fragment the file to reduces search space and speed up processing

${STATS_SCRIPT_LOC}/make_counter.py 100000 $pws_len > $LOC_TEMP_DIR/$HYB_NAME.counter.txt

# divide the job across available CPUs

split -a 4 -d -l $num_cpu $LOC_TEMP_DIR/$HYB_NAME.counter.txt $LOC_TEMP_DIR/$HYB_NAME.counter_

echo "counter made and running while loop on HYB.gtf and PWS"

for x in $LOC_TEMP_DIR/$HYB_NAME.counter_* ;
do 
  while read i ;
  do
    echo "processing " $LOC_TEMP_DIR/$HYB_NAME.gtf.sort$i.gtf $LOC_TEMP_DIR/$PWS_NAME.sort$i.pws
    ${STATS_SCRIPT_LOC}/hyb_gtf2max_read_pws.py $LOC_TEMP_DIR/$HYB_NAME.gtf.sort$i.gtf $LOC_TEMP_DIR/$PWS_NAME.sort$i.pws &
  done < $x
wait
done

wait

echo "while loop finished"

#put the file fragments back together and remove header line (ok)
cat $LOC_TEMP_DIR/$HYB_NAME.gtf.sort*.gtf_correlate.txt | awk ' !/#/ {print} ' | sort -k5,5 -k7n,7 -k9,9 > $LOC_TEMP_DIR/$HYB_NAME.correlate.txt

echo "files reassembled"

#reassemble hybrid halfs based on ID (ok)
awk ' fnd[$1] {print fnd[$1] "\t"  $0 } !fnd[$1] {fnd[$1] = $0 } ' $LOC_TEMP_DIR/$HYB_NAME.correlate.txt > $LOC_TEMP_DIR/$HYB_NAME.correlate.reassembled.txt

echo "hybrid halfs reassembled"

# make a header

echo -e "ID\tmax_reads\thyb_count\trna_class\tchromo\tname\tstart\tend\tstrand\tID\tmax_reads\thyb_count\trna_class\tchromo\tname\tstart\tend\tstrand\tmin\tmax" > $MAX_READS

#add min and max fields
awk '{ $2 > $11 ? max = $2 : max = $11 ; $2 > $11 ? min = $11 : min = $2 ; print $0 "\t" min "\t" max }' $LOC_TEMP_DIR/$HYB_NAME.correlate.reassembled.txt >> $MAX_READS

echo "min and max variables added"
echo "time for a tidy up"

#tidy up step to remove intermediate files
rm $LOC_TEMP_DIR/$HYB_NAME.counter*
rm $LOC_TEMP_DIR/$HYB_NAME.gtf
rm $LOC_TEMP_DIR/$HYB_NAME.gtf.sort
rm $LOC_TEMP_DIR/$HYB_NAME.gtf.sort*.gtf
rm $LOC_TEMP_DIR/$PWS_NAME.sort
rm $LOC_TEMP_DIR/$PWS_NAME.sort*.pws
rm $LOC_TEMP_DIR/$HYB_NAME.gtf.sort*.gtf_correlate.txt
rm $LOC_TEMP_DIR/$HYB_NAME.correlate.txt
rm $LOC_TEMP_DIR/$HYB_NAME.correlate.reassembled.txt

