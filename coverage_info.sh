# /usr/local/bin/bash
inbam=$1
OUT=$2 
tmpfile=tmpfile
touch $tmpfile

samtools flagstat $inbam > $tmpfile

echo -e "Path_file,File,Total_reads,Supplementary_reads,Duplicate_reads,Paired_reads" > $OUT
TOTAL=$(grep 'total' $tmpfile | cut -f1 -d' ')
SUPPLEMENTARY=$(grep 'supplementary' $tmpfile | cut -f1 -d' ')
DUPLICATES=$(grep 'duplicates' $tmpfile | cut -f1 -d' ')
PAIRED=$(grep 'paired in sequencing' $tmpfile | cut -f1 -d' ')
prefix=$(basename $inbam)
echo -e "${inbam},${prefix},${TOTAL},${SUPPLEMENTARY},${DUPLICATES},${PAIRED}" >> $OUT
rm -rf $tmpfile
