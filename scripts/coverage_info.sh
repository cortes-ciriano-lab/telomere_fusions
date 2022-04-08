# /usr/local/bin/bash
inbam=$1
OUT=$2 
uuid=$(uuidgen)
tmpfile=tmpfile_${uuid}
touch ${tmpfile}_${uuid}

samtools flagstat $inbam > ${tmpfile}_${uuid}

echo -e "Path_file,File,Total_reads,Supplementary_reads,Duplicate_reads,Paired_reads" > $OUT
TOTAL=$(grep 'total' ${tmpfile}_${uuid} | cut -f1 -d' ')
SUPPLEMENTARY=$(grep 'supplementary' ${tmpfile}_${uuid} | cut -f1 -d' ')
DUPLICATES=$(grep 'duplicates' ${tmpfile}_${uuid} | cut -f1 -d' ')
PAIRED=$(grep 'paired in sequencing' ${tmpfile}_${uuid} | cut -f1 -d' ')
prefix=$(basename $inbam)
echo -e "${inbam},${prefix},${TOTAL},${SUPPLEMENTARY},${DUPLICATES},${PAIRED}" >> $OUT
rm -rf ${tmpfile}_${uuid}
