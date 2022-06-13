UsageHelp="Usage: ${0} -I <gemone index> -1 <fsatq1> [-2 <fastq2>] [-p <out prefix>] [-o <out dir>] [-t <threads>]"
while getopts "I:p:1:2:t:o:" opt; do
    case $opt in 
        I) INDEX=$OPTARG ;;
	p) PREFIX=$OPTARG ;;
	1) FQ1=$OPTARG ;;
	2) FQ2=$OPTARG ;;
	t) THREAD=$OPTARG ;;
	o) OUTDIR=$OPTARG ;;
	[?]) echo "$UsageHelp" && exit 1 ;;
    esac
done

if [ -z "$INDEX" ] || [ -z "$FQ1" ] 
then
    echo "$UsageHelp"
    exit 1
fi

if [ -z "$PREFIX" ]
then
    PREFIX='out'
fi

if [ -z "$OUTDIR" ]
then
    OUTDIR='./'
fi

if [ -z "$THREAD" ]
then
    THREAD='1'
fi

RawSortBam="${OUTDIR}/${PREFIX}.sort.bam"
MarkDupBam="${OUTDIR}/${PREFIX}.sort.markdup.bam"
FilterSortBam="${OUTDIR}/${PREFIX}.sort.markdup.filter.bam"
RmDupBam="${OUTDIR}/${PREFIX}.sort.rmdup.filter.bam"
BigWig="${OUTDIR}/${PREFIX}.RPKM.bigwig"
Command=""
if [ -z "$FQ2" ]
then
    bwa mem -t $THREAD $INDEX $FQ1 | samtools view -@ $THREAD -Sbh | samtools sort -@ $THREAD -o $RawSortBam
    sambamba markdup -t $THREAD $RawSortBam $MarkDupBam
    samtools flagstat $MarkDupBam > ${OUTDIR}/${PREFIX}.flagstat.txt
    samtools view -F 260 -q 20 -@ $THREAD -bh -o $FilterSortBam $MarkDupBam
else
    bwa mem -t $THREAD $INDEX $FQ1 $FQ2 | samtools view -@ $THREAD -Sbh | samtools sort -@ $THREAD -o $RawSortBam
    sambamba markdup -t $THREAD $RawSortBam $MarkDupBam
    samtools flagstat $MarkDupBam > ${OUTDIR}/${PREFIX}.flagstat.txt
    samtools view -F 268 -q 20 -@ $THREAD -bh -o $FilterSortBam $MarkDupBam
fi
samtools index -@ $THREAD $FilterSortBam
samtools rmdup -S $FilterSortBam $RmDupBam
samtools index -@ $THREAD $RmDupBam
bamCoverage -b $FilterSortBam -o $BigWig --binSize 10 --normalizeUsing RPKM -p $THREAD --ignoreDuplicates
