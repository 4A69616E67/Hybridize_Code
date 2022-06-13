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

RAWBam="${OUTDIR}/${PREFIX}.bam"
FilterSortBam="${OUTDIR}/${PREFIX}.filter.sort.bam"
UniqBam="${OUTDIR}/${PREFIX}.filter.sort.rmdup.bam"
BigWig="${OUTDIR}/${PREFIX}.10bp.RPKM.bigwig"
Command=""

if [ ! -d $OUTDIR ]
then
    mkdir $OUTDIR
fi

if [ -z "$FQ2" ]
then
    hisat2 -p $THREAD -x $INDEX -U $FQ1 | samtools view -@ $THREAD -Sbh -o $RAWBam
#    bwa mem -t $THREAD $INDEX $FQ1 | samtools view -@ $THREAD -Sbh -o $RAWBam
#    samtools flagstat $RAWBam > ${OUTDIR}/${PREFIX}.flagstat.txt
#    samtools view -F 260 -q 20 -@ $THREAD -bh $RAWBam | samtools sort -@ $THREAD -o $FilterSortBam
else
    hisat2 -p $THREAD -x $INDEX -1 $FQ1 -2 $FQ2 | samtools view -@ $THREAD -Sbh -o $RAWBam
#    bwa mem -t $THREAD $INDEX $FQ1 $FQ2 | samtools view -@ $THREAD -Sbh -o $RAWBam
#    samtools flagstat $RAWBam > ${OUTDIR}/${PREFIX}.flagstat.txt
#    samtools view -F 268 -q 20 -@ $THREAD -bh $RAWBam | samtools sort -@ $THREAD -o $FilterSortBam
fi
samtools flagstat $RAWBam > ${OUTDIR}/${PREFIX}.flagstat.txt
samtools view -F 260 -q 20 -@ $THREAD -bh $RAWBam | samtools sort -@ $THREAD -o $FilterSortBam
samtools index -@ $THREAD $FilterSortBam
sambamba markdup -r -t $THREAD $FilterSortBam $UniqBam
bamCoverage -b $UniqBam -o $BigWig --binSize 10 --normalizeUsing RPKM -p $THREAD
