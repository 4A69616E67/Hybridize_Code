########primer options
Samples=(s9ctcf-1 s9ctcf-2 s9Inp s9k27ac-1 s9k27ac-2 s9k27me3-1 s9k27me3-2 s9k4me1-1 s9k4me1-2 s9k4me3-1 s9k4me3-2 s9p53-1 s9p53-2 s9pol-1 s9pol-2 s9rad21-1 s9rad21-2 telsctcf-1 telsctcf-2 telsInp telsk27ac-1 telsk27ac-2 telsk27me3-1 telsk27me3-2 telsk4me1-1 telsk4me1-2 telsk4me3-1 telsk4me3-2 telsp53-1 telsp53-2 telspol-1 telspol-2 telsrad21-1 telsrad21-2)
SArrays={s9ctcf-1,s9ctcf-2,s9Inp,s9k27ac-1,s9k27ac-2,s9k27me3-1,s9k27me3-2,s9k4me1-1,s9k4me1-2,s9k4me3-1,s9k4me3-2,s9p53-1,s9p53-2,s9pol-1,s9pol-2,s9rad21-1,s9rad21-2,telsctcf-1,telsctcf-2,telsInp,telsk27ac-1,telsk27ac-2,telsk27me3-1,telsk27me3-2,telsk4me1-1,telsk4me1-2,telsk4me3-1,telsk4me3-2,telsp53-1,telsp53-2,telspol-1,telspol-2,telsrad21-1,telsrad21-2}
Prefix=SZY-Hybridize-ChIP
DataDir=01.data
MappingDir=02.mapping
ResultDir=03.result
StatDir=0a.stat
gtfFile=~/data/ref_data/annotation/XenTr_10.0_XenLa_10.1_GCF.gff3
XtrGtfFile="~/data/ref_data/annotation/XENTR_10.0_NCBI.clean.gtf"
XlaGtfFile="~/data/ref_data/annotation/XENLA_10.1_GCF.gtf"

mkdir ../$DataDir ../$MappingDir ../$ResultDir ../$StatDir

##Note: remove prefix and suffix of raw data file name before processing
##========================data quality control and trim========================
for i in ${Samples[*]}
do
	fastqc -t 2 -o ../$DataDir/$i ../$DataDir/$i/${i}_R1_001.fastq.gz ../$DataDir/$i/${i}_R2_001.fastq.gz
	fastp -i ../$DataDir/$i/${i}_R1_001.fastq.gz -I ../$DataDir/$i/${i}_R2_001.fastq.gz -o ../$DataDir/$i/${i}_trim_R1.fq.gz -O ../$DataDir/$i/${i}_trim_R2.fq.gz -f 10 -F 10 --html ../$DataDir/$i/$i.trim.fastp.html --json ../$DataDir/$i/$i.trim.fastp.json -w 16 -c
done
##=============================================================================
##=======mapping to fusion genome and extract data to different species========
for i in ${Samples[*]}
do
	Align.sh -I ~/data/ref_data/index/XenTr_v10.0_XenLa_v10.1_clean.bwa_index -1 ../$DataDir/$i/${i}_trim_R1.fq.gz -2 ../$DataDir/$i/${i}_trim_R2.fq.gz -p $i -o ../$MappingDir/$i -t 16
	#########extract data aligned to Xenopus tropicalis
	samtools view -b -@ 4 -o ../$MappingDir/$i/$i.sort.rmdup.filter.Xtr.bam ../$MappingDir/$i/$i.sort.rmdup.filter.bam Chr{1..10}
	samtools index ../$MappingDir/$i/$i.sort.rmdup.filter.Xtr.bam
	bamCoverage -b ../$MappingDir/$i/$i.sort.rmdup.filter.Xtr.bam -o ../$MappingDir/$i/$i.RPKM.Xtr.bigwig --binSize 10 --normalizeUsing RPKM -p 8 --ignoreDuplicates
	#########extract data aligned to Xenopus laevis
	samtools view -b -@ 4 -o ../$MappingDir/$i/$i.sort.rmdup.filter.Xla.bam ../$MappingDir/$i/$i.sort.rmdup.filter.bam Chr{{1..8},9_10}S Chr{{1..8},9_10}L
	samtools index ../$MappingDir/$i/$i.sort.rmdup.filter.Xla.bam
	bamCoverage -b ../$MappingDir/$i/$i.sort.rmdup.filter.Xla.bam -o ../$MappingDir/$i/$i.RPKM.Xla.bigwig --binSize 10 --normalizeUsing RPKM -p 8 --ignoreDuplicates
done
##=============================================================================
##============================statistic of alignment result====================
##correlation
multiBamSummary bins -b ../$MappingDir/$SArrays/*.sort.rmdup.filter.bam -o ../$StatDir/$Prefix.multiBamSummary.npz -p 10 --labels ${Samples[*]} ; plotCorrelation -in ../$StatDir/$Prefix.multiBamSummary.npz -c pearson -p heatmap -o ../$StatDir/$Prefix.cor.pearson.pdf --skipZeros --outFileCorMatrix ../$StatDir/$Prefix.cor.pearson.matrix --plotFileFormat pdf
##region plot Xtr
computeMatrix scale-regions -R ~/data/ref_data/annotation/XENTR_10.0_NCBI.gtf -S ../$MappingDir/$SArrays/*.RPKM.Xtr.bigwig --skipZeros -b 1000 -a 1000 -m 5000 -p 12 --missingDataAsZero -o ../$StatDir/$Prefix.deeptools.Xtr.matrix.gz --samplesLabel ${Samples[*]}
plotProfile -m ../$StatDir/$Prefix.deeptools.Xtr.matrix.gz -o ../$StatDir/$Prefix.deeptools.Xtr.profile.pdf --plotFileFormat pdf
plotHeatmap -m ../$StatDir/$Prefix.deeptools.Xtr.matrix.gz -out ../$StatDir/$Prefix.deeptools.Xtr.profile.heatmap.pdf --colorList white,red --plotFileFormat pdf
##fingerprint Xtr
plotFingerprint -b ../$MappingDir/$SArrays/*.sort.rmdup.filter.Xtr.bam --labels ${Samples[*]} -bs 100 --skipZeros -T "Fingerprints of different samples" -o ../$StatDir/$Prefix.deeptools.Xtr.fingerprint.pdf --outRawCounts ../$StatDir/$Prefix.deeptools.Xtr.fingerprint.tab -p 10
##region plot Xla
computeMatrix scale-regions -R ~/data/ref_data/annotation/XENTR_10.0_NCBI.gtf -S ../$MappingDir/$SArrays/*.RPKM.Xla.bigwig --skipZeros -b 1000 -a 1000 -m 5000 -p 12 --missingDataAsZero -o ../$StatDir/$Prefix.deeptools.Xla.matrix.gz --samplesLabel ${Samples[*]}
plotProfile -m ../$StatDir/$Prefix.deeptools.Xla.matrix.gz -o ../$StatDir/$Prefix.deeptools.Xla.profile.pdf --plotFileFormat pdf
plotHeatmap -m ../$StatDir/$Prefix.deeptools.Xla.matrix.gz -out ../$StatDir/$Prefix.deeptools.Xla.profile.heatmap.pdf --colorList white,red --plotFileFormat pdf
##fingerprint Xla
plotFingerprint -b ../$MappingDir/$SArrays/*.sort.rmdup.filter.Xla.bam --labels ${Samples[*]} -bs 100 --skipZeros -T "Fingerprints of different samples" -o ../$StatDir/$Prefix.deeptools.Xla.fingerprint.pdf --outRawCounts ../$StatDir/$Prefix.deeptools.Xla.fingerprint.tab -p 10
##=============================================================================
##====================================call peak================================
##Xenopus tropicalis wild type
for i in {s9ctcf-{1,2},s9k27ac-{1,2},s9k27me3-{1,2},s9k4me1-{1,2},s9k4me3-{1,2},s9p53-{1,2},s9pol-{1,2},s9rad21-{1,2}}
do
	macs2 callpeak -t ../$MappingDir/${i}/$i.sort.rmdup.filter.Xtr.bam -c ../$MappingDir/s9Inp/s9Inp.sort.rmdup.filter.Xtr.bam -f BAMPE -g 1.4e9 -n $i.macs2 --outdir ../$ResultDir/$i -q 0.05
done
##Xenopus tropicalis hybridize
for i in {telsctcf-{1,2},telsk27ac-{1,2},telsk27me3-{1,2},telsk4me1-{1,2},telsk4me3-{1,2},telsp53-{1,2},telspol-{1,2},telsrad21-{1,2}}
do
	macs2 callpeak -t ../$MappingDir/${i}/$i.sort.rmdup.filter.Xtr.bam -c ../$MappingDir/telsInp/telsInp.sort.rmdup.filter.Xtr.bam -f BAMPE -g 1.4e9 -n $i.macs2 --outdir ../$ResultDir/$i -q 0.05
done
###############call peak use all samples
##Xenopus tropicalis wild type
for i in {s9ctcf,s9k27ac,s9k27me3,s9k4me1,s9k4me3,s9p53,s9pol,s9rad21}
do
	macs2 callpeak -t ../$MappingDir/{${i}-1,${i}-2}/*.sort.rmdup.filter.Xtr.bam -c ../$MappingDir/s9Inp/s9Inp.sort.rmdup.filter.Xtr.bam -f BAMPE -g 1.4e9 -n $i.macs2 --outdir ../$ResultDir/$i -q 0.05
done
##Xenopus tropicalis hybridize
for i in {telsctcf,telsk27ac,telsk27me3,telsk4me1,telsk4me3,telsp53,telspol,telsrad21}
do
	macs2 callpeak -t ../$MappingDir/{${i}-1,${i}-2}/*.sort.rmdup.filter.Xtr.bam -c ../$MappingDir/telsInp/telsInp.sort.rmdup.filter.Xtr.bam -f BAMPE -g 1.4e9 -n $i.macs2 --outdir ../$ResultDir/$i -q 0.05
done
##=============================================================================
##=========bedtools annotation and get high confident peaks (Xtr)==============
for i in {s9ctcf,s9k27ac,s9k27me3,s9k4me1,s9k4me3,s9p53,s9pol,s9rad21,telsctcf,telsk27ac,telsk27me3,telsk4me1,telsk4me3,telsp53,telspol,telsrad21}
do
	#####get information about merged peaks overlap with peaks from individual sample
	bedtools annotate -i ../$ResultDir/$i/$i.macs2_peaks.narrowPeak -files ../$ResultDir/${i}-1/${i}-1.macs2_peaks.narrowPeak ../$ResultDir/${i}-2/${i}-2.macs2_peaks.narrowPeak -names ${i}-1 ${i}-2 > ../$ResultDir/$i/$i.raw_consensus.peaks.bed
	#####remain high confident peaks
	awk '{if(\$11>0.5 && \$12>0.5){printf(\"%s\t%d\t%d\t%s\n\",\$1,\$2,\$3,\$4)}}' ../$ResultDir/$i/$i.raw_consensus.peaks.bed | sed "s/$i.macs2_//" > ../$ResultDir/$i/$i.consensus.peaks.bed
	##motif enrichment
	findMotifsGenome.pl ../$ResultDir/$i/$i.consensus.peaks.bed /work/bio-jiangh/data/ref_data/fasta/XENTR_10.0_genome.clean.fasta ../$ResultDir/$i/motif -mask -p 16
	##peak annotation
	annotatePeaks.pl ../$ResultDir/$i/$i.consensus.peaks.bed ~/data/ref_data/fasta/XENTR_10.0_genome.clean.fasta -m ../$ResultDir/$i/motif/homerResults/motif1.motif ../$ResultDir/$i/motif/homerResults/motif2.motif -gtf $XtrGtfFile -annStats ../$ResultDir/$i/motif/$i.motif.anno.stat > ../$ResultDir/$i/motif/$i.motif.anno.txt
done
##=============================================================================
##=============================================================================