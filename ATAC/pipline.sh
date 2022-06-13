##############primer options
Samples=(P53Ls-Rep1 P53Ls-Rep2 TeLs-Rep1 TeLs-Rep2 P53Ts-Rep1 P53Ts-Rep2 WTXT-Rep1 WTXT-Rep2)
SArrays={P53Ls-Rep1,P53Ls-Rep2,TeLs-Rep1,TeLs-Rep2,P53Ts-Rep1,P53Ts-Rep2,WTXT-Rep1,WTXT-Rep2}
Prefix=Hybri-ATAC
DataDir=01.data
MappingDir=02.mapping
ResultDir=03.result
StatDir=0a.stat
gtfFile=~/data/ref_data/annotation/XenTr_10.0_XenLa_10.1_GCF.gff3
XtrGtfFile="~/data/ref_data/annotation/XENTR_10.0_NCBI.clean.gtf"
XlaGtfFile="~/data/ref_data/annotation/XENLA_10.1_GCF.gtf"

mkdir ../$DataDir ../$MappingDir ../$ResultDir ../$StatDir


##========================data quality control and trim========================
for i in ${Samples[*]}
do
	fastqc -t 2 -o ../$DataDir/$i ../$DataDir/$i/${i}_L002_R1_001.fastq.gz ../$DataDir/$i/${i}_L002_R2_001.fastq.gz
	fastp -i ../$DataDir/$i/${i}_L002_R1_001.fastq.gz -I ../$DataDir/$i/${i}_L002_R2_001.fastq.gz -o ../$DataDir/$i/${i}_trim_R1.fq.gz -O ../$DataDir/$i/${i}_trim_R2.fq.gz -f 15 -F 15 --html ../$DataDir/$i/$i.trim.fastp.html --json ../$DataDir/$i/$i.trim.fastp.json -w 8 -c 
done
##=============================================================================
##=======mapping to fusion genome and extract data to different species========
for i in ${Samples[*]}
do
    Align.sh -I ~/data/ref_data/index/XenTr_v10.0_XenLa_v10.1_clean.bwa_index -1 ../$DataDir/$i/${i}_trim_R1.fq.gz -2 ../$DataDir/$i/${i}_trim_R2.fq.gz -p $i -o ../$MappingDir/$i -t 20
	#########extract data aligned to Xenopus tropicalis
	samtools view -b -@ 4 -o ../$MappingDir/$i/$i.sort.rmdup.filter.Xtr.bam ../$MappingDir/$i/$i.sort.rmdup.filter.bam Chr{1..10}
	samtools index ../$MappingDir/$i/$i.sort.rmdup.filter.Xtr.bam
	bamCoverage -b ../$MappingDir/$i/$i.sort.rmdup.filter.Xtr.bam -o ../$MappingDir/$i/$i.split.Xtr.bigwig --normalizeUsing RPKM --ignoreDuplicates -p 5 --binSize 10
	bigWigToBedGraph ../$MappingDir/$i/$i.split.Xtr.bigwig ../$MappingDir/$i/$i.split.Xtr.bedgraph
	#########extract data aligned to Xenopus laevis
	samtools view -b -@ 4 -o ../$MappingDir/$i/$i.sort.rmdup.filter.Xla.bam ../$MappingDir/$i/$i.sort.rmdup.filter.bam Chr{{1..8},9_10}S Chr{{1..8},9_10}L
	samtools index ../$MappingDir/$i/$i.sort.rmdup.filter.Xla.bam
	bamCoverage -b ../$MappingDir/$i/$i.sort.rmdup.filter.Xla.bam -o ../$MappingDir/$i/$i.split.Xla.bigwig --normalizeUsing RPKM --ignoreDuplicates -p 5 --binSize 10
	bigWigToBedGraph ../$MappingDir/$i/$i.split.Xla.bigwig ../$MappingDir/$i/$i.split.Xla.bedgraph
done
##=============================================================================
##============================statistic of alignment result====================
####insert size
bamPEFragmentSize -hist ../$StatDir/$Prefix.fragmentSize.pdf --plotFileFormat pdf -T "Fragment size" --maxFragmentLength 1000 -b ../$MappingDir/*/*.sort.rmdup.filter.bam --table ../$StatDir/$Prefix.fragmentSize.table --outRawFragmentLengths ../$StatDir/$Prefix.fragmentSize.data -p 24
###correlation
multiBamSummary bins -b ../$MappingDir/*/*.sort.rmdup.filter.bam -o ../$StatDir/$Prefix.multiBamSummary.npz -p 10
plotCorrelation -in ../$StatDir/$Prefix.multiBamSummary.npz -c pearson -p heatmap -o ../$StatDir/$Prefix.cor.pearson.pdf --skipZeros --outFileCorMatrix ../$StatDir/$Prefix.cor.pearson.matrix --plotFileFormat pdf
###region plot
echo computeMatrix scale-regions -R $XtrGtfFile $XlaGtfFile -S ../$MappingDir/*/$SArrays.RPKM.bigwig --skipZeros -b 1000 -a 1000 -m 5000 -p 10 -o ../$StatDir/$Prefix.deeptools.matrix.gz --samplesLabel ${Samples[*]}
plotProfile -m ../$StatDir/$Prefix.deeptools.matrix.gz -o ../$StatDir/$Prefix.deeptools.profile.pdf --plotFileFormat pdf --regionsLabel Xtr Xla
####Fingerprint
plotFingerprint -b ../$MappingDir/$SArrays/*.sort.rmdup.filter.Xtr.bam --smartLabels -bs 100 --skipZeros -T "Fingerprints of different samples in Xtr" -o ../$StatDir/$Prefix.deeptools.fingerprint.Xtr.pdf --outRawCounts ../$StatDir/$Prefix.deeptools.fingerprint.Xtr.tab -p 10 --labels $SArrays
plotFingerprint -b ../$MappingDir/$SArrays/*.sort.rmdup.filter.Xla.bam --smartLabels -bs 100 --skipZeros -T "Fingerprints of different samples in Xla" -o ../$StatDir/$Prefix.deeptools.fingerprint.Xla.pdf --outRawCounts ../$StatDir/$Prefix.deeptools.fingerprint.Xla.tab -p 10 --labels $SArrays
##=============================================================================
##====================================call peak================================
for i in ${Samples[*]}
do
	#######Xenopus tropicalis
	macs2 callpeak -t ../$MappingDir/$i/$i.sort.rmdup.filter.Xtr.bam -f BAM -g 1.4e9 -n $i.split.Xtr.macs2 --outdir ../$ResultDir/${i} -q 0.1 --nomodel --shift -70 --extsize 140 --keep-dup all
	sort -k8,8nr ../$ResultDir/${i}/$i.split.Xtr.macs2_peaks.narrowPeak > ../$ResultDir/${i}/$i.split.Xtr.macs2_peaks.sort.narrowPeak
	######Xenopus laevis
	macs2 callpeak -t ../$MappingDir/$i/$i.sort.rmdup.filter.Xla.bam -f BAM -g 2.6e9 -n $i.split.Xla.macs2 --outdir ../$ResultDir/${i} -q 0.1 --nomodel --shift -70 --extsize 140 --keep-dup all
	sort -k8,8nr ../$ResultDir/${i}/$i.split.Xla.macs2_peaks.narrowPeak > ../$ResultDir/${i}/$i.split.Xla.macs2_peaks.sort.narrowPeak	
done
#################  call merged peak
for i in {P53Ls-Rep,TeLs-Rep,P53Ts-Rep,WTXT-Rep}
do
	############Xenopus tropicalis
	macs2 callpeak -t ../$MappingDir/${i}1/${i}1.sort.rmdup.filter.Xtr.bam ../$MappingDir/${i}2/${i}2.sort.rmdup.filter.Xtr.bam -f BAM -g 1.4e9 -n $i.split.Xtr.merged.macs2 --outdir ../$ResultDir/merged -q 0.1 --nomodel --shift -70 --extsize 140 --keep-dup all
	sort -k8,8nr ../$ResultDir/merged/$i.split.Xtr.merged.macs2_peaks.narrowPeak > ../$ResultDir/merged/$i.split.Xtr.merged.macs2_peaks.sort.narrowPeak
	############Xenopus laevis
	macs2 callpeak -t ../$MappingDir/${i}1/${i}1.sort.rmdup.filter.Xla.bam ../$MappingDir/${i}2/${i}2.sort.rmdup.filter.Xla.bam -f BAM -g 2.6e9 -n $i.split.Xla.merged.macs2 --outdir ../$ResultDir/merged -q 0.1 --nomodel --shift -70 --extsize 140 --keep-dup all
	sort -k8,8nr ../$ResultDir/merged/$i.split.Xla.merged.macs2_peaks.narrowPeak > ../$ResultDir/merged/$i.split.Xla.merged.macs2_peaks.sort.narrowPeak
done
##=============================================================================
##=============bedtools annotation and get high confident peaks================
############Xenopus tropicalis
for i in {P53Ls,P53Ls-Rep,TeLs,TeLs-Rep,P53Ts,P53Ts-Rep,WTXT,WTXT-Rep}
do
	#####get information about merged peaks overlap with peaks from individual sample
	bedtools annotate -i ../$ResultDir/merged/$i.split.Xtr.merged.macs2_peaks.sort.narrowPeak -files ../$ResultDir/${i}1/${i}1.split.Xtr.macs2_peaks.sort.narrowPeak ../$ResultDir/${i}2/${i}2.split.Xtr.macs2_peaks.sort.narrowPeak -names ${i}1 ${i}2 > ../$ResultDir/merged/$i.split.Xtr.merged.raw_consensus.peaks.bed
	#####remain high confident peaks
	awk '{if($11>0.5 && $12>0.5){printf("%s\t%d\t%d\t%s\n",$1,$2,$3,$4)}}' ../$ResultDir/merged/$i.split.Xtr.merged.raw_consensus.peaks.bed | sed "s/$i.split.Xtr.merged.macs2_//" > ../$ResultDir/merged/$i.split.Xtr.merged.consensus.peaks.bed
	####motif enrichment for high confident peaks
	findMotifsGenome.pl ../$ResultDir/merged/$i.split.Xtr.merged.consensus.peaks.bed /work/bio-jiangh/data/ref_data/fasta/XENTR_10.0_genome.clean.fasta ../04.anno/Consensus/$i.xtr -mask -p 16
done
############Xenopus laevis
for i in {P53Ls,P53Ls-Rep,TeLs,TeLs-Rep,P53Ts,P53Ts-Rep,WTXT,WTXT-Rep}
do
	#####get information about merged peaks overlap with peaks from individual sample
	bedtools annotate -i ../$ResultDir/merged/$i.split.Xla.merged.macs2_peaks.sort.narrowPeak -files ../$ResultDir/${i}1/${i}1.split.Xla.macs2_peaks.sort.narrowPeak ../$ResultDir/${i}2/${i}2.split.Xla.macs2_peaks.sort.narrowPeak -names ${i}1 ${i}2 > ../$ResultDir/merged/$i.split.Xla.merged.raw_consensus.peaks.bed
	#####remain high confident peaks
	awk '{if($11>0.5 && $12>0.5){printf("%s\t%d\t%d\t%s\n",$1,$2,$3,$4)}}' ../$ResultDir/merged/$i.split.Xla.merged.raw_consensus.peaks.bed | sed "s/$i.split.Xla.merged.macs2_//" > ../$ResultDir/merged/$i.split.Xla.merged.consensus.peaks.bed
	####motif enrichment of high confident peaks
	findMotifsGenome.pl ../$ResultDir/merged/$i.split.Xla.merged.consensus.peaks.bed /work/bio-jiangh/data/ref_data/fasta/XENLA_10.1_genome.clean.fa ../04.anno/Consensus/$i.xla -mask -p 16
done
##=============================================================================
##========================BindDiff find different peaks========================
##use R script "Diif_peak_analysis.R"

##=============================================================================
##=====================motif enrichment of different peaks=====================
for i in {P53Ls-WT,P53Ls-P53Ts,TeLs-WT,P53Ts-WT,TeLs-P53Ls}
do
	#######accessibility up peaks
	grep 'up' ../04.anno/$i.peaks.DiffBind.tsv | awk '{printf("%s\t%d\t%d\n",$1,$2,$3)}' > ../04.anno/$i.peaks.DiffBind.up.bed
	findMotifsGenome.pl ../04.anno/$i.peaks.DiffBind.up.bed /work/bio-jiangh/data/ref_data/fasta/XENTR_10.0_genome.clean.fasta ../04.anno/$i.up -mask -p 16
	#######accessibility down peaks
	grep 'down' ../04.anno/$i.peaks.DiffBind.tsv | awk '{printf("%s\t%d\t%d\n",$1,$2,$3)}' > ../04.anno/$i.peaks.DiffBind.down.bed
	findMotifsGenome.pl ../04.anno/$i.peaks.DiffBind.down.bed /work/bio-jiangh/data/ref_data/fasta/XENTR_10.0_genome.clean.fasta ../04.anno/$i.down -mask -p 16
done
##=============================================================================
