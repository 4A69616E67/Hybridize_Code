##Primer options
######this code just a example for RNA-Seq analysis
Samples=(TeLs-s8-0 TeLs-s8-1 TeLs-s9-0 TeLs-s9-1 TeLs-s10-6_5h-0 TeLs-s10-6_5h-1 TeLs-s11-7_5h-0 TeLs-s11-7_5h-1 TeTs-s8-0 TeTs-s8-1 TeTs-s9-0 TeTs-s9-1 TeTs-s10-0 TeTs-s10-1 TeTs-s11-0 TeTs-s11-1)
$Sarrays={TeLs-s8-0,TeLs-s8-1,TeLs-s9-0,TeLs-s9-1,TeLs-s10-6_5h-0,TeLs-s10-6_5h-1,TeLs-s11-7_5h-0,TeLs-s11-7_5h-1,TeTs-s8-0,TeTs-s8-1,TeTs-s9-0,TeTs-s9-1,TeTs-s10-0,TeTs-s10-1,TeTs-s11-0,TeTs-s11-1} 
Prefix=SZY-Hybri-RNA
DataDir=01.data
MappingDir=02.mapping
CountDir=03.count
AssemblyDir=04.assembly
StatDir=0a.stat
gtfFile=~/data/ref_data/annotation/XenTr_10.0_XenLa_10.1_GCF.gff3

mkdir ../$DataDir ../$MappingDir ../ $CountDir ../$AssemblyDir ../$StatDir
##========================data quality control and trim========================
for i in ${Samples[*]}
do
    fastp -i $i/${i}_R1.fq.gz -I $i/${i}_R2.fq.gz -o $i/${i}_trim_R1.fq.gz -O $i/${i}_trim_R2.fq.gz -f 12 -F 12 --html $i/$i.trim.fastp.html --json $i/${i}.fastp.json -w 8 -c
done
##=============================================================================
##=======mapping to fusion genome and extract data to different species========
for i in ${Samples[*]}
do
    RNA-Align.sh -I ~/ref_data/index/XenTr_v10.0_XenLa_v10.1_clean.hisat2_index -1 ../$DataDir/$i/${i}_trim_R1.fq.gz -2 ../$DataDir/$i/${i}_trim_R2.fq.gz -p $i -o ../$MappingDir/$i -t 16
	bamCoverage -b ../$MappingDir/$i/$i.filter.sort.rmdup.bam -o ../$MappingDir/$i/$i.RPKM.bigwig --binSize 10 --normalizeUsing RPKM -p 4 --ignoreDuplicates
	##extract data in Xenopus tropicalis
	samtools view -b -@ 4 -o ../$MappingDir/$i/$i.filter.sort.rmdup.Xtr.bam ../$MappingDir/$i/$i.filter.sort.rmdup.bam Chr{1..10}
    samtools index ../$MappingDir/$i/$i.filter.sort.rmdup.Xtr.bam
    bamCoverage -b ../$MappingDir/$i/$i.filter.sort.rmdup.Xtr.bam -o ../$MappingDir/$i/$i.split.Xtr.bigwig --binSize 10 --normalizeUsing RPKM -p 4 --ignoreDuplicates
	##extract data in Xenopus laevis
    samtools view -b -@ 4 -o ../$MappingDir/$i/$i.filter.sort.rmdup.Xla.bam ../$MappingDir/$i/$i.filter.sort.rmdup.bam Chr{{1..8},9_10}{S,L}
    samtools index ../$MappingDir/$i/$i.filter.sort.rmdup.Xla.bam
    bamCoverage -b ../$MappingDir/$i/$i.filter.sort.rmdup.Xla.bam -o ../$MappingDir/$i/$i.split.Xla.bigwig --binSize 10 --normalizeUsing RPKM -p 4 --ignoreDuplicates
done

#######gene coverage in Xtr#########
computeMatrix scale-regions -S ../$MappingDir/$Sarrays/*.RPKM.bigwig -R ~/data/ref_data/annotation/XENTR_10.0_NCBI.gtf --metagene -o Hybridize.deeptools.xtr.gene_coverage.matrix.gz -p 24
plotProfile -m Hybridize.deeptools.xtr.gene_coverage.matrix.gz -o Hybridize.deeptools.xtr.gene_coverage.pdf --plotFileFormat pdf --outFileNameData Hybridize.deeptools.xtr.gene_coverage.tab
#######gene coverage in Xla#########
computeMatrix scale-regions -S ../$MappingDir/*/*.Xla.bigwig -R ~/data/ref_data/annotation/XENLA_10.1_GCF.gtf --metagene -o Hybridize.deeptools.xla.gene_coverage.matrix.gz -p 24
plotProfile -m Hybridize.deeptools.xla.gene_coverage.matrix.gz -o Hybridize.deeptools.xla.gene_coverage.pdf --plotFileFormat pdf --outFileNameData Hybridize.deeptools.xla.gene_coverage.tab
##########samples correlation############
multiBigwigSummary bins -b ../$MappingDir/*.bigwig -o ../result/Hybri.deeptools.bw.npz --labels $Samples -p 10
plotCorrelation -in Hybri.deeptools.bw.npz -c pearson -p heatmap -o Hybri.cor.pdf --skipZeros --removeOutliers --outFileCorMatrix Transparent.cor.tab --plotFileFormat pdf

##########Xtr gene Count 
featureCounts -T 10 -p -g gene -a ~/data/ref_data/annotation/XENTR_10.0_GCF.gff3 -o ../$CountDir/Hybridize.featureCount.Xtr ../$MappingDir/*/*.filter.sort.rmdup.bam

###Xla gene Count
featureCounts -T 10 -p -g gene -a ~/data/ref_data/annotation/XENLA_10.1_GCF.gff3 -o ../$CountDir/Hybridize.featureCount.Xla ../$MappingDir/*/*.filter.sort.rmdup.bam

#####file format
perl featureCount2matrix.pl Hybridize.featureCount.Xtr > Hybridize.featureCount.Xtr.matrix
perl featureCount2matrix.pl Hybridize.featureCount.Xla > Hybridize.featureCount.Xla.matrix

#######################
####merge gene count in Xtr and Xla using R script
####Different gene expression processed using R script