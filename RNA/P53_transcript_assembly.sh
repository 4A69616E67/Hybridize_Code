###################################tp53 transcript assembly using stringtie
Samples=(TeLs-s8-0 TeLs-s8-1 TeLs-s9-0 TeLs-s9-1 TeLs-s10-6_5h-0 TeLs-s10-6_5h-1 TeLs-s11-7_5h-0 TeLs-s11-7_5h-1 TeTs-s8-0 TeTs-s8-1 TeTs-s9-0 TeTs-s9-1 TeTs-s10-0 TeTs-s10-1 TeTs-s11-0 TeTs-s11-1)
for i in ${Samples[*]}
do
    samtools view -b ../02.mapping/$i/$i.filter.sort.rmdup.Xtr.bam Chr3:153630000-153642000 > ../02.mapping/$i/$i.filter.sort.rmdup.Xtr.tp53.bam
	stringtie ../02.mapping/$i/$i.filter.sort.rmdup.Xtr.tp53.bam -o ../02.mapping/$i/$i.stringtie_assembly.tp53.gtf -p 10 
done

Samples=(TeLs-s8-0 TeLs-s8-1 TeLs-s9-0 TeLs-s9-1 TeLs-s10-6_5h-0 TeLs-s10-6_5h-1 TeLs-s11-7_5h-0 TeLs-s11-7_5h-1 TeTs-s8-0 TeTs-s8-1 TeTs-s9-0 TeTs-s9-1 TeTs-s10-0 TeTs-s10-1 TeTs-s11-0 TeTs-s11-1)
for i in ${Samples[*]}
do
    samtools view -b ../02.mapping/$i/$i.filter.sort.rmdup.Xla.bam Chr3S:131,802,000-131,824,000 > ../02.mapping/$i/$i.filter.sort.rmdup.Xla.tp53.bam
	stringtie ../02.mapping/$i/$i.filter.sort.rmdup.Xla.tp53.bam -o ../02.mapping/$i/$i.stringtie_assembly.Xla.tp53.gtf -p 10 
done

